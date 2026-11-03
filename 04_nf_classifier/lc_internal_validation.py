"""Patient-level internal LC validation of the MSI-selected five-feature panel.

The script never selects MSI features from LC abundance. It reads the final
MSI panel produced by panel_selection/analyze.py, joins the reviewed MSI-to-LC
peak correspondence table, and stops if any selected MSI feature lacks an LC
peak assignment. The matched LC peaks are then evaluated by nested
leave-one-patient-out multinomial ridge regression.
"""

from pathlib import Path
import hashlib
import json
import os

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.optimize import minimize


HERE = Path(__file__).resolve().parent
PANEL_FILE = Path(os.environ.get(
    "MSI_PANEL_FILE", HERE / "panel_selection" / "final_selection.json"
))
MAPPING_FILE = Path(os.environ.get(
    "LC_PEAK_MAPPING_FILE", HERE / "lc_peak_correspondence.csv"
))
LC_WORKBOOK = Path(os.environ.get(
    "LC_INPUT_WORKBOOK", HERE / "LC_input_snapshot.xlsx"
))
LC_SHEET = os.environ.get("LC_INPUT_SHEET", "LC_processed_matrix")
OUT = Path(os.environ.get(
    "LC_VALIDATION_OUT", HERE / "lc_internal_validation_outputs"
))

CLASSES = ["NF-GH", "NF-ACTH", "NF-GnH"]
SAMPLES = [f"{group}-{i}" for group in CLASSES for i in range(1, 5)]
Y = np.repeat(np.arange(len(CLASSES)), 4)
RIDGE_GRID = [0.01, 0.1, 1.0, 10.0]


def read_sheet(path, sheet):
    workbook = load_workbook(path, read_only=True, data_only=True)
    if sheet not in workbook.sheetnames:
        workbook.close()
        raise ValueError(f"LC worksheet not found: {sheet}")
    rows = workbook[sheet].values
    header = next(rows)
    data = pd.DataFrame(rows, columns=header)
    workbook.close()
    return data


def standardize_train(x_train, x_test):
    center = x_train.mean(axis=0)
    scale = x_train.std(axis=0, ddof=0)
    scale[~np.isfinite(scale) | (scale < 1e-12)] = 1.0
    return (x_train - center) / scale, (x_test - center) / scale, center, scale


def softmax_reference(logits_nonref):
    logits = np.column_stack([logits_nonref, np.zeros(len(logits_nonref))])
    logits -= logits.max(axis=1, keepdims=True)
    prob = np.exp(logits)
    return prob / prob.sum(axis=1, keepdims=True)


def fit_ridge(x, y, ridge_lambda):
    n_features = x.shape[1]
    n_nonref = len(CLASSES) - 1

    def unpack(theta):
        coef = theta[:n_features * n_nonref].reshape(n_features, n_nonref)
        intercept = theta[n_features * n_nonref:]
        return coef, intercept

    def objective(theta):
        coef, intercept = unpack(theta)
        prob = softmax_reference(x @ coef + intercept)
        nll = -np.log(np.maximum(prob[np.arange(len(y)), y], 1e-15)).mean()
        penalty = 0.5 * ridge_lambda * np.sum(coef ** 2)
        return nll + penalty

    initial = np.zeros(n_features * n_nonref + n_nonref)
    result = minimize(objective, initial, method="L-BFGS-B")
    if not result.success:
        raise RuntimeError(f"Ridge optimization failed: {result.message}")
    return unpack(result.x)


def predict_ridge(x_train, y_train, x_test, ridge_lambda):
    z_train, z_test, center, scale = standardize_train(x_train, x_test)
    coef, intercept = fit_ridge(z_train, y_train, ridge_lambda)
    probability = softmax_reference(z_test @ coef + intercept)
    return probability, center, scale, coef, intercept


def metrics(probability, truth):
    prediction = probability.argmax(axis=1)
    recall = [np.mean(prediction[truth == k] == k) for k in range(len(CLASSES))]
    balanced_accuracy = float(np.mean(recall))
    logloss = float(-np.log(
        np.maximum(probability[np.arange(len(truth)), truth], 1e-15)
    ).mean())
    return balanced_accuracy, logloss


def choose_lambda(x, y, training_indices):
    rows = []
    for ridge_lambda in RIDGE_GRID:
        probabilities = []
        truths = []
        for held_out in training_indices:
            inner_train = training_indices[training_indices != held_out]
            probability, *_ = predict_ridge(
                x[inner_train], y[inner_train], x[[held_out]], ridge_lambda
            )
            probabilities.append(probability[0])
            truths.append(y[held_out])
        probability = np.asarray(probabilities)
        truth = np.asarray(truths)
        balanced_accuracy, logloss = metrics(probability, truth)
        rows.append({
            "ridge_lambda": ridge_lambda,
            "balanced_accuracy": balanced_accuracy,
            "logloss": logloss,
        })
    scores = pd.DataFrame(rows)
    best = scores.sort_values(
        ["balanced_accuracy", "logloss", "ridge_lambda"],
        ascending=[False, True, True],
    ).iloc[0]
    return float(best.ridge_lambda), scores


def main():
    for path in [PANEL_FILE, MAPPING_FILE, LC_WORKBOOK]:
        if not path.exists():
            raise FileNotFoundError(path)
    OUT.mkdir(parents=True, exist_ok=True)

    panel = json.loads(PANEL_FILE.read_text())
    if "5" not in panel or "features" not in panel["5"]:
        raise ValueError("final_selection.json does not contain a five-feature panel")
    features = list(panel["5"]["features"])
    if len(features) != 5 or len(set(features)) != 5:
        raise ValueError(f"Expected five unique MSI features, found: {features}")

    mapping = pd.read_csv(MAPPING_FILE)
    if mapping.feature_id.duplicated().any():
        raise ValueError("LC correspondence table contains duplicated feature_id values")
    selected = pd.DataFrame({"feature_id": features}).merge(
        mapping, on="feature_id", how="left", validate="one_to_one"
    )
    missing = selected.loc[selected.LC_ID.isna(), "feature_id"].tolist()
    if missing:
        raise RuntimeError(
            "LC validation cannot start. Add reviewed LC peak correspondences for: "
            + ", ".join(missing)
        )
    selected.to_csv(OUT / "selected_lc_peak_correspondence.csv", index=False)

    lc = read_sheet(LC_WORKBOOK, LC_SHEET)
    required_columns = {"Feature_ID", *SAMPLES}
    absent_columns = sorted(required_columns.difference(lc.columns))
    if absent_columns:
        raise ValueError("LC table is missing columns: " + ", ".join(absent_columns))
    lc = lc.set_index("Feature_ID")
    missing_ids = [x for x in selected.LC_ID if x not in lc.index]
    if missing_ids:
        raise RuntimeError("LC peak IDs not found in workbook: " + ", ".join(missing_ids))

    rows = lc.loc[selected.LC_ID, SAMPLES]
    if rows.index.duplicated().any() or len(rows) != 5:
        raise ValueError("Each selected LC peak must map to exactly one workbook row")
    raw = rows.to_numpy(dtype=float).T
    if not np.isfinite(raw).all() or (raw < 0).any():
        raise ValueError("LC peak areas must be finite and non-negative")
    x = np.log2(np.maximum(raw, 1.0))
    pd.DataFrame(x, index=SAMPLES, columns=features).to_csv(
        OUT / "lc_log2_inputs.csv"
    )

    outer_rows = []
    inner_rows = []
    all_indices = np.arange(len(SAMPLES))
    for held_out in all_indices:
        outer_train = all_indices[all_indices != held_out]
        ridge_lambda, inner_scores = choose_lambda(x, Y, outer_train)
        inner_scores.insert(0, "held_out_patient", SAMPLES[held_out])
        inner_rows.append(inner_scores)
        probability, *_ = predict_ridge(
            x[outer_train], Y[outer_train], x[[held_out]], ridge_lambda
        )
        prediction = int(probability[0].argmax())
        outer_rows.append({
            "sample_id": SAMPLES[held_out],
            "reference_class": CLASSES[Y[held_out]],
            "predicted_class": CLASSES[prediction],
            "correct": prediction == Y[held_out],
            "selected_ridge_lambda": ridge_lambda,
            **{f"probability_{c}": probability[0, k]
               for k, c in enumerate(CLASSES)},
        })

    predictions = pd.DataFrame(outer_rows)
    predictions.to_csv(OUT / "lc_outer_predictions.csv", index=False)
    pd.concat(inner_rows, ignore_index=True).to_csv(
        OUT / "lc_inner_scores.csv", index=False
    )
    probability_columns = [f"probability_{c}" for c in CLASSES]
    probability = predictions[probability_columns].to_numpy()
    balanced_accuracy, logloss = metrics(probability, Y)
    summary = pd.DataFrame([{
        "n_patients": len(SAMPLES),
        "n_correct": int(predictions.correct.sum()),
        "accuracy": float(predictions.correct.mean()),
        "balanced_accuracy": balanced_accuracy,
        "logloss": logloss,
        **{
            f"recall_{c}": float(predictions.loc[
                predictions.reference_class == c, "correct"
            ].mean())
            for c in CLASSES
        },
    }])
    summary.to_csv(OUT / "lc_validation_summary.csv", index=False)

    final_lambda, final_inner_scores = choose_lambda(x, Y, all_indices)
    final_inner_scores.to_csv(OUT / "final_lambda_inner_scores.csv", index=False)
    z, _, center, scale = standardize_train(x, x)
    coefficient, intercept = fit_ridge(z, Y, final_lambda)
    np.savez(
        OUT / "final_lc_ridge_model.npz",
        coefficient=coefficient,
        intercept=intercept,
        center=center,
        scale=scale,
        features=np.asarray(features),
        LC_IDs=selected.LC_ID.to_numpy(),
        classes=np.asarray(CLASSES),
        reference_class=CLASSES[-1],
        ridge_lambda=final_lambda,
    )

    design = {
        "MSI_panel_file": str(PANEL_FILE),
        "LC_mapping_file": str(MAPPING_FILE),
        "LC_workbook_sha256": hashlib.sha256(LC_WORKBOOK.read_bytes()).hexdigest(),
        "LC_sheet": LC_SHEET,
        "features": features,
        "LC_IDs": selected.LC_ID.tolist(),
        "input": "log2(max(LC peak area, 1)); train-fold standardization",
        "validation": "outer leave-one-patient-out; ridge lambda selected by inner leave-one-patient-out",
        "ridge_grid": RIDGE_GRID,
        "selection_rule": "maximize inner balanced accuracy, then minimize logloss",
        "scope": "exploratory internal validation; MSI coefficients are not transferred to LC",
    }
    (OUT / "lc_validation_design.json").write_text(json.dumps(design, indent=2))
    print(summary.to_string(index=False))
    print("LC validation output:", OUT)


if __name__ == "__main__":
    main()
