import argparse
import json
import pandas as pd

drug_map = {
    "amikacin": "AMI",
    "bedaquiline": "BDQ",
    "capreomycin": "CAP",
    "clofazimine": "CFZ",
    "delamanid": "DLM",
    "ethambutol": "EMB",
    "ethionamide": "ETH",
    "isoniazid": "INH",
    "kanamycin": "KAN",
    "levofloxacin": "LEV",
    "linezolid": "LZD",
    "moxifloxacin": "MXF",
    "pyrazinamide": "PZA",
    "rifampicin": "RIF",
    "streptomycin": "STM",
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--jsons",
        type=str,
        nargs="+",
        required=True,
    )
    args = parser.parse_args()
    effects_table = {
        "ENA_RUN_ACCESSION": [],
        "DRUG": [],
        "GENE": [],
        "MUTATION": [],
        "confidence": [],
        "PREDICTION": [],
    }
    predictions_table = {
        "ENA_RUN_ACCESSION": [],
        "DRUG": [],
        "PREDICTION": [],
    }

    values = "RFUS"
    items = {}
    for j in args.jsons:
        data = json.load(open(j, "r", encoding="utf-8"))
        sample_id = data.get("id")
        drug_table = data.get("drug_table", [])
        for entry in drug_table:
            drug = drug_map.get(entry.get("drug"))
            gene = entry.get("gene")
            mutation = entry.get("change")
            confidence = entry.get("confidence")
            if not drug:
                continue

            if confidence == "Uncertain significance":
                prediction = "U"
            elif confidence in ["Assoc w R", "Assoc w R - Interim"]:
                prediction = "R"
            elif confidence in ["", "Not assoc w R", "Not assoc w R - Interim"]:
                prediction = "S"
            else:
                raise ValueError("Unknown confidence value")

            seen_prediction = items.get((sample_id, drug))
            if seen_prediction:
                if values.index(prediction) > values.index(seen_prediction):
                    # This prection is less important than the one we already have,
                    # so we skip it
                    continue

            items[(sample_id, drug)] = prediction
            if gene and mutation:
                effects_table["ENA_RUN_ACCESSION"].append(sample_id)
                effects_table["DRUG"].append(drug)
                effects_table["GENE"].append(gene)
                effects_table["MUTATION"].append(mutation)
                effects_table["confidence"].append(confidence)
                effects_table["PREDICTION"].append(prediction)

    for (ena_run_accession, drug), prediction in items.items():
        predictions_table["ENA_RUN_ACCESSION"].append(ena_run_accession)
        predictions_table["DRUG"].append(drug)
        predictions_table["PREDICTION"].append(prediction)

    df = pd.DataFrame(predictions_table)
    df.to_csv("dat/tbprofiler_PREDICTIONS.csv", index=False)
    effects_df = pd.DataFrame(effects_table)
    effects_df.to_csv("dat/tbprofiler_EFFECTS.csv", index=False)
