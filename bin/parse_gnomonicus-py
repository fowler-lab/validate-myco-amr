import pandas
import pathlib
import json

from collections import defaultdict

predictions = defaultdict(list)
effects = defaultdict(list)

for folder in ["mgit", "ukmyc"]:

    path = pathlib.Path("dat/outputs/")

    for i in (path / folder).glob("*.gnomonicus-out.json"):

        # Exclude the five samples with read naming issues
        if i.stem.split(".")[0] in [
            "ERR4796519",
            "ERR4796408",
            "ERR4796312",
            "ERR4796311",
            "ERR4796303",
        ]:
            print(f"Skipping {i.stem} due to read naming issues")
            continue

        f = open(i)

        sample = i.stem.split(".")[0]

        data = json.load(f)

        for drug, prediction in data["data"]["antibiogram"].items():
            predictions["ENA_RUN_ACCESSION"].append(sample)
            predictions["DRUG"].append(drug)
            predictions["PREDICTION"].append(prediction)

        for drug in data["data"]["effects"]:
            for j in data["data"]["effects"][drug]:
                if "phenotype" in j.keys():
                    continue
                effects["ENA_RUN_ACCESSION"].append(sample)
                effects["DRUG"].append(drug)
                effects["GENE"].append(j["gene"])
                effects["MUTATION"].append(j["mutation"])
                effects["PREDICTION"].append(j["prediction"])
                epistasis = False
                if "expert_rule" in j["evidence"].keys():
                    if "epistasis" in j["evidence"]["expert_rule"]:
                        epistasis = True
                effects["EPISTASIS"].append(epistasis)

predictions = pandas.DataFrame.from_dict(predictions)
effects = pandas.DataFrame.from_dict(effects)


def assign_booleans(row):
    minor_call = False
    is_null = False
    if ":" in row.MUTATION:
        minor_call = True
    if "x" in row.MUTATION:
        is_null = True
    elif "X" in row.MUTATION:
        is_null = True
    return pandas.Series([minor_call, is_null])


effects[["IS_MINOR_ALLELE", "IS_NULL"]] = effects.apply(assign_booleans, axis=1)

effects.set_index(
    ["ENA_RUN_ACCESSION", "DRUG", "GENE", "MUTATION"],
    inplace=True,
    verify_integrity=True,
)
predictions.set_index(
    ["ENA_RUN_ACCESSION", "DRUG"], inplace=True, verify_integrity=True
)

effects.to_csv("dat/RAW_EFFECTS.csv")
predictions.to_csv("dat/RAW_PREDICTIONS.csv")
