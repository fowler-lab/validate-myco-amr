import numpy, pandas
from scipy.stats import sem
from tqdm import tqdm


def create_predictions_table(effects, who_drugs):

    samples = []
    table = {"SET": [], "ENA_RUN_ACCESSION": [], "DRUG": [], "PREDICTION": []}
    added_to_table = False

    for idx, row in tqdm(effects.iterrows(), total=effects.shape[0]):
        uid = idx[0]
        drug = idx[1]

        # setup for a new sample
        if uid not in samples:

            # first check if we need to add the previous sample to the table
            if samples:
                for i in predictions:
                    if predictions[i] == "E":
                        value = "S"
                    else:
                        value = predictions[i]
                    table["SET"].append(row.SET)
                    table["ENA_RUN_ACCESSION"].append(working_uid)
                    table["DRUG"].append(i)
                    table["PREDICTION"].append(value)
                    added_to_table = True

            predictions = {drug: "S" for drug in who_drugs}
            added_to_table = False
            samples.append(uid)
            working_uid = uid

        # only update if the predicted effect is more severe
        if row.PREDICTION == "R" and predictions[drug] in ["F", "U", "S"]:
            predictions[drug] = "R"
        elif row.PREDICTION == "F" and predictions[drug] in ["U", "S"]:
            predictions[drug] = "F"
        elif row.PREDICTION == "U" and predictions[drug] == "S":
            predictions[drug] = "U"
        elif row.PREDICTION == "S" and row.EPISTASIS:
            predictions[drug] = "E"

    if not added_to_table:
        for i in predictions:
            table["SET"].append(row.SET)
            table["ENA_RUN_ACCESSION"].append(uid)
            table["DRUG"].append(i)
            table["PREDICTION"].append(predictions[i])

    results = pandas.DataFrame(table)
    results.set_index(["SET", "ENA_RUN_ACCESSION", "DRUG"], inplace=True)
    return results


def build_exact_table(
    table_rows,
    df,
    drug_list,
    quality_list,
    set,
    method,
    debug=False,
    include_fails=False,
):

    for quality in quality_list:
        for drug in drug_list:
            sensitivity = []
            specificity = []
            ppv = []
            if quality == "HIGH":
                foo = df[(df.DRUG == drug) & (df.PHENOTYPE_QUALITY == "HIGH")]
            else:
                foo = df[(df.DRUG == drug)]
            table = pandas.crosstab(
                foo.BINARY_PHENOTYPE, foo.PREDICTION, margins=True, dropna=False
            )
            if include_fails:
                sensitivity = table["R"]["R"] / (
                    table["R"]["R"]
                    + table["S"]["R"]
                    + table["U"]["R"]
                    + table["F"]["R"]
                )
                specificity = (table["S"]["S"] + table["U"]["S"] + table["F"]["S"]) / (
                    table["S"]["S"]
                    + table["R"]["S"]
                    + table["U"]["S"]
                    + table["F"]["S"]
                )
            else:
                sensitivity = table["R"]["R"] / (
                    table["R"]["R"] + table["S"]["R"] + table["U"]["R"]
                )
                specificity = (table["S"]["S"] + table["U"]["S"]) / (
                    table["S"]["S"] + table["R"]["S"] + table["U"]["S"]
                )
            ppv = (table["R"]["R"]) / (table["R"]["R"] + table["R"]["S"])
            if debug:
                print(drug)

            row = [
                set,
                drug,
                method,
                quality,
                "entire",
                100 * numpy.mean(sensitivity),
                0,
                100 * numpy.mean(specificity),
                0,
                100 * numpy.mean(ppv),
                0,
            ]
            row.append(table["R"]["R"])
            row.append(table["S"]["R"])
            row.append(table["U"]["R"])
            row.append(table["R"]["S"])
            row.append(table["S"]["S"])
            row.append(table["U"]["S"])
            row.append(table["All"]["All"])

            if debug:
                print(
                    "%.1f %.1f %.1f"
                    % (
                        100 * numpy.mean(sensitivity),
                        100 * numpy.mean(specificity),
                        100 * numpy.mean(ppv),
                    )
                )
                print(table)
                print()
            table_rows.append(row)

    return table_rows


def build_bootstrapped_table(
    table_rows,
    df,
    drug_list,
    quality_list,
    set,
    method,
    n_bootstraps=50,
    n_sample=500,
    debug=False,
    include_fails=False,
):

    for quality in quality_list:

        for drug in drug_list:

            sensitivity = []
            specificity = []
            ppv = []
            for i in range(n_bootstraps):
                if quality == "HIGH":
                    foo = df[
                        (df.DRUG == drug) & (df.PHENOTYPE_QUALITY == "HIGH")
                    ].sample(n=n_sample, random_state=42 + i, replace=True)
                else:
                    foo = df[(df.DRUG == drug)].sample(
                        n=n_sample, random_state=42 + i, replace=True
                    )

                table = pandas.crosstab(
                    foo.BINARY_PHENOTYPE, foo.PREDICTION, margins=True, dropna=False
                )
                if include_fails:
                    sensitivity.append(
                        table["R"]["R"]
                        / (
                            table["R"]["R"]
                            + table["S"]["R"]
                            + table["U"]["R"]
                            + table["F"]["R"]
                        )
                    )
                    specificity.append(
                        (table["S"]["S"] + table["U"]["S"] + table["F"]["S"])
                        / (
                            table["S"]["S"]
                            + table["R"]["S"]
                            + table["U"]["S"]
                            + table["F"]["S"]
                        )
                    )
                else:
                    sensitivity.append(
                        table["R"]["R"]
                        / (table["R"]["R"] + table["S"]["R"] + table["U"]["R"])
                    )
                    specificity.append(
                        (table["S"]["S"] + table["U"]["S"])
                        / (table["S"]["S"] + table["R"]["S"] + table["U"]["S"])
                    )
                ppv.append((table["R"]["R"]) / (table["R"]["R"] + table["R"]["S"]))

            if debug:
                print(drug)

            row = [
                set,
                drug,
                method,
                quality,
                "bootstrapped50",
                100 * numpy.mean(sensitivity),
                1.96 * 100 * sem(sensitivity),
                100 * numpy.mean(specificity),
                1.96 * 100 * sem(specificity),
                100 * numpy.mean(ppv),
                1.96 * 100 * sem(ppv),
            ]
            for i in range(7):
                row.append(None)

            if debug:
                print(
                    "%.1f %.1f %.1f %.1f %.1f %.1f"
                    % (
                        100 * numpy.mean(sensitivity),
                        1.96 * 100 * sem(sensitivity),
                        100 * numpy.mean(specificity),
                        1.96 * 100 * sem(specificity),
                        100 * numpy.mean(ppv),
                        1.96 * 100 * sem(ppv),
                    )
                )
                print()

            table_rows.append(row)

    return table_rows