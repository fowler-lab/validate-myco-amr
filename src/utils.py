import numpy
import pandas
from scipy.stats import sem
from tqdm import tqdm

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import matplotlib

import seaborn as sns

sns.set_theme(style="whitegrid", palette="muted")
matplotlib.rcParams.update({"font.size": 7})
plt.rcParams.update({"font.size": 7})


def plot_dilution_boxplot(df, filename=None, savefig=False, exclude_fails=False):

    MIC_VALUES = df.sort_values("DILUTION").MIC.unique()
    df["MIC"] = pandas.Categorical(df["MIC"], categories=MIC_VALUES, ordered=True)
    df = df[["MIC", "OUTCOME", "ENA_RUN_ACCESSION"]].groupby(["OUTCOME", "MIC"]).count()
    df.columns = ["NUMBER"]
    df = df[df.NUMBER > 0]
    df["PROP"] = df["NUMBER"] / df.groupby("OUTCOME")["NUMBER"].transform("sum")
    df.reset_index(inplace=True)
    df.OUTCOME = pandas.Categorical(
        df.OUTCOME, categories=["(S+U)S", "(S+U)R", "RR", "RS"], ordered=True
    )
    df.OUTCOME = df.OUTCOME.sort_values()

    plt.figure(figsize=(3, 1.2))
    palette = {
        "(S+U)S": "#bbbbbb",
        "RR": "#bbbbbb",
        "(S+U)R": "#e41a1c",
        "RS": "#fc9272",
    }
    ax = sns.scatterplot(
        data=df,
        x="MIC",
        y="OUTCOME",
        size="PROP",
        alpha=0.8,
        hue="OUTCOME",
        sizes=(50, 500),
        palette=palette,
    )
    for idx, row in df.iterrows():
        ax.text(
            row["MIC"],
            row["OUTCOME"],
            row["NUMBER"],
            ha="center",
            va="center",
            fontsize=7,
        )
    ax.grid(False)
    ax.get_legend().set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    x_axis = ax.axes.get_xaxis()
    x_axis.label.set_visible(False)
    ax.set_ylim(3.5, -0.5)
    y_axis = ax.axes.get_yaxis()
    y_axis.label.set_visible(False)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(7)
    if filename is not None and savefig:
        plt.savefig(
            "pdf/mic/mic-" + filename + ".pdf", bbox_inches="tight", transparent=True
        )
    plt.close()


def plot_growth_boxplot(df, filename=None, savefig=False, exclude_fails=False):
    plt.rcParams.update({"font.size": 7})
    plt.figure(figsize=(3, 1.2))
    palette = {
        "(S+U)S": "#bbbbbb",
        "RR": "#bbbbbb",
        "(S+U)R": "#e41a1c",
        "RS": "#fc9272",
    }
    ax = sns.boxplot(
        data=df,
        x="POS_AVG_GROWTH",
        y="OUTCOME",
        order=["(S+U)S", "(S+U)R", "RR", "RS"],
        hue="OUTCOME",
        palette=palette,
        fill=False,
    )
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    x_axis = ax.axes.get_xaxis()
    x_axis.label.set_visible(False)
    ax.set_xlim(-2, 102)
    y_axis = ax.axes.get_yaxis()
    y_axis.label.set_visible(False)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(7)
    if filename is not None and savefig:
        plt.savefig("pdf/growth/" + filename, bbox_inches="tight", transparent=True)
    plt.close()


def plot_truthtables(
    results, qualities, filestem=None, savefig=False, exclude_fails=False
):

    for quality in qualities:

        df = results[(results.quality == quality) & (results.dataset == "entire")]

        for idx, row in df.iterrows():

            fig = plt.figure(figsize=(1, 1))
            axes = plt.gca()

            axes.add_patch(
                Rectangle((0, 0), 1, 1, fc="#e41a1c", alpha=0.7, edgecolor=None)
            )
            axes.add_patch(
                Rectangle((0, 1), 1, 1, fc="#bbbbbb", alpha=0.3, edgecolor=None)
            )
            # fc="#4daf4a", alpha=0.1))
            axes.add_patch(
                Rectangle((1, 1), 1, 1, fc="#fc9272", alpha=0.7, edgecolor=None)
            )
            axes.add_patch(
                Rectangle((1, 0), 1, 1, fc="#bbbbbb", alpha=0.3, edgecolor=None)
            )

            axes.set_xlim([0, 2])
            axes.set_ylim([0, 2])

            axes.grid(False)

            axes.set_xticks([0.5, 1.5], labels=["R", "S"], fontsize=7)
            axes.set_yticks([0.5, 1.5], labels=["S+U", "R"], fontsize=7)

            axes.text(0.5, 0.5, int(row["SR"] + row["UR"]), ha="center", va="center")
            axes.text(1.5, 0.5, int(row["SS"] + row["US"]), ha="center", va="center")
            axes.text(0.5, 1.5, int(row["RR"]), ha="center", va="center")
            axes.text(1.5, 1.5, int(row["RS"]), ha="center", va="center")

            if filestem is not None and savefig:
                if exclude_fails:
                    fig.savefig(
                        "pdf/"
                        + "/truthtables/"
                        + filestem
                        + row.quality.lower()
                        + "-"
                        + row.drug
                        + "-exclude-fails.pdf",
                        bbox_inches="tight",
                        transparent=True,
                    )
                else:
                    fig.savefig(
                        "pdf/"
                        + "/truthtables/"
                        + filestem
                        + row.quality.lower()
                        + "-"
                        + row.drug
                        + ".pdf",
                        bbox_inches="tight",
                        transparent=True,
                    )
            plt.close()


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
            try:
                table["U"]
            except KeyError:
                # If there's no unknown predictions, add it with zeros
                table["U"] = 0
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
                try:
                    table["U"]
                except KeyError:
                    # If there's no unknown predictions, add it with zeros
                    table["U"] = 0
                if include_fails:
                    current_sensitivity = table["R"]["R"] / (
                        table["R"]["R"]
                        + table["S"]["R"]
                        + table["U"]["R"]
                        + table["F"]["R"]
                    )

                    current_sensitivity = (
                        table["S"]["S"] + table["U"]["S"] + table["F"]["S"]
                    ) / (
                        table["S"]["S"]
                        + table["R"]["S"]
                        + table["U"]["S"]
                        + table["F"]["S"]
                    )

                else:
                    current_sensitivity = table["R"]["R"] / (
                        table["R"]["R"] + table["S"]["R"] + table["U"]["R"]
                    )

                    current_specificity = (table["S"]["S"] + table["U"]["S"]) / (
                        table["S"]["S"] + table["R"]["S"] + table["U"]["S"]
                    )

                current_ppv = (table["R"]["R"]) / (table["R"]["R"] + table["R"]["S"])

                sensitivity.append(current_sensitivity)
                specificity.append(current_specificity)
                ppv.append(current_ppv)

                row = [
                    set,
                    drug,
                    method,
                    quality,
                    "bootstrap-" + str(i),
                    100 * current_sensitivity,
                    None,
                    100 * current_specificity,
                    None,
                    100 * current_ppv,
                    None,
                ]
                table_rows.append(row)

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
