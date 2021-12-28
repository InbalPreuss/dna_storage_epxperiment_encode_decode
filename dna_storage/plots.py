import gzip
import json
import re
import uuid
from pathlib import Path
from textwrap import wrap
import os

import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
sns.set_theme(style="ticks")

Path("C:/Users/Inbal/PycharmProjects/DnaStorage/data/testing/plots").mkdir(parents=True, exist_ok=True)


def delete_double_gz() -> pd.DataFrame:
    output_dir = Path("data/testing")
    run_dirs = list(f for f in output_dir.iterdir() if f.name.startswith("["))
    for run_dir in run_dirs:
        files = list(run_dir.iterdir())
        for file in files:
            suffixes = file.suffixes
            if len(suffixes) > 1 and suffixes[-1] == ".gz" and suffixes[-2] == ".gz":
                with gzip.open(file, "rb") as f_in, open(file.with_suffix(""), "wb") as f_out:
                    f_out.writelines(f_in)
                os.remove(file)


def load_sorted_oligos_to_df() -> pd.DataFrame:
    import csv
    output_dir = Path("data/testing")
    # run_dirs = list(f for f in output_dir.iterdir() if f.name.startswith("["))
    # sorted_files = []
    #
    # for run_dir in run_dirs:
    #     try:
    #         sorted_file = [f for f in run_dir.iterdir() if "sort_oligo_results_file" in f.name][0]
    #     except IndexError:
    #         continue
    #     sorted_files.append(sorted_file)
    # print(f'sort_oligo_results_file')
    # info_list = []

    with open("sorted_oligos.csv", 'w') as csvfile:
        fieldnames = ["uid", "barcode", "read_len"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for sorted_file in sorted_files:
            uid = uuid.uuid4()
            try:
                with gzip.open(sorted_file, "rb") as f:
                    for line in f:
                        line = line.strip()
                        writer.writerow({
                            "uid": uid,
                            "barcode": line[:12].decode(),
                            "read_len": len(line[12:]),
                        })
            except:
                print(sorted_file)
                print(line)

    print(f'info_list')

    df_reads = pd.read_csv("C:/Users/Inbal/PycharmProjects/DnaStorage/sorted_oligos.csv")
    # df_reads = pd.read_csv("../sorted_oligos.csv")
    return df_reads


def draw_reads_histograms(df: pd.DataFrame):
    plt.subplots()
    fig = sns.histplot(df["read_len"], kde=False)
    fig.set_ylabel("count")
    fig = plt.gcf()
    fig.savefig(Path("C:/Users/Inbal/PycharmProjects/DnaStorage/data/testing/plots") / "read_len_count")

    groups = df.groupby(["uid", "barcode"])
    reads_per_barcode = []
    for _, group in groups:
        reads_per_barcode.append({"reads": len(group)})

    plt.subplots()
    # how many lines (sequences) was read per each barcode.
    # (That's why we have concentrations around 20, 50, 100, 1000)
    df_reads_per_barcode = pd.DataFrame(reads_per_barcode)
    fig = sns.histplot(df_reads_per_barcode["reads"], kde=False, binwidth=1)
    fig.set_ylabel("count")
    fig.set_xlabel("number of reads per barcode")
    fig = plt.gcf()
    fig.savefig(Path("data/testing/plots") / "number of reads per barcode")


def load_data_to_df() -> pd.DataFrame:
    output_dir = Path("C:/Users/Inbal/PycharmProjects/DnaStorage/data/testing")

    run_dirs = list(f for f in output_dir.iterdir() if f.name.startswith("["))

    json_files = []

    for run_dir in run_dirs:
        try:
            json_file = [f for f in run_dir.iterdir() if f.suffix == ".json"][0]
        except IndexError:
            continue
        json_files.append(json_file)

    info_list = []

    for json_file in json_files:
        with json_file.open("rb") as f:
            info = json.load(f)
        info_list.append(info)

    df = pd.DataFrame(info_list)
    df["output_dir"] = df["output_dir"].apply(lambda x: x.split("data/testing/")[1].split(" trial")[0])
    return df


def draw_boxplots(df: pd.DataFrame, percentage: bool = False):
    trials_group = df.groupby([
        'number_of_oligos_per_barcode',
        'number_of_sampled_oligos_from_file',
        'size',
        'bits_per_z'
    ])

    errors = ["substitution_error", "deletion_error", "insertion_error"]
    y_values = ["levenshtein_distance", "levenshtein_distance_sigma_before_rs", "levenshtein_distance_sigma_after_rs_payload", "levenshtein_distance_sigma_after_rs_wide"]
    for y in y_values:
        if y == "levenshtein_distance":
            palette = "Purples"
        else:
            palette = "Blues"
        for idx, trial_group in trials_group:
            fig, axes = plt.subplots(nrows=3)
            fig.suptitle("\n".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0] + f"\n{y}".replace("_", " "), 71)))
            fig.subplots_adjust(top=0.85, hspace=0.5)
            for ax, error in zip(axes, errors):
                zero_cols = [e for e in errors if e != error]
                df_for_err = trial_group
                for col in zero_cols:
                    df_for_err = df_for_err[df_for_err[col] == 0]
                if percentage:
                    ax = sns.barplot(x=error, y=y, data=df_for_err, estimator=estimate, ax=ax, palette=palette)
                    ax.set(ylabel="D% [levenshtein]")
                else:
                    ax = sns.boxplot(x=error, y=y, data=df_for_err, ax=ax, palette=palette)
                    ax = sns.swarmplot(x=error, y=y, data=df_for_err, color=".25", ax=ax)
                    ax.set(ylabel="D [levenshtein]")
            if percentage:
                fig.savefig(Path("data/testing/plots") / " ".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0].replace("[", " ").replace("]", " ") + f"{y} success", 71)))
                plt.close(fig)
            else:
                fig.savefig(Path("data/testing/plots") / "".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0].replace("[", " ").replace("]", " ") + f"{y}", 71)))
                plt.close(fig)


def draw_boxplots_all_samples(df: pd.DataFrame, percentage: bool = False):
    trials_group = df.groupby([
        'number_of_oligos_per_barcode',
        'size',
        'bits_per_z'
    ])

    errors = ["substitution_error", "deletion_error", "insertion_error"]
    y_values = {"levenshtein_distance": "input_text_len",
                "levenshtein_distance_sigma_before_rs": "input_data_encoder_results_file_len",
                "levenshtein_distance_sigma_after_rs_payload": "input_data_encoder_without_rs_payload_len",
                "levenshtein_distance_sigma_after_rs_wide": "input_data_encoder_without_rs_wide_len"
                }
    for y_value in y_values.items():
        df[y_value[0]] = df.apply(lambda x: x[y_value[0]] / x.get(y_value[1], 1), axis=1)
    for y in y_values:
        if y == "levenshtein_distance":
            palette = "Purples"
        else:
            palette = "Blues"
        for idx, trial_group in trials_group:
            fig, axes = plt.subplots(nrows=3)
            title = trial_group["output_dir"].iloc[0].split("[ errors")[0] + f"\n{y}".replace("_", " ")
            title = re.sub(r'\[ number of oligos sampled after synthesis[^\S\n\t]+\d+[^\S\n\t]+\]', '', title)
            fig.suptitle("\n".join(wrap(title, 71)))
            fig.subplots_adjust(top=0.85, hspace=0.5, right=0.8)
            for ax_idx, (ax, error) in enumerate(zip(axes, errors)):
                zero_cols = [e for e in errors if e != error]
                df_for_err = trial_group
                for col in zero_cols:
                    df_for_err = df_for_err[df_for_err[col] == 0]
                if percentage:
                    ax = sns.barplot(x=error, y=y, hue="number_of_sampled_oligos_from_file", data=df_for_err, estimator=estimate, ax=ax, palette=palette)
                    ax.set(ylabel="D% [levenshtein]")
                else:
                    ax = sns.boxplot(showfliers=True,x=error, y=y, hue="number_of_sampled_oligos_from_file", data=df_for_err, ax=ax, palette=palette)

                    # ax = sns.swarmplot(x=error, y=y, hue="number_of_sampled_oligos_from_file", data=df_for_err, ax=ax, dodge=True, palette=palette)
                    ax.set(ylabel="D [levenshtein]")
                if ax_idx != 0:
                    ax.get_legend().remove()
                else:
                    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            if percentage:
                fig.savefig(Path("C:/Users/Inbal/PycharmProjects/DnaStorage/data/testing/plots") / " ".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0].replace("[", " ").replace("]", " ") + f"{y} success", 71)))
                plt.close(fig)
            else:
                fig.savefig(Path("C:/Users/Inbal/PycharmProjects/DnaStorage/data/testing/plots") / "".join(wrap(trial_group["output_dir"].iloc[0].split("[ errors")[0].replace("[", " ").replace("]", " ") + f"{y}", 71)))
                plt.close(fig)
    x=3

def draw_sampled_vs_error(df: pd.DataFrame):
    fig, ax = plt.subplots()
    ax = sns.boxplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, ax=ax, palette="Purples")
    ax = sns.swarmplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, color=".25", ax=ax)
    ax.set(xlabel="number_of_sampled_oligos_from_file".replace("_", " "), ylabel="D [levenshtein]")
    fig.savefig(Path("data/testing/plots") / "number of sampled oligos and levenshtein distance")
    plt.close(fig)


def estimate(x):
    res = sum(x == 0) * 100.0 / len(x)
    return res


def draw_zero_error_percentage(df: pd.DataFrame):
    fig, ax = plt.subplots()
    sns.barplot(x="number_of_sampled_oligos_from_file", y="levenshtein_distance", data=df, estimator=estimate, ax=ax, palette="Blues")
    ax.set(xlabel="number_of_sampled_oligos_from_file".replace("_", " "), ylabel="D% [levenshtein]")
    fig.savefig(Path("data/testing/plots") / "number of sampled oligos and levenshtein distance success")
    plt.close(fig)


def draw_error_per_number_of_sampled_oligos(df: pd.DataFrame):
    for idx, row in df.iterrows():
        if row["substitution_error"] != 0:
            df.loc[idx, "error_type"] = "substitution_error"
            df.loc[idx, "error"] = row["substitution_error"]
        elif row["deletion_error"] != 0:
            df.loc[idx, "error_type"] = "deletion_error"
            df.loc[idx, "error"] = row["deletion_error"]
        elif row["insertion_error"] != 0:
            df.loc[idx, "error_type"] = "insertion_error"
            df.loc[idx, "error"] = row["insertion_error"]
        else:
            df.loc[idx, "error_type"] = "no_error"
            df.loc[idx, "error"] = 0

    fig, ax = plt.subplots()
    ax = sns.catplot(
        data=df,
        x="number_of_sampled_oligos_from_file",
        y="levenshtein_distance",
        row="error_type",
        col="error",
        palette="Purples",
        kind="bar",
    )
    ax.set(xlabel="error_per_number_of_sampled_oligos".replace("_", " "), ylabel="D [levenshtein]")
    fig.savefig(Path("data/testing/plots") / "error_per_number_of_sampled_oligos")

    fig = px.bar(
        df, x="number_of_sampled_oligos_from_file", y="levenshtein_distance",
        facet_row="error_type", facet_col="error",
        barmode="group",
    )
    fig.show()

def main():
    plt.ion()
    # df_reads = load_sorted_oligos_to_df()
    # draw_reads_histograms(df=df_reads)
    
    df = load_data_to_df()
    # draw_zero_error_percentage(df=df)
    # draw_boxplots(df=df)
    # draw_boxplots(df=df, percentage=True)
    draw_boxplots_all_samples(df=df)
    draw_boxplots_all_samples(df=df, percentage=True)
    # draw_sampled_vs_error(df=df)
    # draw_error_per_number_of_sampled_oligos(df=df)
    # input("Hit enter to terminate")


if __name__ == '__main__':
    main()
