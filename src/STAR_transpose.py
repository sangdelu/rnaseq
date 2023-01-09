import argparse
import numpy as np
import pandas as pd
import subprocess

AIM_TAB = "/ReadsPerGene.out.tab"


def get_tab(star_dir):
    """

    :param star_dir: String. Path to STAR's results. Should start from root directory.
    :return: list of string. Paths of result table ReadsPerGene.out.tab of all samples.
    """
    table_cmd = "ls " + star_dir + "*" + AIM_TAB
    process = subprocess.Popen(table_cmd, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    tab_lst = output.decode()[:-1].split("\n")
    return tab_lst


def check_strandedness(table):
    table = table.astype({"strand1":'int32',"strand2":'int32'})
    check_table = table[(table["strand1"] >= 20) | (table["strand2"] >= 20)]
    strand1_sum = sum(check_table["strand1"].tolist())
    strand2_sum = sum(check_table["strand2"].tolist())
    #print(strand1_sum)
    #print(strand2_sum)
    if strand1_sum/strand2_sum > 10:
        return "strand1"
    elif strand2_sum/strand1_sum > 10:
        return "strand2"
    else:
        return "unstranded"

def trans_tab(tab_lst, strandness):
    """

    :param tab_lst: list of string. Paths of result table ReadsPerGene.out.tab of all samples.
    :return: Data.Frame. Count matrix containing all genes and samples.
    """
    #strand_dict = {0: "unstranded", 1: "strand1", 2: "strand2"}
    #selected_strand = strand_dict[strandness]
    first_table = pd.read_csv(tab_lst[0],sep="\t", skiprows=4, index_col=0, names=['', 'unstranded', 'strand1', 'strand2'])
    suggest_strand = check_strandedness(first_table)
    if suggest_strand != strandness:
        print("Please double check the strandedness. Looks the one you choosed is not correct.")
        print("Possible correct strand could be: " + suggest_strand)
        exit()
    else:
        res = first_table[[strandness]]
        res = res.rename(columns={strandness: tab_lst[0].split("/")[-2]})
        for t in tab_lst[1:]:
            temp = pd.read_csv(t, sep="\t",skiprows=4, index_col=0, names=['', 'unstranded', 'strand1', 'strand2'])
            temp = temp.rename(columns={strandness: t.split("/")[-2]})
            res = res.merge(temp[[t.split("/")[-2]]], left_index=True, right_index=True, how="outer")
    return res


def main(star_dir, out_path, strandness):
    """

    :param star_dir: String. Path to STAR's results. Should start from root directory.
    :param out_path: String. Path of result count matrix.
    :return:
    """
    if not star_dir.endswith("/"):
        star_dir = star_dir + "/"
    tab_lst = get_tab(star_dir)
    res = trans_tab(tab_lst, strandness)
    res.to_csv(out_path, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="STAR count table transpose")
    parser.add_argument("star_dir")
    parser.add_argument("res_path")
    parser.add_argument("strandness")
    args = parser.parse_args()
    main(args.star_dir, args.res_path, args.strandness)
