from Code.Enum.FileType import FileType
from Code.Extraction.DataExtracter import DataExtracter
from Code.Files.FileReader import FileReader
from Code.Files.FileMaker import FileMaker
from Code.Utils.Ensembl import Ensembl


def splitDict(d):
    split_index = len(d) // 2
    items = list(d.items())
    d1 = items[:split_index]
    d2 = items[split_index:]
    return d1, d2


def filter_out_duplicated_genes(sent_file_path, sent_path_name, in_file_path, in_file_name, out_file_name):
    filtered_lst = []
    sent_genes = FileReader(sent_file_path, sent_path_name).from_file_to_list(False)
    lst = FileReader(in_file_path, in_file_name).from_file_to_list(False)
    # to filter out self-duplicates
    lst_without_duplicates = list(set(lst))
    for gene in lst_without_duplicates:
        if gene not in sent_genes:
            filtered_lst.append(gene)
    FileMaker().fromListToFile(filtered_lst, out_file_name)


def priti_request(file_path, file_name1, file_name2, file_name3):
    data = {}
    worm_ids1 = FileReader(file_path, file_name1, fileType=FileType.CSV).from_file_to_dict()
    print("Done reading worm ids1, items:", len(worm_ids1))
    worm_ids2 = FileReader(file_path, file_name2, fileType=FileType.CSV).from_file_to_dict()
    print("Done reading worm ids2, items:", len(worm_ids2))
    worm_ids3 = FileReader(file_path, file_name3, fileType=FileType.CSV).from_file_to_dict()
    print("Done reading worm ids3, items:", len(worm_ids3))
    worm_ids_lists = [worm_ids1, worm_ids2, worm_ids3]
    for worm_id_list in worm_ids_lists:
        for worm_id in worm_id_list:
            gene_name = Ensembl.get_gene_name_by_gene_id(worm_id)
            gene_description = DataExtracter.get_c_elegans_description_for_gene_id(worm_id)
            dtc = True if "dtc" in gene_description or "distal tip cell" in gene_description else False
            data[worm_id] = (gene_name, gene_description, dtc)
        print("Done with list")

    L3N2 = open("L3N2-updated", "w")
    for gene_id in worm_ids1:
        if gene_id in data:
            L3N2.write(gene_id + "\t" + worm_ids1[gene_id] + "\t" + "\t".join([str(item) for item in data[gene_id]]) + "\n")
    L3N2.close()
    L4N2 = open("L4N2-updated", "w")
    for gene_id in worm_ids2:
        if gene_id in data:
            L4N2.write(gene_id + "\t" + worm_ids2[gene_id] + "\t" + "\t".join([str(item) for item in data[gene_id]]) + "\n")
    L4N2.close()
    L3L4 = open("L3L4-updated", "w")
    for gene_id in worm_ids3:
        if gene_id in data:
            L3L4.write(gene_id + "\t" + worm_ids3[gene_id] + "\t" + "\t".join([str(item) for item in data[gene_id]]) + "\n")
    L3L4.close()
    # common = set(worm_ids1) & set(worm_ids2)
    # print("common ids:", *common, sep="\n")
    # specific_to_one = set(worm_ids1) - set(worm_ids2)
    # print("only in L3&N2:", *specific_to_one, sep="\n")
    # specific_to_two = set(worm_ids2) - set(worm_ids1)
    # print("only in L4&N2:", *specific_to_two, sep="\n")


# priti_request(file_path=r"C:\Users\Liran\PycharmProjects\Research\Code\Utils\priti",
#               file_name1=r"\l3_vs_n2_new.csv.csv",
#               file_name2=r"\l4_vs_n2_new.csv.csv",
#               file_name3=r"\l3_vs_l4_new.csv.csv")
