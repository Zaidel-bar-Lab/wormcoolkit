from Code.Enum.FileMode import FileMode
from Code.Enum.FileType import FileType
from Code.Extraction.DataExtracter import DataExtracter
from Code.Files.FileMaker import FileMaker
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.BioPython import BioPython
from Code.Utils.Ensembl import Ensembl
from Code.Utils.Strings import Strings
from Test.TestFunctions import TestFunctions
#import pandas as pd

import socket


class executor:

    def __init__(self):
        self.mmp_data_by_gene_name = FileReader(FileReader.research_path + r"\Data",
                                                r"\mmp_mut_strains_data_Mar14.txt"). \
            from_MMP_file_to_dict_with_listed_values(9, [11, 18, 12], delete_first=True)

    # creating a file to restore the c.elegans gene id (WB...) and the number of its conserved domains
    @staticmethod
    def c_elegans_conserved_domains_file():
        print("def: c_elegans_conserved_domains_file")
        c_elegans_genes_id_WB = FileReader(FileReader.research_path + r"\Executors",
                                           r"\data-230619-fixed-ratio-with-C-elegans-phenotypes",
                                           FileType.TSV).get_c_elegans_genes_from_output(1)
        # c_elegans_genes_id_wb_number_dict = BioPython.get_genes_id(c_elegans_genes_id_WB)
        # FileMaker().fromDictToUniqueFile(c_elegans_genes_id_wb_number_dict, "relevant-c-elegans-genes-and-their-id-110619")

        de = DataExtracter()
        c_elegans_genes_and_conserved_domains = de.get_conserved_domains(c_elegans_genes_id_WB)
        print("The length of the genes_id_domain is", len(c_elegans_genes_and_conserved_domains))

        fm = FileMaker()
        fm.from_dict_to_file(c_elegans_genes_and_conserved_domains, "c-elegans-genes-and-conserved-domains-230619")

    # creating a file to restore the c.elegans gene id (WB...) and the number of its conserved domains
    @staticmethod
    def human_conserved_domains_file():
        print("def: human_conserved_domains_file")
        fr2 = FileReader(FileReader.research_path + r"\Executors",
                         r"\genes-orthologs-phenotypes-filteredBySize-100619",
                         FileType.TSV)
        genes_names = fr2.get_genes_list()

        de = DataExtracter()
        genes_id = de.from_human_genes_names_to_human_genes_id(genes_names)
        print("The length of the genes ids list is", len(genes_id))
        print("human genes ids for example:", genes_id[0], genes_id[1])

        human_genes_and_conserved_domains = de.get_conserved_domains(genes_id)

        fm = FileMaker()
        fm.from_dict_to_file(human_genes_and_conserved_domains, "human-genes-and-conserved-domains-110619")

    # using the (1) C.elegans genes conserved domains file, (2) human genes conserved domains file, (3) most
    # recent data file, and (4) human gene id-gene name file, extracts the info needed, conducts the calculation
    # needed and creates a file of data that includes the ratio of conserved domains
    @staticmethod
    def add_domains_score_info_to_file():
        print("def: add_domains_score_info_to_file")
        # first we need the C.elegans genes-number of domains dictionary
        fr1 = FileReader(FileReader.research_path + r"\Executors",
                         r"\c-elegans-genes-and-conserved-domains-110619",
                         FileType.TSV)
        c_elegans_genes_domains_dic = fr1.from_file_to_dict(0, 1)
        tf1 = TestFunctions("from_file_to_dict to read c.elegans genes and domains",
                            dictionary=c_elegans_genes_domains_dic)
        tf1.checkSize()
        tf1.print_first_lines_in_dict(2)

        # second we need the human genes-number of domains dictionary
        fr2 = FileReader(FileReader.research_path + r"\Executors",
                         r"\human-genes-and-conserved-domains-110619",
                         FileType.TSV)
        human_genes_domains_dic = fr2.from_file_to_dict(0, 1)
        tf2 = TestFunctions("from_file_to_dict to read human genes and domains",
                            dictionary=human_genes_domains_dic)
        tf2.checkSize()
        tf2.print_first_lines_in_dict(2)

        # then we'll need the dictionary to convert human gene name to its c.elegans orthologs
        fr3 = FileReader(FileReader.research_path + r"\Executors",
                         r"\genes-orthologs-phenotypes-filteredBySize-090619",
                         FileType.TSV)
        orthologs = fr3.from_file_to_dict(0, 1)
        tf3 = TestFunctions("from_file_to_dict to read human genes and orthologs",
                            dictionary=orthologs)
        tf3.checkSize()
        tf3.print_first_lines_in_dict(2)

        # and last we need a dictionary to convert human gene name to human gene id
        fr4 = FileReader(FileReader.research_path + r"\Data",
                         r"\human_gene_id-gene_name.txt",
                         FileType.TSV)
        human_genes_names_ids_dict = fr4.from_file_to_dict(1, 0, True)
        tf4 = TestFunctions("from_file_to_dict to read human genes names and ids",
                            dictionary=human_genes_names_ids_dict)
        tf4.checkSize()
        tf4.print_first_lines_in_dict(2)

        human_genes_id_orthologs_dic = DataExtracter.add_conserved_domains_info(c_elegans_genes_domains_dic,
                                                                                human_genes_domains_dic,
                                                                                orthologs,
                                                                                human_genes_names_ids_dict)
        tf5 = TestFunctions("add_conserved_domains_info",
                            dictionary=human_genes_id_orthologs_dic)
        tf5.checkSize()
        tf5.print_first_lines_in_dict(7)

        # now lets write all this back to a file, but for that we need two things:
        # first, a converter from gene id to gene name:
        human_genes_ids_names_dict = fr4.from_file_to_dict(0, 1, True)
        tf6 = TestFunctions("from_file_to_dict to read human genes ids and names",
                            dictionary=human_genes_ids_names_dict)
        tf6.checkSize()
        tf6.print_first_lines_in_dict(2)

        # and second, a dictionary of human genes' names and conditions
        fr5 = FileReader(FileReader.research_path + r"\Executors",
                         r"\genes-orthologs-phenotypes-filteredBySize-100619",
                         FileType.TSV)
        genes_names_conditions = fr5.from_file_to_dict(0, 2)
        tf7 = TestFunctions("from_file_to_dict for dictionary of human genes' names and conditions",
                            dictionary=genes_names_conditions)
        tf7.checkSize()
        tf7.print_first_lines_in_dict(2)

        fm = FileMaker()
        fm.from_two_dict_of_different_keys_to_one_file(human_genes_id_orthologs_dic,
                                                       genes_names_conditions,
                                                       human_genes_ids_names_dict,
                                                       "genes_id-genes_names-orthologs-conditions-110619")

    @staticmethod
    def get_only_genes_with_human_ortholog():
        print("def: get_only_genes_with_human_ortholog ")
        print("*** Function \"get_only_genes_with_human_ortholog\" has started ***")

        accession_numbers_and_hit_ids, hit_ids_and_hsps = DataExtracter.get_hit_ids_for_accession_numbers(
            31,
            FileReader.research_path + r"\Test\Files\accessions-and-sequences",
            r"\results-part-")
        print("accessions numbers and hit ids and hsp have been obtained")

        # hit_ids_and_accessions = DataExtracter.from_plural_valued_dict_to_plural_valued_reversed_dict(accession_numbers_and_hit_ids)
        # hit_ids = hit_ids_and_accessions.keys()
        # print("There are " + str(len(hit_ids_and_accessions)) + " hit ids")
        # hit_ids_and_gene_names = BioPython.get_gene_name_from_protein_accession(hit_ids)
        #
        # print("Gene names for the hit ids have been obtained!")
        # tf1 = TestFunctions("get_gene_name_from_protein_accession",
        #                     dictionary=hit_ids_and_gene_names)
        # tf1.checkSize()
        # tf1.print_first_lines_in_dict(2)

        # FileMaker().from_dict_to_file(hit_ids_and_gene_names, "blast_hit_ids_and-c-elegans-gene_names0-30")
        # FileMaker().from_dict_to_file(hit_ids_and_hsps, "blast_hit_ids_and_hsp_scores0-30")

        hit_ids_and_gene_names = FileReader(FileReader.research_path + r"\Executors",
                                            r"\blast-hit-ids-and-c-elegans-gene-names0-30",
                                            FileType.TSV).from_file_to_dict(0, 1)

        tf1 = TestFunctions("get_gene_name_from_protein_accession",
                            dictionary=hit_ids_and_gene_names)
        tf1.checkSize()
        tf1.print_first_lines_in_dict(2)

        gene_ids_and_accession_numbers = FileReader(FileReader.research_path + r"\Executors",
                                                    r"\extra_genes_ids_number_and_chosen_accessions-120619",
                                                    FileType.TSV).from_file_to_dict(0, 1, False)
        tf2 = TestFunctions("from_file_to_dict to get genes id and the accession numbers",
                            dictionary=gene_ids_and_accession_numbers)
        tf2.checkSize()
        tf2.print_first_lines_in_dict(2)

        gene_names_WB_and_gene_ids = FileReader(FileReader.research_path + r"\Executors",
                                                r"\fixed-relevant-c-elegans-genes-and-their-id-110619",
                                                FileType.TSV).from_file_to_dict(0, 1)

        tf3 = TestFunctions("from_file_to_dict to get genes names and the gene ids",
                            dictionary=gene_names_WB_and_gene_ids)
        tf3.checkSize()
        tf3.print_first_lines_in_dict(2)

        # last one - gene names to their homoloug

        gene_ids_WB_and_human_names = FileReader(FileReader.research_path + r"\Executors",
                                                 r"\genes_id-genes_names-orthologs-conditions-110619",
                                                 FileType.TSV).get_celegans_id_to_human_name_dict()

        tf4 = TestFunctions("get_celegans_id_to_human_name_dict",
                            dictionary=gene_ids_WB_and_human_names)
        tf4.checkSize()
        tf4.print_first_lines_in_dict(2)

        gene_ids_WB_and_true_homologs = {}

        WBGenes = gene_ids_WB_and_human_names.keys()
        print("There are", len(WBGenes), "WBGenes")
        for WBGene in WBGenes:
            gene_name = gene_ids_WB_and_human_names[WBGene]
            print("human gene name for", WBGene, ":", gene_name)
            try:
                gene_id = gene_names_WB_and_gene_ids[WBGene]
                print("c-elegans gene id for", gene_name, ":", gene_id)
            except Exception as e:
                print("Exception in get_only_genes_with_human_ortholog:", e)
                print("couldn't find the gene id for", gene_name)
                continue
            try:
                accession_number = gene_ids_and_accession_numbers[gene_id]
            except Exception as e:
                print("Exception in get_only_genes_with_human_ortholog:", e)
                print("couldn't find the accession number for", gene_id)
                continue
            print("c-elegans accession number for", gene_id, ":", accession_number)
            try:
                hit_ids = accession_numbers_and_hit_ids[accession_number]
                print("hit ids for " + accession_number, ":", hit_ids)
            except KeyError:
                print("couldn't find hit ids for", accession_number)
                continue

            for hit_id in hit_ids:
                try:
                    blasted_gene_name = hit_ids_and_gene_names[hit_id]
                except KeyError:
                    print("couldn't find", hit_id, "in the hit_ids_and_gene_names_dict")
                    continue
                print(blasted_gene_name, gene_name)
                if blasted_gene_name == gene_name:
                    gene_ids_WB_and_true_homologs[WBGene] = gene_name

        tf5 = TestFunctions("unplaced function to find true homologs",
                            dictionary=gene_ids_WB_and_true_homologs)
        tf5.checkSize()
        tf5.print_first_lines_in_dict(2)

        FileMaker().from_dict_to_file(gene_ids_WB_and_true_homologs, "true-homologs0-30")

        human_genes_and_c_elegans_WB_id = FileReader(FileReader.research_path + r"\Executors",
                                                     r"\true-homologs0-30",
                                                     FileType.TSV).from_file_to_dict_with_plural_values(1, 0)
        print("There are", len(human_genes_and_c_elegans_WB_id), "human genes with homologs")

    # uses the true homologs files created by the former function to filter the list of genes
    @staticmethod
    def filterGenesAccordingToReversedBlast(results_file_names: list, data_file_name: str, new_data_file_name: str):
        print("def: filterGenesAccordingToReversedBlast")
        # first we make a dictionary of human genes as keys and their true C.elegans orthologs as list of values
        human_genes_and_c_elegans_orthologs = {}
        for file_name in results_file_names:
            f = open(file_name)
            for line in f:
                genes = line.rstrip("\n").split("\t")
                c_elegans_gene = genes[0]
                human_gene = genes[1]
                DataExtracter.add_to_dictionary(human_genes_and_c_elegans_orthologs, human_gene, c_elegans_gene)
            f.close()
        FileMaker().from_plural_valued_dict_to_file(human_genes_and_c_elegans_orthologs, "human-and-c-elegans-genes")

        # now we filter the data file
        data_file = open(data_file_name)
        new_data_file = open(new_data_file_name, FileMode.WRITE.value)
        for row in data_file:
            orthologs_to_be_written = []
            line = row.rstrip("\n").split("\t")
            human_gene_id = line[0]
            human_gene_name = line[1]
            c_elegans_orthologs = FileReader.fromStringToTuplesList(line[2])
            phenotype = line[3]
            if human_gene_name in human_genes_and_c_elegans_orthologs:
                true_orthologs = human_genes_and_c_elegans_orthologs[human_gene_name]
                for gene_tuple in c_elegans_orthologs:
                    maybe_ortholog = gene_tuple[0]  # the c_elegans_name
                    if maybe_ortholog in true_orthologs:
                        orthologs_to_be_written.append(gene_tuple)
            else:
                print("human gene name", human_gene_name, "wasn't found in the dictionary of true orthologs :(")
                print("meaning no one of his claimed orthologs was found to be real... check it")
                continue
            new_data_file.write(human_gene_id + "\t" +
                                human_gene_name + "\t" +
                                str(orthologs_to_be_written).strip('[]') + "\t" +
                                phenotype + "\n")
        data_file.close()
        new_data_file.close()

    @staticmethod
    def get_human_genes_and_true_orthologs():
        print("def: get_human_genes_and_true_orthologs")
        human_genes_and_c_elegans_WB_id = FileReader(FileReader.research_path + r"\Executors",
                                                     r"\true-homologs0-18",
                                                     FileType.TSV).from_file_to_dict_with_plural_values(1, 0)
        print("There are", len(human_genes_and_c_elegans_WB_id), "human genes with homologs")
        FileMaker().from_plural_valued_dict_to_file(human_genes_and_c_elegans_WB_id, "human-genes-and-true-orthologs")

    @staticmethod
    def check_extra_genes():
        print("def: check_extra_genes")
        # obtaining the new C.elegans genes id (number) list
        c_elegans_current_genes_id_number = FileReader(FileReader.research_path + r"\Executors",
                                                       r"\fixed-relevant-c-elegans-genes-and-their-id-110619",
                                                       FileType.TSV).get_genes_list(1)

        # obtaining the former C.elegans genes id (number) list
        c_elegans_former_genes_id_number = FileReader(FileReader.research_path + r"\Test\Files",
                                                      r"\relevant-c-elegans-genes-and-their-id",
                                                      FileType.TSV).get_genes_list(1)

        print("length of now-genes is:", len(c_elegans_current_genes_id_number), "and length of former genes is:",
            len(c_elegans_former_genes_id_number))
        extra_genes = []
        for gene in c_elegans_current_genes_id_number:
            if gene not in c_elegans_former_genes_id_number:
                print(gene + " does not exist in the former list")
                extra_genes.append(gene.rstrip("\n"))
        print("There are", len(extra_genes), "extra genes ids_number that needs a blast")

        fm = FileMaker()
        fm.fromListToFile(extra_genes, "extra-c-elegans-gened-id-number-120619")

    @staticmethod
    def from_accessions_to_blast():
        print("def: from_accessions_to_blast")
        # extra_genes_id_number_and_accessions_dict = FileReader(FileReader.research_path + r"\Data",
        #                                                        r"\extra-gene-ids-number-accession-numbers.txt",
        #                                                        FileType.TSV).from_file_to_dict_with_plural_values(0,1)
        # print("Size of genes and accessions is " + str(len(extra_genes_id_number_and_accessions_dict)))
        # extra_genes_id_number_and_chosen_accessions = DataExtracter.from_multiple_accessions_to_one(
        #     extra_genes_id_number_and_accessions_dict)
        # print("Size of genes and chosen accessions is " + str(len(extra_genes_id_number_and_chosen_accessions)))
        # FileMaker().from_dict_to_file(extra_genes_id_number_and_chosen_accessions, "extra_genes_ids_number_and_
        # chosen_accessions-120619")

        # accessions = FileReader(FileReader.research_path + r"\Executors",
        #                         r"\extra_genes_ids_number_and_chosen_accessions-120619",
        #                         FileType.TSV).get_genes_list(1)
        # print("Got " + str(len(accessions)) + " accessions! let the work begin...")
        #
        # bp = BioPython()
        # accessions_and_seqs = bp.make_accession_number_and_seq_dict(accessions, "ORIGIN", "translation=")
        # FileMaker().from_dict_to_file(accessions_and_seqs, "extra_chosen_accessions_and_sequences")

        accessions_and_seqs = FileReader(FileReader.research_path + r"\Test\Files\accessions-and-sequences",
                                       r"\accessions-and-sequences-part-helper",
                                       FileType.TSV).from_file_to_dict(0, 1)
        BioPython().blastp_by_accessions("blastp", "nr", accessions_and_seqs)

    # receives (1) data file path, (2) data file name, (3) url address, (4) new file name, (5) boolean value to indicate
    # whether or not we need to disregard the first line of the data, it retrieves the worm genes names from the data,
    # pass it to function that returns a dictionary of said genes and phenotypes, and copies the data with phenotype
    # to each worm gene to the new file
    @staticmethod
    def add_c_elegans_phenotypes(data_file_path, data_file_name, url, new_file, delete_first_line: bool = False):
        print("def: add_c_elegans_phenotypes")
        # first we achieve a list of all relevant C.elegans genes
        cElegansGenes = FileReader(data_file_path, data_file_name, FileType.TSV).get_genes_list(2)
        print("List of genes with", len(cElegansGenes), "have been obtained:", cElegansGenes[0], cElegansGenes[1])

        genes_and_phenotypes = DataExtracter().add_c_elegans_phenotypes(cElegansGenes, url)

        in_file = open(data_file_path + data_file_name)
        out_file = open(new_file, FileMode.WRITE.value)
        if delete_first_line:
            in_file.readline()
        for line in in_file:
            info = line.rstrip("\n").split("\t")
            c_elegans_gene = info[2]
            if c_elegans_gene in genes_and_phenotypes:
                phenotypes: list = genes_and_phenotypes[c_elegans_gene]
            else:
                phenotypes = ["Not Found"]
            out_file.write(line.rstrip("\n") + "\t" + ", ".join(phenotypes) + "\n")
        in_file.close()
        out_file.close()

    # receives (1) a dictionary of human-worm pairs, if by file or if by a human-worm-id dictionary and for each pair,
    # and the species whose genes we wish to blast, and returns whether the reversed blast confirms this pair as
    # orthologous.
    @staticmethod
    def pair_pipeline(list_of_pairs_file_path="", list_of_pair_file_name="", dic_of_optional_orthologs: dict = None,
                      key_species="C.elegans"):
        print("def: pair_pipeline")
        print(list_of_pairs_file_path, list_of_pair_file_name, dic_of_optional_orthologs,key_species)
        if not dic_of_optional_orthologs:  # read dict from file
            gene_name_to_gene_name_dict = FileReader(list_of_pairs_file_path, list_of_pair_file_name,
                                                     FileType.TSV).from_file_to_dict_with_plural_values(0, 1, True)
        else:
            reversed_dic = DataExtracter().from_plural_valued_dict_to_plural_valued_reversed_dict(dic_of_optional_orthologs)
            print("reversed dic:", reversed_dic)
            gene_name_to_gene_name_dict = DataExtracter.convert_dic(reversed_dic, True, key_species)

        print("before blast we got those ",len(gene_name_to_gene_name_dict),"results:",gene_name_to_gene_name_dict)
        print("List of all pairs have been obtained! we have", len(gene_name_to_gene_name_dict), "genes to check")

        true_matches = {}
        de, bp = DataExtracter(), BioPython()
        #print("***BLAST***")
        #bp.build_genes_hit_ids_dictionary([gene_name for gene_name in gene_name_to_gene_name_dict], key_species)
        #if not bp.blast_results:
        #    return true_matches
        for key_gene_name in gene_name_to_gene_name_dict:
            ortholog_genes_names: list = gene_name_to_gene_name_dict[key_gene_name]
            for ortholog_gene_name in ortholog_genes_names:

                c_elegans_gene_name = key_gene_name if key_species == "C.elegans" else ortholog_gene_name
                human_gene_name = ortholog_gene_name if key_species == "C.elegans" else key_gene_name
                print("Now working on", key_gene_name, "and", ortholog_gene_name)

                #result = de.check_reversed_blast_hit_ids(bp, key_gene_name, ortholog_gene_name, key_species)
                result = True
                if result:
                    print("first row of the if")
                    genes_tuple = (ortholog_gene_name, key_gene_name)
                    print(ortholog_gene_name , key_gene_name , human_gene_name)
                    status_tuple = de.get_status_tuple(c_elegans_gene_name, human_gene_name, key_species)
                    print(status_tuple)
                    #if status_tuple[0] != -1:
                    true_matches[genes_tuple] = status_tuple
                    print("The domains ratio, number of sources, human and C.elegans gene length: ", status_tuple)
        return true_matches , gene_name_to_gene_name_dict

    @staticmethod
    def get_result_list(true_matches_dictionary, false_matches_dictionary):
        print("def: get_result_list")
        print(true_matches_dictionary)
        print(false_matches_dictionary)
        true_result_list, false_result_list = [], []
        if true_matches_dictionary:
            for genes_tuple in true_matches_dictionary:
                data = true_matches_dictionary[genes_tuple]
                true_result_list.append([genes_tuple[0], genes_tuple[1], data[0], data[1], data[2][0], data[2][1], data[3]])
        if false_matches_dictionary:
            for gene in false_matches_dictionary:
                false_result_list.append([gene, false_matches_dictionary[gene]])
        print(true_result_list, false_result_list)
        return true_result_list, false_result_list

    @staticmethod
    def filterGenesTest():
        print("def: filterGenesTest")
        DataExtracter.filter_genes_by_name(FileReader.research_path + r"\Executors",
                                           r"\data-250619-fixedDomains-short-lethalFiltered",
                                           FileReader.research_path + r"\Data",
                                           r"\data-010719-sorted",
                                           "data-010719-filteredGenesByUnwantedGenes")

    # receives (1) file type to know whether to print to a file or to console, read all variants and extracts for each
    # variant the sequence of its human gene and its C.elegans ortholog, runs pairwise alignment and returns data
    # regarding the alignments of the two sequences. also filter out all genes for which the amino acid mutated is not
    # conserved
    # not used method, outdated 150919
    @staticmethod
    def get_variants_dict(file_type=FileType.CONSOLE):
        print("def: get_variants_dict")
        human_genes_and_sequences = {}
        human_genes_and_variants = FileReader(FileReader.research_path + r"\Data\variants",
                                              r"\all-missense-nonsense").fromHGMDtoDict()
        data = FileReader(FileReader.research_path + r"\Data",
                          r"\data-250619-fixed-domains").readData(1, False)
        genes_names = human_genes_and_variants.keys()

        orthologs_dic = FileReader(FileReader.research_path + r"\Data",
                                   r"\data-250619-fixed-domains").from_file_to_dict_with_plural_values(1, 2, False)
        c_elegans_id_and_accessions = FileReader(
            FileReader.research_path + r"\Test\Files",
            r"\c-elegans-genes-and-longest-accession_number_complete").from_file_to_dict(0, 1, False)
        # new file: FileReader.research_path + r"\Test\Files\c-elegans-gene-ids-and-accession-numbers"

        accessions_and_sequences = FileReader(FileReader.research_path + r"\Test\Files",
                                              r"\all-accession-numbers-and-sequences").from_file_to_dict(0, 1)
        # new file: FileReader.research_path + r"\Extraction\c-elegans-accession-numbers-and-sequences""

        mmp_data_by_gene_name = FileReader(FileReader.research_path + r"\Data",
                                           r"\mmp_mut_strains_data_Mar14.txt"). \
            from_MMP_file_to_dict_with_listed_values(9, [11, 18, 12], delete_first=True)

        f = open("data-with-variants-230919", FileMode.WRITE.value) if file_type != FileType.CONSOLE else ""
        print(human_genes_and_variants)
        print("Number of genes:", len(human_genes_and_variants))
        for human_gene_name in genes_names:
            human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
            human_seq = HttpRequester().get_human_protein_from_uniprot(human_gene_id)
            human_genes_and_sequences[human_gene_name] = human_seq
            mmp_data = mmp_data_by_gene_name[human_gene_name] if human_gene_name in mmp_data_by_gene_name \
                else "Not mention in MMP"
            print("Human gene name:", human_gene_name, ", seq:", human_seq)
            orthologs_id_WB = orthologs_dic[human_gene_name]
            print("Ortholog gene id WB:", orthologs_id_WB)
            for ortholog_id_WB in orthologs_id_WB:
                ortholog_id_number = Ensembl.get_ncbi_id_by_gene_id(ortholog_id_WB)
                if not ortholog_id_number:
                    ortholog_id_number = BioPython.get_gene_id_by_entrez(ortholog_id_WB)
                    if not ortholog_id_number:
                        print("Couldn't find C.elegans' id (number), moving on to the next gene")
                        continue
                print("id:", ortholog_id_number)
                try:
                    c_elegans_accession_number = c_elegans_id_and_accessions[ortholog_id_number]
                except Exception as e:
                    print("Exception in get_variants_dict:", e)
                    print("Couldn't find C.elegans' accession number, moving on to the next gene")
                    # right now no way to get accession number other that DAVID
                    continue
                print("accession: " + c_elegans_accession_number)
                try:
                    c_elegans_seq = accessions_and_sequences[c_elegans_accession_number]
                except:
                    c_elegans_seq = BioPython.get_aa_seq_from_entrez(c_elegans_accession_number)
                print("C.elegans seq: " + c_elegans_seq)

                variants = human_genes_and_variants[human_gene_name]
                print("variants:", ", ".join(variants))
                print("executor line 494")
                for variant in variants:
                    print("For variant: " + variant)
                    former_aa, place, current_aa = Strings.from_variant_string_to_tuple(variant)
                    result, count, c_elegans_location, alignment_conservation_score = \
                        BioPython.pairwise_alignment_inspector(human_seq,
                                                               c_elegans_seq,
                                                               Strings.from_name_to_symbol(former_aa),
                                                               place)
                    if result.startswith("not conserved") or result.startswith("not in the human"):
                        # filter if amino acid is not conserved
                        continue
                    print("For gene:", human_gene_name, "with variant", variant + ", the amino acid is", result)
                    line = data[human_gene_name] + "\t" + variant + "\t" + result + "\t" + str(
                        c_elegans_location) + "\t" + str(alignment_conservation_score) + "\t" + str(count) + "\t" + \
                           mmp_data + "\n"
                    FileMaker().write_to(file_type, line, f)
        if file_type != FileType.CONSOLE:
            f.close()

    # receives (1) file type to know whether to print to a file or to console, human genes and variants dictionary, read
    # all variants and extracts for each variant the sequence of its human gene and its C.elegans ortholog, runs
    # pairwise alignment and returns data regarding the alignments of the two sequences.
    def get_variants_data_for_server(self, human_genes_names_and_variants, sources_bar=2, length_range=(0.5, 2)):
        print("def: get_variants_data_for_server")
        if not executor.is_connected():
            return None, "No Internet Connection"

        failed_genes = []
        output = []
        bp = BioPython()
        print("Human genes:", list(human_genes_names_and_variants.keys()))
        for human_gene_name in human_genes_names_and_variants:
            human_seq = bp.get_human_aa_seq_by_name(human_gene_name)
            if not human_seq:
                failed_genes.append([human_gene_name, "Couldn't find gene sequence"])
                continue
            mmp_data = self.mmp_data_by_gene_name[human_gene_name] if human_gene_name in self.mmp_data_by_gene_name \
                else "No mention in MMP"

            results, _  = executor.find_me_orthologs_for_human(human_genes=[human_gene_name],
                                                              sources_bar=sources_bar,
                                                              length_range=length_range)
            true_matches_pairs, false_matches = results
            for gene_list in false_matches:
                failed_genes.append([gene_list[0], "No orthologs were found"])
            orthologs_names = []
            for lst in true_matches_pairs:
                orthologs_names.append(lst[1])

            orthologs_id_WB = Ensembl.get_c_elegans_genes_ids_by_genes_names(orthologs_names)
            print("Human gene name:", human_gene_name, ", seq:", human_seq, ", orthologs gene ids WB:",
                   orthologs_id_WB)

            for ortholog_id_WB in orthologs_id_WB:
                print("for ortholog:", ortholog_id_WB)
                c_elegans_seq = BioPython().get_c_elegans_aa_seq_by_id(ortholog_id_WB)
                if not c_elegans_seq:
                    continue
                print("C.elegans seq:", c_elegans_seq)

                variants = human_genes_names_and_variants[human_gene_name]
                print("variants: ", ", ".join(variants))
                print("executor line 556")
                for variant in variants:
                    print("For variant: " + variant)
                    former_aa, place = Strings.from_variant_string_to_tuple(variant)
                    print(former_aa , place)
                    print("inputs:", human_seq , c_elegans_seq , Strings.from_name_to_symbol(former_aa) , place)
                    result, count, c_elegans_location, alignment_conservation_score = \
                        BioPython.pairwise_alignment_inspector(human_seq,
                                                               c_elegans_seq,
                                                               Strings.from_name_to_symbol(former_aa),
                                                               place)
                    print("For gene:", human_gene_name, "with variant", variant + ", the amino acid is", result)
                    line = [human_gene_name, variant, Ensembl.get_gene_name_by_gene_id(ortholog_id_WB), result,
                            c_elegans_location, alignment_conservation_score, count, mmp_data]
                    output.append(line)

        return output, failed_genes

    # keys = human gene names
    @staticmethod
    def add_mmp_record_to_data(data_file_path, data_file_name, key_index, value_indexes, new_file_name):
        print("def: add_mmp_record_to_data")
        data = FileReader(data_file_path, data_file_name).readData(1, True)
        mmp_data_by_gene_name = FileReader(FileReader.research_path + r"\Data",
                                           r"\mmp_mut_strains_data_Mar14.txt"). \
            from_MMP_file_to_dict_with_listed_values(key_index, value_indexes, delete_first=True)
        f = open(new_file_name, FileMode.WRITE.value)
        gene_names = data.keys()
        for gene_name in gene_names:
            if gene_name in mmp_data_by_gene_name:
                f.write(data[gene_name] + "\t" + "".join(mmp_data_by_gene_name[gene_name]) + "\n")
            else:
                f.write(data[gene_name] + "\t" + "Not mentioned in MMP" + "\n")
        f.close()

    @staticmethod
    def add_failed_genes(failed_genes_dic, original_genes_list, current_genes, reason):
        print("def add_failed_genes")
        #print("function inputs: " , failed_genes_dic, original_genes_list, current_genes, reason)
        for gene_id in original_genes_list:
            print(gene_id)
            gene_name = Ensembl.get_gene_name_by_gene_id(gene_id)
            #if gene_name is None:
                #gene_name = DataExtracter.get_c_elegans_gene_name_for_gene_id(gene_id)
            print(gene_name)
            if gene_id not in current_genes and gene_name not in current_genes:
                if gene_name not in failed_genes_dic and gene_id not in failed_genes_dic:
                    if gene_name:
                        failed_genes_dic[gene_name] = reason
                    else:  # gene id didn't have gene name
                        failed_genes_dic[gene_id] = reason

    @staticmethod
    def find_me_orthologs_for_human(human_genes,
                                    genes_in_names: bool = True,
                                    sources_bar: int = 2,
                                    length_range: tuple = (0.5, 2),
                                    domains_range: tuple = (0.5, 3),
                                    species="Human"):
        print("def: find_me_orthologs_for_human ")
        true_matches = {}
        failed_genes = {}
        if not executor.is_connected():
            return None, "No Internet Connection"
        print(human_genes, genes_in_names, species, failed_genes)
        human_genes_ids, error = DataExtracter.get_genes_ids(human_genes, genes_in_names, species, failed_genes)
        print(human_genes_ids,error)
        if not human_genes_ids:
            print(error)
            return executor.get_result_list(true_matches, failed_genes), error

        # from human gene id to C.elegans gene id dictionary
        orthologs_dic = FileReader(FileReader.research_path + r"\Data",
                                   r"\ortholist_master",
                                   FileType.TSV).from_file_to_dict_with_plural_values(4, 0, True)
        relevant_orthologs_dic = DataExtracter.get_filtered_dic_of_orthologs(human_genes_ids, orthologs_dic)
        executor.add_failed_genes(failed_genes, human_genes_ids, relevant_orthologs_dic, "Couldn't find orthologs")
        if DataExtracter.is_dict_empty(relevant_orthologs_dic):
            return executor.get_result_list(true_matches, failed_genes), "No orthologs found"
        print("Orthologs: ", relevant_orthologs_dic)

        # now we have the ortholist-orthologs to our human genes.
        # next step - filtration by sources
        # sources dic is a dictionary with a tuple-keys of (human gene id, c.elegans gene id) and values of number of
        # sources supporting the pair is homologous
        sources_dic = FileReader(FileReader.research_path + r"\Data",
                                 r"\ortholist_master",
                                 FileType.TSV).from_file_to_dict_with_tuple_key(4, 0, 6, True)
        filtered_by_sources_orthologs = DataExtracter.filter_dic_by_sources(relevant_orthologs_dic,
                                                                            sources_dic,
                                                                            sources_bar)
        executor.add_failed_genes(failed_genes, human_genes_ids, filtered_by_sources_orthologs,
                                  "Failed at sources filtration")
        if DataExtracter.is_dict_empty(filtered_by_sources_orthologs, "after sources filtration"):
            return executor.get_result_list(true_matches, failed_genes), "No orthologs left after sources filtration"
        print(len(filtered_by_sources_orthologs), "genes are left:", filtered_by_sources_orthologs,
              "after filtration by sources")

        # now orthologs are filtered by sources
        # next step - filtration by length ratio
        print("line 1")
        filtered_by_length_orthologs = DataExtracter.filter_genes_by_length_differences(filtered_by_sources_orthologs,
                                                                                        length_range)
        print("line 2!")
        executor.add_failed_genes(failed_genes, human_genes_ids, filtered_by_length_orthologs, "Failed at length filtration")
        print("line 3!")
        if DataExtracter.is_dict_empty(filtered_by_length_orthologs, "after length filtration"):
            print("we are inside the if!!")
            print(executor.get_result_list(true_matches, failed_genes), "No orthologs left after length filtration")
            return executor.get_result_list(true_matches, failed_genes), "No orthologs left after length filtration"
        print(len(filtered_by_length_orthologs), "genes are left:", filtered_by_length_orthologs,
              "after filtration by length")

        # now we have genes filtered by size and sources.
        # next step: by domains ratio
        filtered_by_conserved_domains = DataExtracter().filter_by_conserved_domains(filtered_by_length_orthologs,
                                                                                    domains_range,
                                                                                    key_organism="Human")
        executor.add_failed_genes(failed_genes, human_genes_ids, filtered_by_conserved_domains, "Failed at domains filtration")
        if DataExtracter.is_dict_empty(filtered_by_conserved_domains, "after domains filtration"):
            return executor.get_result_list(true_matches, failed_genes), "No orthologs left after ratio of conserved domains filtration"
        print(len(filtered_by_conserved_domains), "genes are left", filtered_by_conserved_domains,
              "after filtration by domains")

        # now we have genes filtered by size, sources and domains ratio.
        # next step: by opposite blast
        true_matches, dict_before_blast = executor.pair_pipeline(dic_of_optional_orthologs=filtered_by_conserved_domains)

        print("true matches:\n", true_matches)
        executor.add_failed_genes(failed_genes, human_genes_ids, [tup[0] for tup in true_matches], "Failed at reversed BLAST")
        return executor.get_result_list(true_matches, failed_genes), None

    #this function is built to answer Ronen's list of genes:
    #it is blackened because it uses pandas, and this package raises errors in A2HOSTING
    # @staticmethod
    # def check_if_gene_has_ortholog_for_worm(file_path, file_name, sheet_name='kinase'):
    #      print("def: check_if_gene_has_ortholog_for_worm")
    #      fd = FileReader(file_path, file_name)
    #      genes = fd.get_list_from_excel_using_pandas('WormBase Gene ID', sheet_name)
    #      print(genes)
    #      #genes_and_ortholog_data = fd.get_dictionary_from_excel_using_pandas('Public Name', 'Human Ortholog',
    #      #                                                                    sheet_name=sheet_name)
    #      result_dictionary, _  = executor.find_me_orthologs_for_worm(genes, False, sources_bar=2)
    #      true_results, failed_dict = result_dictionary
    #      #print("before_blast_lst:" , before_blast_dict)
    #      print("true_results: ", true_results)
    #      print("failed result: ", failed_dict)
    #      df = pd.DataFrame(true_results, columns=["Your Gene","Ortholog","# Sources","Conserved Domains Ratio(C.elegans\Human)"
    #                                               ,"Human Gene Length",	"C.elegans Gene Length","Ortholog Description"])
    #      print(df)
         # df = df[["Your Gene","Ortholog"]]
         # df1 = df.astype(str).groupby('Your Gene').agg(','.join).reset_index()
         # df1['number_of_results_after_blast'] = df1['Ortholog'].apply(lambda x: x.count(',')+1 if "," in x else 1)
         # print(df1)
         # for failed_gene in failed_dict:
         #     print(failed_gene)
         #     if failed_gene[1] == 'Failed at reversed BLAST':
         #         print("IN")
         #         failed_gene_new = failed_gene.extend(["0"])
         #         df1.loc[len(df1)] = failed_gene_new
         #         failed_dict.remove(failed_gene)
         # print(failed_dict)
         # df1 = df1.sort_values(['Your Gene'])
         # print(df1)
         # #df1['number_of_results_before_blast'] = before_blast_lst
         # print(df1)
         # for failed_gene in failed_dict:
         #     failed_gene.extend([""])
         #     df1.loc[len(df1)] = failed_gene
         # print(df1)
         # gene_dict ={}
         # for key, value in before_blast_dict.items():
         #    for name in value:
         #        if name in gene_dict:
         #            gene_dict[name] +=1
         #        else:
         #            gene_dict[name] = 1
         # print(gene_dict)
         #df.to_csv("all_togther_49.csv")

         # parameters =
         # orthologous_genes = [lst[0] for lst in true_results]
         # if not orthologous_genes:
         #     print("Something went wrong, orthologous genes is empty")
         #     return None
         # print("orthologous_genes:", orthologous_genes)
         # count = 0
         # for worm_gene_name in genes_and_ortholog_data:
         #     res = 1 if worm_gene_name in orthologous_genes else 0
         #     count += res
         #     print(worm_gene_name, "input:", genes_and_ortholog_data[worm_gene_name], "output:", res)
         # print(count, "had orthologs out of", len(genes_and_ortholog_data))


    #this function is built to answer Ronen's list of genes:
    #it is blackened because it uses pandas, and this package raises errors in A2HOSTING
    #@staticmethod
    # def check_if_gene_has_ortholog_for_human(file_path, file_name, sheet_name='kinase'):
    #      print("def: check_if_gene_has_ortholog_for_human")
    #      fd = FileReader(file_path, file_name)
    #      genes = fd.get_list_from_excel_using_pandas('Ensembl ID', sheet_name)
    #      print(genes)
    #      #genes_and_ortholog_data = fd.get_dictionary_from_excel_using_pandas('Public Name', 'Human Ortholog',
    #      #                                                                    sheet_name=sheet_name)
    #      result_dictionary, _ , before_blast_dict = executor.find_me_orthologs_for_human(genes, False, sources_bar=2)
    #      true_results, failed_dict = result_dictionary
    #      print("before_blast_lst:" , before_blast_dict)
    #      print("true_results: ", true_results)
    #      print("failed result: ", failed_dict)
    #      df = pd.DataFrame(true_results, columns=["Your Gene","Ortholog","# Sources","Conserved Domains Ratio(C.elegans\Human)"
    #                                               ,"Human Gene Length",	"C.elegans Gene Length","C.elegans Gene Description"])
    #      print(df)
         #########################################################
         # those lines are for checking whether the blast step is changing the resultes
         #########################################################

         # df = df[["Your Gene","Ortholog"]]
         # df1 = df.astype(str).groupby('Your Gene').agg(','.join).reset_index()
         # df1['number_of_results_after_blast'] = df1['Ortholog'].apply(lambda x: x.count(',')+1 if "," in x else 1)
         # print(df1)
         # failed_dict_copy = []
         # for failed_gene in failed_dict:
         #     print(failed_gene)
         #     if failed_gene[1] == 'Failed at reversed BLAST':
         #         print("IN: Failed at reversed BLAST")
         #         failed_gene.extend(["0"])
         #         print(failed_gene)
         #         df1.loc[len(df1)] = failed_gene
         #         #failed_dict.remove(failed_gene)
         #     else:
         #         failed_dict_copy.append(failed_gene)
         # print(failed_dict_copy)
         # print(df1)
         # df1 = df1.sort_values(['Your Gene'])
         # print(df1)
         # #df1['number_of_results_before_blast'] = before_blast_lst
         # print(df1)
         # for failed_gene in failed_dict_copy:
         #     print(failed_gene)
         #     print(failed_dict_copy)
         #     failed_gene.extend([""])
         #     df1.loc[len(df1)] = failed_gene
         # print(df1)
         # gene_dict ={}
         # for key, value in before_blast_dict.items():
         #    for name in value:
         #        if name in gene_dict:
         #            gene_dict[name] +=1
         #        else:
         #            gene_dict[name] = 1
         # print(gene_dict)
         # gene_dict1 = {k: v for k, v in sorted(gene_dict.items(), key=lambda item: item[0])}
         # print(gene_dict1)
         # lst_before = list(gene_dict1.values())
         # print(lst_before)
         # for i in range(len(failed_dict_copy)):
         #    lst_before.append(0)
         # print(lst_before)
         # #df1.to_csv("results_of_set_of_confirmed_orthologs8.csv")
         # df1['number_of_results_before_blast'] = lst_before
         # print(df1)
         # df1.to_csv("results_of_set_of_confirmed_orthologs4.csv")
         ###########################################################

         # parameters =
         # orthologous_genes = [lst[0] for lst in true_results]
         # if not orthologous_genes:
         #     print("Something went wrong, orthologous genes is empty")
         #     return None
         # print("orthologous_genes:", orthologous_genes)
         # count = 0
         # for worm_gene_name in genes_and_ortholog_data:
         #     res = 1 if worm_gene_name in orthologous_genes else 0
         #     count += res
         #     print(worm_gene_name, "input:", genes_and_ortholog_data[worm_gene_name], "output:", res)
         # print(count, "had orthologs out of", len(genes_and_ortholog_data))


    # pipeline that provides you with humans genes that are orthologous for your worm ones
    @staticmethod
    def find_me_orthologs_for_worm(worm_genes,
                                   genes_in_names: bool = True,
                                   sources_bar: int = 2,
                                   length_range: tuple = (0.5, 2),
                                   domains_range: tuple = (0.5, 2),
                                   species="C.elegans"):
        print("def: find_me_orthologs_for_worm")
        true_matches = {}
        failed_genes = {}
        if not executor.is_connected():
            return None, "No Internet Connection"
        worm_genes_ids, error = DataExtracter.get_genes_ids(worm_genes, genes_in_names, species, failed_genes)
        print("worm_genes_ids", worm_genes_ids)
        if not worm_genes_ids:
            print("error" , error)
            return executor.get_result_list(true_matches, failed_genes), error

        # from C.elegans gene id to human gene id dictionary
        orthologs_dic = FileReader(FileReader.research_path + r"\Data",
                                   r"\ortholist_master",
                                   FileType.TSV).from_file_to_dict_with_plural_values(0, 4, True)
        ## this dict (orthologs_dic) contains all the worm ids genes and their homologs from the ortholist

        relevant_orthologs_dic = DataExtracter.get_filtered_dic_of_orthologs(worm_genes_ids, orthologs_dic)
        print("relevant_orthologs_dic" , relevant_orthologs_dic)     ## only the genes that are relevant for our lookup

        executor.add_failed_genes(failed_genes, worm_genes_ids, relevant_orthologs_dic, "Couldn't find orthologs")

        if DataExtracter.is_dict_empty(relevant_orthologs_dic):
            return executor.get_result_list(true_matches, failed_genes), "No orthologs found"
        print("Orthologs: ", relevant_orthologs_dic)

        # now we have the ortholist-orthologs to our human genes.
        # next step - filtration by sources
        # sources dic is a dictionary with a tuple-keys of (human gene id, c.elegans gene id) and values of number of
        # sources supporting the pair is homologous
        sources_dic = FileReader(FileReader.research_path + r"\Data",
                                 r"\ortholist_master",
                                 FileType.TSV).from_file_to_dict_with_tuple_key(0, 4, 6, True)
        filtered_by_sources_orthologs = DataExtracter.filter_dic_by_sources(relevant_orthologs_dic,
                                                                            sources_dic,
                                                                            sources_bar)
        executor.add_failed_genes(failed_genes, worm_genes_ids, filtered_by_sources_orthologs, "Failed at sources filtration")
        if DataExtracter.is_dict_empty(filtered_by_sources_orthologs, "after sources filtration"):
            return executor.get_result_list(true_matches, failed_genes), "No orthologs left after sources filtration"
        print(len(filtered_by_sources_orthologs), "genes are left", filtered_by_sources_orthologs,
              "after filtration by sources")

        # now orthologs are filtered by sources
        # next step - filtration by length ratio
        filtered_by_length_orthologs = DataExtracter.filter_genes_by_length_differences(filtered_by_sources_orthologs,
                                                                                        length_range,
                                                                                        key_gene="c_elegans")
        executor.add_failed_genes(failed_genes, worm_genes_ids, filtered_by_length_orthologs, "Failed at length filtration")
        if DataExtracter.is_dict_empty(filtered_by_length_orthologs, "after length filtration"):
            return executor.get_result_list(true_matches, failed_genes), "No orthologs left after length filtration"
        print(len(filtered_by_length_orthologs), "genes are left", filtered_by_length_orthologs,
              "after filtration by length")

        # now we have genes filtered by size and sources.
        # next step: by domains ratio
        filtered_by_conserved_domains = DataExtracter().filter_by_conserved_domains(filtered_by_length_orthologs,
                                                                                    domains_range,
                                                                                    key_organism="C.elegans")
        executor.add_failed_genes(failed_genes, worm_genes_ids, filtered_by_conserved_domains, "Failed at domains filtration")
        if DataExtracter.is_dict_empty(filtered_by_conserved_domains, "after domains filtration"):
            return executor.get_result_list(true_matches, failed_genes), "No orthologs left after ratio of conserved domains filtration"
        print(len(filtered_by_conserved_domains), "genes are left", filtered_by_conserved_domains,
              "after filtration by domains")
        #dict_before_blast = {}
        #for key, value in filtered_by_conserved_domains.items():
            #print (key , len(value))
            #dict_before_blast[key] = len(value)
            #list_before_blast.append(len(value))
        print("len before blast step: ",len(filtered_by_conserved_domains))

        # now we have genes filtered by size, sources and domains ratio.
        # next step: by opposite blast
        true_matches , dict_before_blast = executor.pair_pipeline(dic_of_optional_orthologs=filtered_by_conserved_domains,
                                              key_species="Human")
        print(true_matches)
        executor.add_failed_genes(failed_genes, worm_genes_ids, [tup[0] for tup in true_matches], "Failed at reversed BLAST")

        return executor.get_result_list(true_matches, failed_genes), None

    # parses input for the flask server's variants function
    @staticmethod
    def parse_input(input_genes):
        print("def: parse_input")
        # input_genes = "TCP1:[Asn284Ser,Ala453Glu],DAP3:[Leu138Phe,Glu369Lys]"
        dic = {}
        while input_genes:
            gene = input_genes[:input_genes.find(":")]
            input_genes = input_genes[len(gene)+1:]
            variations_str = input_genes[input_genes.find("[")+1:input_genes.find("]")]
            variations = variations_str.split(",")
            dic[gene] = variations
            input_genes = input_genes[len(variations_str)+3:]
        return dic

    # parses output for the flask's server finding functions (works only for dictionaries)
    @staticmethod
    def dictionary_output_parser(results):
        print("def: dictionary_output_parser")
        if isinstance(results, dict):
            output = ''
            for key in results:
                if isinstance(key, tuple):
                    output += ", ".join(key)
                else:
                    output += key
                if results.get(key) and results.get(key) != [""]:
                    output += ": " + ", ".join([str(item) for item in results.get(key)]) + "\n"
            return output
        return results

    @staticmethod
    def is_connected():
        print("def: is_connected")
        try:
            # connect to the host -- tells us if the host is actually reachable
            socket.create_connection(("www.google.com", 80))
            return True
        except OSError:
            pass
        print("No Internet Connection")
        return False

# irrelevant function
    # def get_shinjini_data(self, human_genes_names, c_elegans_genes_names):
    #     de = DataExtracter()
    #     for i in range(len(human_genes_names)):
    #         human_gene_name = human_genes_names[i]
    #         c_elegans_gene_name = c_elegans_genes_names[i]
    #         conserved_domains_value = de.get_conserved_domains_ratio_of_pair(c_elegans_gene_name,
    #                                                                          human_gene_name)
    #         human_seq = BioPython().get_human_protein_from_uniprot_by_gene_name(human_gene_name)
    #         try:
    #             c_elegans_id_wb = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
    #         except:
    #             print("Couldn't find id for", c_elegans_gene_name)
    #             continue
    #         c_elegans_seq = BioPython().\
    #             get_c_elegans_aa_seq(c_elegans_id_wb)
    #         conservation_score = BioPython.get_conservation_score(human_seq, c_elegans_seq)
    #         print(human_gene_name, c_elegans_gene_name, conservation_score, conserved_domains_value)


