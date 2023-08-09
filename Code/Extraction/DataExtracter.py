from Code.Enum.FileMode import FileMode
from Code.Enum.FileType import FileType
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.Ensembl import Ensembl
from Code.Utils.BioPython import BioPython
from Bio import Entrez


class DataExtracter:

    def __init__(self):
        self.accession_numbers_to_hit_ids = FileReader(FileReader.research_path + r"\Executors",
                                                       r"\accession numbers and hit ids").\
                                                       from_file_to_dict_with_lists_as_values(0, 1)
        self.hit_ids_to_gene_names = FileReader(FileReader.research_path + r"\Executors",
                                                r"\blast-hit-ids-and-human-gene-names0-30",
                                                FileType.TSV).from_file_to_dict(0, 1)

        self.c_elegans_domains_dic = FileReader(FileReader.research_path + r"\Data",
                                                      r"\c-elegans-genes-and-conserved-domains-230619",
                                                FileType.TSV).from_file_to_dict(0, 1)

        self.human_domains_dic = FileReader(FileReader.research_path + r"\Data",
                                            r"\human-genes-and-conserved-domains-230619",
                                            FileType.TSV).from_file_to_dict(0, 1)

        self.orthologs = FileReader(FileReader.research_path + r"\Data",
                         r"\ortholist_master",
                         FileType.TSV).from_file_to_dict_with_plural_values(0, 4, True)

        self.sources = FileReader(FileReader.research_path + r"\Data",
                                  r"\ortholist_master",
                                  FileType.TSV).from_file_to_dict_with_tuple_key(0, 4, 6, True)

        self.humans_id_cd_length = FileReader(FileReader.research_path + r"\Data",
                                              r"\human_cd_length.txt").get_genes_cd_length(1, 2, True)

        self.c_elegans_cd_length = FileReader(FileReader.research_path + r"\Data",
                                              r"\c_elegans_genes_cd_length.txt").get_genes_cd_length(1, 2, True)


    # receives (1) a dictionary of c.elegans genes id (number) as keys and accession numbers as values, enact the
    # function that chooses the longest isoform of all accession number, and returns a dictionary with c.elegans gene
    # ids and the chosen longest accession number
    @staticmethod
    def from_multiple_accessions_to_one(genes_and_accessions: dict):
        print("def: from_multiple_accessions_to_one")
        genes_and_chosen_accessions = {}
        keys = genes_and_accessions.keys()
        count = 0
        for gene_id in keys:
            count += 1
            accessions_set = genes_and_accessions[gene_id]
            chosen_isoform = BioPython.get_longest_accession(accessions=accessions_set)
            genes_and_chosen_accessions[gene_id] = chosen_isoform
            if count % 100 == 0:
                print("got to", str(count), "with", gene_id)
        return genes_and_chosen_accessions

    # receives (1) list of human/c.elegans genes id (WB\ENSG) and (2) a term to search the needed element, make an
    # http request to NCBI to extract the number of domains for the gene, and returns a dictionary of the gene id
    #  and its number of domains
    @staticmethod
    def get_conserved_domains(genes: list):
        print("def: get_conserved_domains")
        genes_and_conserved_domain = {}
        parsed_number = DataExtracter().get_conserved_domain_per_gene_id(genes)
        for gene in genes:
            if parsed_number is not None:
                genes_and_conserved_domain[gene] = parsed_number
        return genes_and_conserved_domain

    # receives (1) a gene id(WB/ENSG) and (2),(3) terms to extract only the needed number out of the line, makes an
    # http request to extract the domains of each gene, parses it to an integer and returns it
    @staticmethod
    def get_conserved_domain_per_gene_id(gene_id):
        print("def: get_conserved_domain_per_gene_id")
        hr = HttpRequester("https://www.ncbi.nlm.nih.gov/gene/?term=" + gene_id)
        data = hr.make_request()
        number = DataExtracter.extract_number_from_data(data)
        try:
            parsed_number = int(number)
            return parsed_number
        except Exception as e:
            print("Exception in get_conserved_domain_per_gene_id:", e)
            print("for gene id:", gene_id, "letter extracted is", number)
            ncbi_id = Ensembl.get_ncbi_id_by_gene_id(gene_id)
            if ncbi_id:
                return DataExtracter().get_conserved_domains_per_ncbi_id(ncbi_id)
        return None

    # extracts number of conserved domains from data given by http request
    @staticmethod
    def extract_number_from_data(data, start_term: str = "Conserved Domains (", end_term: str = ")"):
        print("def: extract_number_from_data")
        if isinstance(data, str):  # request was successful
            index = data.find(start_term)
            shorter_data = data[index + len(start_term):]
            print("shorter_data")
            #print(shorter_data)
            end_index = shorter_data.find(end_term)
            number = shorter_data[:end_index]
            print("the number is: " , number)
            return number
        return None

    @staticmethod
    def get_conserved_domains_per_ncbi_id(ncbi_id):
        print("def: get_conserved_domains_per_ncbi_id")
        hr = HttpRequester("https://www.ncbi.nlm.nih.gov/gene/" + ncbi_id)
        data = hr.make_request()
        print("lets see why conserved domain is empty")
        #print(data)
        number = DataExtracter.extract_number_from_data(data)
        print("the number is: " ,number)
        try:
            parsed_number = int(number)
            return parsed_number
        except Exception as e:
            print("Exception in get_conserved_domains_per_ncbi_id:", e)
            print("for ncbi gene id:", ncbi_id, " letter extracted is", number)
        return None

    # receives a worm gene name and a human gene name, get conserved domains for both genes, and returns their ratio
    def get_conserved_domains_ratio_of_pair(self, c_elegans_gene_name, human_gene_name):
        print("def: get_conserved_domains_ratio_of_pair")
        c_elegans_id = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
        if not c_elegans_id:
            print("no id for gene:", c_elegans_gene_name)
            return None
        human_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
        if not human_id:
            print("no id for gene:", human_gene_name)
            return None
        return self.get_conserved_domains_ratio_of_pair_by_id(c_elegans_id, human_id)

    def get_conserved_domains_ratio_of_pair_by_id(self, c_elegans_gene_id, human_gene_id):
        print("def: get_conserved_domains_ratio_of_pair_by_id")
        c_elegans_domains = self.get_conserved_domain_per_gene_id(c_elegans_gene_id)
        human_domains = self.get_conserved_domain_per_gene_id(human_gene_id)
        try:
            return float(c_elegans_domains) / float(human_domains)
        except Exception as e:
            print("Exception in get_conserved_domains_ratio_of_pair:", e)
            return None

    # receives worm gene name and human gene name, reads the ortholist file to a dictionary with
    # (worm gene id, human gene id) keys and number of sources as values, and returns the number of sources if the tuple
    # is in the dictionary, otherwise returns -1
    def get_sources(self, c_elegans_gene_name, human_gene_name):
        print("def: get_sources")
        c_elegans_gene_id = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
        human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
        if c_elegans_gene_id and human_gene_id:
            if c_elegans_gene_id in self.orthologs:
                if human_gene_id in self.orthologs[c_elegans_gene_id]:
                    print(int(self.sources[(c_elegans_gene_id, human_gene_id)]))
                    return int(self.sources[(c_elegans_gene_id, human_gene_id)])
            print("Couldn't find number of sources for tuple:", c_elegans_gene_name, human_gene_name)
        return -1

    # receives a (1) dictionary with gene name as key and the gene id as value, and (2) a list of gene names, and
    # returns a list of the genes ids
    @staticmethod
    def from_human_genes_names_to_human_genes_id(genes_names: list):
        print("def: from_human_genes_names_to_human_genes_id")
        genes_ids = []
        for gene_name in genes_names:
            gene_id = Ensembl.get_human_gene_id_by_gene_name(gene_name)
            if gene_id:
                genes_ids.append(gene_id)
        return genes_ids

    # receives (1) a dictionary of c.elegans genes and number of domains, (2) a dictionary of human genes ids and
    # number of domains, (3) a dictionary with human genes ids as keys and c.elegans orthologs as values, and
    # returns a new dic with the human genes and their orthologs if they are inside the range of ratio
    def filter_by_conserved_domains(self,
                                    orthologs_dic: dict,
                                    domains_range,
                                    key_organism):
        print("def: filter_by_conserved_domains")
        key_domains_dic = self.human_domains_dic if key_organism == "Human" else self.c_elegans_domains_dic
        orthologs_domains_dic = self.c_elegans_domains_dic if key_organism == "Human" else self.human_domains_dic

        new_orthologs_dic = {}
        for key_gene_id in orthologs_dic:
            orthologs = orthologs_dic[key_gene_id]
            if key_gene_id in key_domains_dic:
                key_domains = key_domains_dic[key_gene_id]
            else:
                key_domains = DataExtracter().get_conserved_domain_per_gene_id(key_gene_id)
            for ortholog in orthologs:
                try:
                    ortholog_domains = orthologs_domains_dic[ortholog]
                except Exception as e:
                    print("Exception in filter_by_conserved_domains:", e)
                    ortholog_domains = DataExtracter().get_conserved_domain_per_gene_id(ortholog)

                try:
                    ratio = float(ortholog_domains) / float(key_domains)
                    if domains_range[0] <= ratio <= domains_range[1]:
                        print("Ratio between", ortholog, "and", key_gene_id, "is", ratio,
                              ", thus passes the domains ratio bar")
                        DataExtracter.add_to_dictionary(new_orthologs_dic, key_gene_id, ortholog)
                    else:
                        print("Ratio between", key_gene_id, "and", ortholog + " is", ratio,
                              ", thus not good enough")
                except Exception as e:
                    print("Exception in filter_by_conserved_domains:", e)
                    print("Problem occurred with dividing", ortholog_domains, "by", key_domains)
                    DataExtracter.add_to_dictionary(new_orthologs_dic, key_gene_id, ortholog)
        return new_orthologs_dic

    # receives (1) a dictionary of c.elegans genes and number of domains, (2) a dictionary of human genes ids and
    # number of domains, (3) a dictionary with human genes names as keys and c.elegans orthologs as values, and (4)
    # a dictionary of human genes names as keys and human genes ids as values, and adds to the orthologs the ratio
    # of domains between C.elegans and humans
    @staticmethod
    def add_conserved_domains_info(c_elegans_domains_dic: dict, human_domains_dic: dict, orthologs_dic: dict,
                                   human_name_id_dic: dict):
        print("def: add_conserved_domains_info")
        new_orthologs_dic = {}
        human_genes = orthologs_dic.keys()
        for human_gene_name in human_genes:
            try:
                human_gene_id = human_name_id_dic[human_gene_name]
            except KeyError:
                print(human_gene_name + " has no equivalent gene id in the dictionary")
                continue
            try:
                c_elegans_orthologs = FileReader.fromStringToTuplesList(orthologs_dic[human_gene_name])
            except KeyError:
                print(human_gene_name + "doesn't exist in the  file of orthologs")
                continue
            if human_gene_id in human_domains_dic:
                human_domains = human_domains_dic[human_gene_id]
                orthologs_with_domains = DataExtracter.add_domains_similarity_score(int(human_domains),
                                                                                    c_elegans_domains_dic,
                                                                                    c_elegans_orthologs)
                new_orthologs_dic[human_gene_id] = orthologs_with_domains
            else:
                print("human gene " + human_gene_id + " doesn't exist in the domains dictionary")
                pass
        return new_orthologs_dic

    # receives (1) number of human gene domains, (2) a dictionary of c.elegans genes ids (WB?) and number of domains as
    # values and (3) list of C.elegans orthologs, and adds the ratio of number of domains between each C.elegans
    # ortholog gene and the human gene domains to the list of orthologs
    @staticmethod
    def add_domains_similarity_score(human_domains: int, c_elegans_domains_dic: dict, orthologs: list):
        print("def: add_domains_similarity_score")
        orthologs_with_domains = []
        for ortholog, sources in orthologs:
            if ortholog in c_elegans_domains_dic:
                ortholog_domains = int(c_elegans_domains_dic[ortholog])
                domainsScore = min(ortholog_domains, human_domains) / max(ortholog_domains, human_domains)
                orthologs_with_domains.append((ortholog, sources, domainsScore))
            else:
                print("c elegans gene " + ortholog + " doesn't exist in dictionary of c.elegans domains")
        return orthologs_with_domains

    # receives (1) the number of results files, (2) the file path, (3) the file name, and (4) an optional start file
    # index, goes through the files and extract for each accession number the given hit ids and their stats, and returns
    # adictionary of each accession number and their hit ids, and each hit id and it stats
    @staticmethod
    def get_hit_ids_for_accession_numbers(end_file: int, file_path, file_name, start_file: int = 0):
        print("def: get_hit_ids_for_accession_numbers")
        accession_numbers_and_hit_ids = {}
        hit_ids_and_hsp = {}
        fd = FileReader(file_path, file_name, FileType.TSV)
        accession_number = "default"
        hit_id = "default"
        for result_index in range(start_file, end_file):
            f = fd.read_results_file()
            print("created a file to extract the results!")
            line = f.readline()  # accession-number...
            while line:  # line is not an empty string, thus we haven't reached the end of the file yet
                if line.startswith("Accession Number"):
                    accession_number = line.rstrip()[len("Accession Number: "):]
                    print("now working on accession number: " + accession_number)
                    line = f.readline()
                elif line.startswith("****Alignment****"):
                    next_line = f.readline().rstrip()
                    hit_id = next_line[next_line.find("|") + 1: next_line.find("|", next_line.find("|") + 1)]
                    extra_letter = next_line[next_line.find(hit_id) + len(hit_id) + 1:
                    next_line.find(hit_id) + len(hit_id) + 2]
                    if extra_letter:
                        hit_id = hit_id + "_" + extra_letter
                    DataExtracter.add_to_dictionary(accession_numbers_and_hit_ids, accession_number, hit_id)
                    line = f.readline()
                elif line.startswith("****HSP****"):
                    e_value = f.readline().rstrip()[len("e value: "):]
                    score = f.readline().rstrip()[len("hsp score: "):]
                    DataExtracter.add_to_dictionary(hit_ids_and_hsp, hit_id, e_value + "+" + score)
                    line = f.readline()
                else:
                    line = f.readline()
            f.close()
        return accession_numbers_and_hit_ids, hit_ids_and_hsp

    # receives a (1) dictionary, (2) a key, and (3) a value, and added the value to a list of values for the
    # given key
    @staticmethod
    def add_to_dictionary(d: dict, key: str, value):
        print("def: add_to_dictionary")
        if key in d:
            if value not in d[key]:
                d[key].append(value)
        else:
            d[key] = [value]

    # receives (1) a dictionary with list of values for each key, and return a reversed dictionary
    @staticmethod
    def from_plural_valued_dict_to_plural_valued_reversed_dict(d: dict):
        print("def: from_plural_valued_dict_to_plural_valued_reversed_dict")
        reversed_dict = {}
        for key in d:
            value_list = d[key]
            for value in value_list:
                DataExtracter().add_to_dictionary(reversed_dict, value, key)
        return reversed_dict

    # receives (1) the data dictionary of orthologs and (2) a domain multiplying factor, creates a score for each
    # C.elegans gene and returns a dictionary of C.elegans genes and their orthology score
    @staticmethod
    def scoring_function(orthologs_dictionary: dict, domain_factor):
        print("def: scoring_function")
        c_elegans_genes_and_scores = {}
        human_genes = list(orthologs_dictionary.keys())
        for human_gene in human_genes:
            c_elegans_orthologs = FileReader.fromStringToTuplesList(orthologs_dictionary[human_gene])
            for homolog in c_elegans_orthologs:
                gene_name = homolog[0]
                sources = homolog[1]
                domain_equality = homolog[2]
                score = int(sources) + domain_factor * float(domain_equality)
                c_elegans_genes_and_scores[gene_name] = score

    @staticmethod
    def add_description_from_worm_base(data_file_path, new_data_file_path, column, empty_description, ortholog_column,
                                       start_term: str = "\"text\":\"", end_term: str = "\",\"evidence\""):
        print("def: add_description_from_worm_base")
        in_file = open(data_file_path)
        out_file = open(new_data_file_path, FileMode.WRITE.value)
        for line in in_file:
            lst = line.strip("\n").split("\t")
            description = lst[column]
            if description == empty_description:
                record = HttpRequester("http://rest.wormbase.org/rest/widget/gene/" +
                                       lst[ortholog_column] + "/overview").make_request()
                info_index = record.find("concise_description")
                shorter_info = record[info_index + len("concise_description"):]
                start_index = shorter_info.find(start_term)
                end_index = shorter_info.find(end_term)
                description = shorter_info[start_index + len(start_term): end_index]

                for i in range(column):
                    out_file.write(lst[i] + "\t")
                out_file.write(description + "\n")
            else:
                out_file.write(line)
        in_file.close()
        out_file.close()

    # receives (1) list of genes, (2) url address to the site where we can find phenotypes, and (3),(4) terms to extract
    # only the relevant phenotypes from the line given, and returns a dictionary of genes and phenotypes of variants
    @staticmethod
    def add_c_elegans_phenotypes(list_of_genes: list, url, search_term: str = "\"class\":\"phenotype\"",
                                 label_term: str = "\"label\":"):
        print("def: add_c_elegans_phenotypes")
        genes_and_phenotypes = {}
        for gene in list_of_genes:
            print("gene: " + gene)
            info = HttpRequester(url + gene + "/" + "phenotype").make_request()
            if info:
                search_index = info.find(search_term)
                while -1 < search_index < info.find("\"phenotype_not_observed\":") and \
                            search_index < info.find("\"phenotype_by_interaction\""):
                    info = info[search_index:]
                    label_index = info.find(label_term)
                    info = info[label_index:]
                    phenotype = info[len(label_term):info.find(",")].strip("\"")
                    print(phenotype)
                    if gene in genes_and_phenotypes:
                        genes_and_phenotypes[gene].append(phenotype)
                    else:
                        genes_and_phenotypes[gene] = [phenotype]
                    search_index = info.find(search_term)

        return genes_and_phenotypes

    @staticmethod
    def get_c_elegans_variants_phenotype(gene_id):
        print("def: get_c_elegans_variants_phenotype")
        url = "http://rest.wormbase.org/rest/widget/gene/"
        search_term: str = "\"class\":\"phenotype\""
        label_term: str = "\"label\":"
        info = HttpRequester(url + gene_id + "/" + "phenotype").make_request()
        if not info:
            return None
        search_index = info.find(search_term)
        while -1 < search_index < info.find("\"phenotype_not_observed\":") and \
                        search_index < info.find("\"phenotype_by_interaction\""):
            info = info[search_index:]
            label_index = info.find(label_term)
            info = info[label_index:]
            phenotype = info[len(label_term):info.find(",")].strip("\"")
        return phenotype

    # receives the C.elegans id and extracts the sequence (id -> ncbi id -> accession numbers -> accession number ->
    # seq), runs a blast search with the sequence and returns the human genes' hit ids
    def get_human_hit_ids(self, c_elegans_gene_id, c_elegans_gene_name):
        print("def: get_human_hit_ids")
        c_elegans_accession_number = None

        # gene id -> ncbi id -> accession number -> hit ids
        c_elegans_ncbi_id = BioPython.get_ncbi_id(c_elegans_gene_id)
        if c_elegans_ncbi_id:
            print("gene ncbi id is: ", c_elegans_ncbi_id)
            c_elegans_accession_number = BioPython().get_c_elegans_accession_number(c_elegans_ncbi_id)
            if c_elegans_accession_number:
                print("Gene's accession number is:", c_elegans_accession_number)
                hit_ids = self.accession_numbers_to_hit_ids.get(c_elegans_accession_number, None)
                if hit_ids:
                    return hit_ids

        # gene id -> amino acid seq -> blast -> hit ids
        c_elegans_seq = BioPython().get_c_elegans_aa_seq_by_id(c_elegans_gene_id)
        if c_elegans_seq:
            hit_ids = BioPython().pipeline_blast_with_seq("blastp", "nr", c_elegans_seq)
            if not hit_ids:
                print("Couldn't get hit ids for", c_elegans_gene_name)
                return None
            print("Length of hit ids list is:", len(hit_ids))
            return hit_ids
        return None

    # receives the human gene id and extracts the sequence, runs a blast search with the sequence and returns the human hit ids
    @staticmethod
    def get_c_elegans_hit_ids(human_gene_id, human_gene_name):
        print("def: get_c_elegans_hit_ids")
        seq = BioPython().get_human_aa_seq_by_id(human_gene_id)
        print("Human Seq:", seq)
        if not seq:
            return None
        print("reversed blast time...")
        hit_ids = BioPython().pipeline_blast_with_seq("blastp", "nr", seq, "Caenorhabditis elegans[organism]",
                                                      "Caenorhabditis elegans")
        if not hit_ids:
            print("Couldn't find hit ids for", human_gene_name)
            return None
        print("Length of hit ids list is: " + str(len(hit_ids)))
        return hit_ids

    # receives (1) blasted gene name, (2) ortholog gene name, (3) whether the blasted gene is
    # human or C.elegans, and runs a reverse blast on the genes to see if we get a match in the form of the ortholog
    # gene scheme: blasted gene id(WB) -> blasted gene id(number) -> accession numbers -> accession number -> hit ids
    # -> ortholog
    # genes names -> check if our ortholog gene is in there
    def check_reversed_blast_hit_ids(self, bp: BioPython, gene_name, ortholog_gene_name, gene_species="C.elegans"):
        print("def: check_reversed_blast_hit_ids")
        # if gene_species == "C.elegans":
        #     hit_ids = DataExtracter().get_human_hit_ids(gene_id, gene_name)
        # else:
        #     hit_ids = DataExtracter().get_c_elegans_hit_ids(gene_id, gene_name)
        # if not hit_ids:
        #     return False
        # for hit_id in hit_ids:

        hit_ids = bp.blast_results.get(gene_name, None)
        if not hit_ids:
            print("No hit ids were found for gene", gene_name)
            return None
        for hit_id in hit_ids:
            hit_gene_name = self.hit_ids_to_gene_names.get(hit_id, None)
            if hit_gene_name:
                if hit_gene_name == ortholog_gene_name:
                    print(gene_name, "and", ortholog_gene_name, "are orthologs indeed!")
                    return True
                else:
                    print("hit human gene name is", hit_gene_name, "thus there is no match")
            else:  # need to extract gene's name
                hit_gene_name = BioPython.get_gene_name_from_protein_accession([hit_id]).get(hit_id, None)
                print("The gene's name for hit id:", hit_id, "is:", hit_gene_name)
                if not hit_gene_name:
                    print("Couldn't find the human gene name to hit id: " + hit_id)
                elif hit_gene_name == ortholog_gene_name:
                    print(gene_name, "and", ortholog_gene_name, "are orthologs indeed!")
                    return True
        print("No match between", gene_name, "and", ortholog_gene_name)
        return False

    # receives (1) list of genes, (2) path of data file, (3) name of data file, (4) name for filtered file,
    # (5) phenotype column and (6) filtering word, and for each gene, if its name is in the genes list, and the
    # filtering word is in its phenotype description, it copies the info line to the new filtered file
    @staticmethod
    def filter_genes(short_list_genes, data_path, data_name, filtered_name, phenotype_column, filtering_word):
        print("def: filter_genes")
        print("So the function works...")
        f = open(data_path + data_name, FileMode.READ.value)
        filtered_file = open(filtered_name, FileMode.WRITE.value)
        f.readline()
        for line in f:
            info = line.strip("\n").split("\t")
            gene_name = info[0]
            if gene_name in short_list_genes:
                print(gene_name + " is in short list!")
                print(info[phenotype_column])
                if filtering_word in info[phenotype_column]:
                    filtered_file.write(line)
                else:
                    continue
            else:
                continue
        f.close()
        filtered_file.close()

    # receives (1) file path of unwanted genes, (2) said file name, (3) path of data, (4) name of data, (5) new data
    # file name, and copies only the data lines for human genes that are not included in the unwanted genes list, and
    # only if the matching C.elegans gene is not included in the unwanted genes for the human genes, to the new file
    @staticmethod
    def filter_genes_by_name(unwanted_genes_file_path, unwanted_genes_file_name, data_path, data_name, new_data_name):
        print("def: filter_genes_by_name")
        unwanted_genes_dic = FileReader(unwanted_genes_file_path,
                                           unwanted_genes_file_name).from_file_to_dict_with_plural_values(0, 2, False)
        data = open(data_path + data_name)
        new_data = open(new_data_name, FileMode.WRITE.value)
        data.readline()
        count = 0
        read = 0
        stayed = 0
        for line in data:
            read += 1
            lst = line.strip("\n").split("\t")
            if lst[0] not in unwanted_genes_dic:
                new_data.write("\t".join(lst) + "\n")
                stayed += 1
            else:
                if lst[2] not in unwanted_genes_dic[lst[0]]:
                    new_data.write("\t".join(lst) + "\n")
                    stayed += 1
                else:
                    print(line.strip("\n"))
                    count += 1
        data.close()
        new_data.close()
        print(count)
        print(read)
        print(stayed)

    # receives (1) list of keys, and (2) dictionary of keys and values and returns a new dictionary containing only the
    # keys and values for the keys given in input (1)
    @staticmethod
    def get_filtered_dic_of_orthologs(list_of_genes, gene_converting_dic):
        print("def : get_filtered_dic_of_orthologs")
        dic = {}
        print(list_of_genes)
        #print("gene_converting_dic ",gene_converting_dic)
        for gene in list_of_genes:
            print(gene)
            try:
                dic[gene] = gene_converting_dic[gene]
            except Exception as e:
                print("Exception in get_filtered_dic_of_orthologs:", e)
                print("Couldn't find orthologs for", gene)
        return dic

    # receives (1) dictionary of human-c.elegans orthologs, (2) dictionary containing tuples of human gene id and
    # c.elegans gene id as keys and number of sources supporting this is a homologous pair as values, and (3) the
    # efficient number of sources bar, and only the pairs that have sufficiently high number of sources are copied into
    # a new filtered dictionary
    @staticmethod
    def filter_dic_by_sources(orthologs_dic, sources_dic, bar):
        print("def: filter_dic_by_sources")
        print(orthologs_dic, bar)
        filtered_dic = {}
        for key in orthologs_dic:
            for ortholog in orthologs_dic[key]:
                try:
                    print(int(sources_dic[(key, ortholog)]))
                    if int(sources_dic[(key, ortholog)]) >= int(bar):
                        print(key, "with its ortholog", ortholog, "have passes the sources bar!")
                        DataExtracter.add_to_dictionary(filtered_dic, key, ortholog)
                except Exception as e:
                    print("Problem has occurred with:", e)
        return filtered_dic

    def get_pair_cd_length(self, human_gene_name, c_elegans_gene_name):
        print("def: get_pair_cd_length")
        human_gene_length = self.humans_id_cd_length.get(human_gene_name, None)
        if not human_gene_length:
            print("Gene", human_gene_name, "'s length wasn't found")
        else:
            human_gene_length = int(human_gene_length)
        c_elegans_gene_length = self.c_elegans_cd_length.get(c_elegans_gene_name, None)
        if not c_elegans_gene_length:
            print("Gene", c_elegans_gene_name, "'s length wasn't found")
        else:
            c_elegans_gene_length = int(c_elegans_gene_length)

        return human_gene_length, c_elegans_gene_length

    # receives (1) a dictionary containing genes name as keys, and orthologs as values, (2) a ratio bar, and (3) what
    # organism are the key genes from: human or C.elegans, after clearing out the length dicts for keys and values, it
    # deletes genes for which the length differences between the human gene and its ortholog
    # are too high
    @staticmethod
    def filter_genes_by_length_differences(d: dict, p: tuple = (0.5, 2), key_gene: str = "Human"):
        print("def: filter_genes_by_length_differences")
        if key_gene == "Human":
            key_id_cd_length = FileReader(FileReader.research_path + r"\Data",
                                          r"\human_cd_length.txt").get_genes_cd_length(0, 2, True)
            value_id_cd_length = FileReader(FileReader.research_path + r"\Data",
                                            r"\c_elegans_genes_cd_length.txt").get_genes_cd_length(0, 2, True)
        else:
            key_id_cd_length = FileReader(FileReader.research_path + r"\Data",
                                          r"\c_elegans_genes_cd_length.txt").get_genes_cd_length(0, 2, True)
            value_id_cd_length = FileReader(FileReader.research_path + r"\Data",
                                            r"\human_cd_length.txt").get_genes_cd_length(0, 2, True)
        #print("the output of filter_genes_by_length_differences is: " , key_id_cd_length,value_id_cd_length)
        genes_and_orthologs = {}
        for key in d:
            try:
                key_gene_length = float(key_id_cd_length[key])
            except Exception as e:
                print("Exception in filter_genes_by_length_differences:", e)
                print("gene:", key, "wasn't found in length dictionary")
                continue
            values = d[key]
            for value_gene in values:
                try:
                    value_gene_length = float(value_id_cd_length[value_gene])
                except Exception as e:
                    print("Exception in filter_genes_by_length_differences:", e)
                    print("gene:", value_gene, "wasn't found in length dictionary")
                    continue
                if (key_gene_length * p[0]) <= value_gene_length <= (key_gene_length * p[1]):
                    print("Gene:", key, "has length of", int(key_gene_length), "bp, and its ortholog", value_gene,
                          "has length of", int(value_gene_length), "bp, so they pass the length bar")
                    DataExtracter.add_to_dictionary(genes_and_orthologs, key, value_gene)
                else:  # genes have too different sizes
                    print("value gene's length is", value_gene_length, "and key gene's length is",
                          key_gene_length, "thus the value gene is", (value_gene_length * 100) / key_gene_length,
                            "percent of the key gene's length")
        print("the length's of the genes are: " ,int(key_gene_length) , int(value_gene_length))
        return genes_and_orthologs

    # receives (1) a possible-plural-valued-dictionary of worm-human orthologs and (2) a boolean value indicating if we
    # are interested in converting id to name or name to id, and returns a new worm-human dictionary with the converted keys and values
    @staticmethod
    def convert_dic(dic, from_id_to_name, key_organism: str = "C.elegans"):
        print("def: convert_dic")
        converted_dic = {}
        for key in dic:
            values = dic[key]
            for value in values:
                if from_id_to_name:
                    key_name = Ensembl.get_gene_name_by_gene_id(key)
                    #if key_organism == "C.elegans" and key_name is None:
                        #key_name = DataExtracter.get_c_elegans_gene_name_for_gene_id(key)
                    value_name = Ensembl.get_gene_name_by_gene_id(value)
                    #if key_organism != "C.elegans" and value_name is None:
                        #value_name = DataExtracter.get_c_elegans_gene_name_for_gene_id(value)
                    if key_name and value_name:
                        DataExtracter.add_to_dictionary(converted_dic, key_name, value_name)
                    else:
                        print("Couldn't find the gene name for either", key, "or/and", value)
                else:
                    if key_organism == "C.elegans":
                        key_id = Ensembl.get_c_elegans_gene_id_by_gene_name(key)
                        value_id = Ensembl.get_human_gene_id_by_gene_name(value)
                    else:
                        key_id = Ensembl.get_human_gene_id_by_gene_name(key)
                        value_id = Ensembl.get_c_elegans_gene_id_by_gene_name(value)
                    if key_id and value_id:
                        DataExtracter.add_to_dictionary(converted_dic, key_id, value_id)
                    else:
                        print("Couldn't find the gene id for either", key, "or/and", value)
        return converted_dic

    # receives (1) a dict and checks if it is empty. if so, it exits the function
    @staticmethod
    def is_dict_empty(dic: dict, step=""):
        print("def: is_dict_empty")
        if not dic:
            print("dictionary has no genes left", step)
            return True
        return False

    @staticmethod
    def get_genes_ids(list_of_genes, genes_in_names, species, failed_genes):
        print("def: get_genes_ids")
        error = None
        if genes_in_names:
            genes_ids, error = Ensembl.convert_from_names_to_ids(list_of_genes, species, failed_genes)
        else:  # genes in ids
            genes_ids = list_of_genes
        return genes_ids, error

    def get_status_tuple(self, c_elegans_gene_name, human_gene_name , key_species):
        print("def: get_status_tuple")
        print("key_species: ",key_species)
        if key_species == "C.elegans":
            return (self.get_sources(c_elegans_gene_name, human_gene_name),
                    self.get_conserved_domains_ratio_of_pair(c_elegans_gene_name, human_gene_name),
                    self.get_pair_cd_length(human_gene_name, c_elegans_gene_name),
                    c_elegans_gene_name + ": " + self.get_c_elegans_description_for_gene_id(
                    Ensembl.get_gene_id_by_gene_name(c_elegans_gene_name, "C.elegans")))
        else:
            print("lets print the output of the pair: ",c_elegans_gene_name, human_gene_name)
            print("1. : ", self.get_sources(c_elegans_gene_name, human_gene_name))
            print("2. : ", self.get_conserved_domains_ratio_of_pair(c_elegans_gene_name, human_gene_name))
            print("3. : ",self.get_pair_cd_length(human_gene_name, c_elegans_gene_name))
            print("4. : ", [BioPython.get_human_gene_id_from_ENTREZ_by_name(human_gene_name)])
            print("5. : ", BioPython.extract_Human_gene_description_from_ENTREZ(
                [BioPython.get_human_gene_id_from_ENTREZ_by_name(human_gene_name)]))

            return (self.get_sources(c_elegans_gene_name, human_gene_name),
                    self.get_conserved_domains_ratio_of_pair(c_elegans_gene_name, human_gene_name),
                    self.get_pair_cd_length(human_gene_name, c_elegans_gene_name),
                    human_gene_name + ": " + BioPython.extract_Human_gene_description_from_ENTREZ([BioPython.get_human_gene_id_from_ENTREZ_by_name(human_gene_name)]))



    ########### irrelevant functions ###########


    @staticmethod
    def fix_conserved_domain_info(data_path, data_name, c_elegans_domains_dic, human_domains_dic, domain_column,
                                  fixed_file, delete_first_line: bool = False):
        print("def: fix_conserved_domain_info")
        in_file = open(data_path + data_name)
        out_file = open(fixed_file, FileMode.WRITE.value)
        if delete_first_line:
            in_file.readline()
        for line in in_file:
            row = line.rstrip("\n").split("\t")
            c_elegans_gene = row[2]
            human_gene = row[0]
            row[domain_column] = str(
                float(c_elegans_domains_dic[c_elegans_gene]) / float(human_domains_dic[human_gene]))
            out_file.write("\t".join(row) + "\n")
        in_file.close()
        out_file.close()

    @staticmethod
    def fix_conserved_domains_file():
        print("def: fix_conserved_domain_info")
        human_id_ENSG = FileReader(FileReader.research_path + r"\Data",
                                   r"\human-genes-and-conserved-domains-230619",
                                      FileType.TSV).from_file_to_dict(0, 1)
        c_elegans_id_WB = FileReader(FileReader.research_path + r"\Data",
                                     r"\c-elegans-genes-and-conserved-domains-230619",
                                        FileType.TSV).from_file_to_dict(0, 1)
        f = open(FileReader.research_path + r"\Executors\data-230619-fixed-ratio-with-C-elegans-phenotypes")
        new_file = open("data-250619-fixed-domains", FileMode.WRITE.value)
        f.readline()
        for line in f:
            info = line.strip("\n").split("\t")
            human_id = info[0]
            c_elegans_id = info[2]
            info[4] = str(float(c_elegans_id_WB[c_elegans_id]) / float(human_id_ENSG[human_id]))
            new_file.write("\t".join(info) + "\n")
        f.close()
        new_file.close()

    @staticmethod
    def get_c_elegans_description_for_gene_id(gene_id, start_term: str = "\"text\":\"", end_term: str = "\",\"evidence\""):
        print("def: get_c_elegans_description_for_gene_id")
        if not gene_id:
            return "No description was found"
        try:
            record = HttpRequester("http://rest.wormbase.org/rest/widget/gene/" +
                                    gene_id + "/overview").make_request()
            print(record)
            if not record:
                return "No description was found"
            info_index = record.find("concise_description")
            shorter_info = record[info_index + len("concise_description"):]
            start_index = shorter_info.find(start_term)
            end_index = shorter_info.find(end_term)
            description = shorter_info[start_index + len(start_term): end_index]
            return description
        except Exception as e:
            print("Exception in get_c_elegans_description_for_gene_id:", e)
            return "no description was found"


    @staticmethod
    def get_c_elegans_gene_name_for_gene_id(gene_id):
        print("def: get_c_elegans_gene_name_for_gene_id")
        if not gene_id:
            return None
        try:
            record = HttpRequester("http://rest.wormbase.org/rest/widget/gene/" +
                                    gene_id + "/overview").make_request()
            #print("the start:")
            #print(record)
            #print("label" , record[record.find("label")])
            if not record:
                return None
            info_index = record.find("legacy_information")
            shorter_info = record[info_index + len("legacy_information"):]
            #print("\n",shorter_info)
            start_index = shorter_info.find("name")
            #print("\n",shorter_info[start_index:])
            next_index = shorter_info[start_index:].find("label")
            #print("\n",shorter_info[start_index:][next_index:])
            end_index = shorter_info[start_index:][next_index:].find('\","class"')
            description = (shorter_info[start_index:][next_index:])[len('"label":'): end_index]
            #print("\n", description)
            return description
        except Exception as e:
            print("Exception in get_c_elegans_description_for_gene_id:", e)
            return None


    # receives (1) a dictionary containing gene name as keys, and orthologs and phenotypes as values, (2) a
    # dictionary for human genes and their lengths, (3) dictionary for c.elegans genes and
    # their lengths and deletes genes for which the length differences between the human gene and its ortholog
    # are to high
    @staticmethod
    def filter_genes_with_size_differences(d: dict, human_genes: dict, c_elegans_genes: dict, p: int):
        print("def: filter_genes_with_size_differences")
        genesAndOrthologs = {}
        genesAndPhenotypes = {}
        d_keys = list(d.keys())
        for key in d_keys:
            try:
                humanGeneLength = float(human_genes[key])
            except Exception as e:
                print("Exception in filter_genes_with_size_differences:", e)
                print("gene " + key + " wasn't found...")
                continue
            d_values = d[key]
            for value in d_values:
                if isinstance(value, tuple):  # item is an ortholog
                    c_elegans_gene = value[0]
                    try:
                        cElegansGeneLength = float(c_elegans_genes[c_elegans_gene])
                    except Exception as e:
                        print("Exception in filter_genes_with_size_differences:", e)
                        print("gene " + c_elegans_gene + " wasn't found...")
                        continue
                    if (humanGeneLength * p) / 100 <= cElegansGeneLength:
                        if key in genesAndOrthologs:
                            genesAndOrthologs[key].append(value)
                        else:
                            genesAndOrthologs[key] = [value]
                    else:  # genes have too different sizes
                        print("C.elegans gene's length is " + str((cElegansGeneLength * 100) / humanGeneLength) +
                              " percent of the human gene's length")
                elif isinstance(value, str):
                    genesAndPhenotypes[key] = value
                else:
                    print("Problem has occured in function: filter_genes_with_size_differences")
                    exit()
        return genesAndOrthologs, genesAndPhenotypes

    # receives two dict: one converting human gene to c.elegans gene and one converting human gene to its variants
    # not sure
    @staticmethod
    def find_genes_variants_and_homologous(homologousDic: dict, variantasDic: dict):
        print("def: find_genes_variants_and_homologous")
        unitedDic = {}
        count = 0
        keys = variantasDic.keys()
        for key in keys:
            humanGene = key[0]
            variantName = key[1]
            if humanGene in homologousDic:
                unitedDic[humanGene, variantName] = [homologousDic[humanGene]] + variantasDic[key]
            else:
                if count < 100:
                    print("key: " + str(key) + str(variantasDic[key]))
                    count += 1
        return unitedDic

    # receivrd a dict of human genes and their conditions as values, and returns a set of all possible conditions
    @staticmethod
    def get_conditions_set(conditionsFilteredDict: dict, condition_column):
        print("def: get_conditions_set")
        s = set()
        values = conditionsFilteredDict.values()
        for value in values:
            s.add(value[condition_column])

        return s

    # receives two dicts: one converting HGNC gene to its conditions and one converting human gene to c.elegans gene
    # and returns a dict with HGNC name for key and orthologs and conditions as values
    @staticmethod
    def find_genes_with_valid_condition_and_homolog(HGNCGenesWConditions: dict, HomologousGenes: dict):
        print("def: find_genes_with_valid_condition_and_homolog")
        genes = {}
        for HCGNGene in HGNCGenesWConditions.keys():  # gene has valid condition
            if HCGNGene in HomologousGenes:  # gene has homoloug gene
                # new dict will have HCGN gene name for key and homoloug genes + condition as value
                genes[HCGNGene] = HomologousGenes[HCGNGene] + [HGNCGenesWConditions[HCGNGene][0]]
        return genes

