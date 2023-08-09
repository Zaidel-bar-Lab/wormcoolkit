import xlrd
from Code.CRISPR.NamedTuples.RestrictionEnzyme import RestrictionEnzyme
from Code.Enum.FileMode import FileMode
from Code.Enum.FileType import FileType
import unicodedata
from sys import platform
#import pandas as pd

chromosomes = [str(i) for i in range(1, 23)] + ['x', 'y']
FILE_TYPES_DELIMETER = {FileType.TSV: "\t", FileType.CSV: ",", FileType.UNCLEAR: "\" \""}


class FileReader:

    #research_path = r"C:\Users\Liran\PycharmProjects\Research"
    research_path = r"C:\Users\sapir\Desktop\new_project"
    # research_path = r"/Users/michaelrokitko/Desktop/Ksenia/SapirProject"

    def __init__(self, path, fileName, fileType: FileType = FileType.TSV):
        self.path = path
        self.name = fileName
        self.type = fileType
        self.size = -1  # not defined yet

    # extracts all genes with confidence in orthology to C.elegans's genes
    # returns a dictionary with key: human gene id, value: c.elegans gene id
    def read_genes_with_orthology_confidence(self):
        print("def: read_genes_with_orthology_confidence")
        orthologous_genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        print("sanity check: " + f.readline())  # sanity check
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])

            if line[10] == '1':  # confidence in orthology
                orthologous_genes[line[0]] = line[4]  # genes
        f.close()
        return orthologous_genes

    # receiving an open file for reading, extracts all genes with alleles and puts them in the dictionary
    # the dictionary contains key: gene id, value:chromosome, start, end, variantName, original bp, alleles, place
    @staticmethod
    def read_genes_with_variants(file, genes_and_variants, variantType):
        print("def: read_genes_with_variants")
        file.readline()  # skip the headlines
        for row in file:
            line = row.rstrip('\n').split('\t')

            if line[1] != "":  # variant exists
                pass
                variant_alleles = line[5]
                if variant_alleles == "COSMIC_MUTATION":
                    alleles = ["cosmic", "cosmic"]
                else:
                    alleles = variant_alleles.split("/")
                    # gene info = chromosome, start, end, original bp, alleles, place, somatic/germ line
                    geneInfo = [line[2], line[3], line[4], alleles[0], alleles[1:], line[7], variantType]
                    genes_and_variants[(line[0], line[1])] = geneInfo

    # receives the path and goes through all number of chromosomes
    def read_all_genes_with_variants(self):
        print("def: read_all_genes_with_variants")
        genes_and_variants = {}
        variant_types = ["somatic", "germline"]
        for variant_type in variant_types:
            for chromosome_number in chromosomes:
                if platform == "darwin":
                    filePath = self.path + "/" + variant_type + self.name + str(chromosome_number) + ".txt"
                else:
                    filePath = self.path + "\\" + variant_type + self.name + str(chromosome_number) + ".txt"
                try:
                    file = open(filePath, FileMode.READ.value)
                    self.read_genes_with_variants(file, genes_and_variants, variant_type)
                except IOError:
                    print("Oops, File " + filePath + " does not exist :O")
        return genes_and_variants

    # half conditions - conditions that are to be filtered only when they are alone
    def read_file_filter_conditions(self, condition_column, conditions):
        print("def: read_file_filter_conditions")
        genes = {}
        if self.type == FileType.XLS:
            wb = xlrd.open_workbook(self.path + self.name)
            sheet = wb.sheet_by_index(0)
        else:
            print("not supported yet")
            exit()
        for row in range(1, sheet.nrows):
            gene_conditions = sheet.cell_value(row, condition_column)
            invalid_conditions = 0
            listed_gene_conditions = FileReader.\
                list_cleaner(gene_conditions.replace("?", ",").replace("|", ",").split(","))
            for gene_condition in listed_gene_conditions:
                for condition in conditions:
                    if condition in gene_condition or condition.lower() in gene_condition:
                        invalid_conditions += 1
                        if "modifier" in condition or "susceptibility" in condition or "risk" in condition:
                            invalid_conditions += 1
            if invalid_conditions < len(listed_gene_conditions):
                genes[sheet.cell_value(row, 0)] = sheet.row_values(row, 1)  # column 0 is the gene name
        return genes

    # receives a list and returns the list without items that are empty or made of spaces
    @staticmethod
    def list_cleaner(l: list):
        print("def: list_cleaner")
        new_list = []
        for item in l:
            if item.strip():
                new_list.append(item)
        return new_list

    # read the ortholist excel page, and returns two dicts: one converting ensembleId of a human gene to its
    # HGNC symbol, and one converting from human gene HGNC symbol to its c-elegans ortholog and number of programs
    # supporting that claim

    def read_ortholist(self):
        print("def: read_ortholist")
        human_genes = {}
        homologous_genes = {}
        wb = xlrd.open_workbook(self.path + self.name)
        sheet = wb.sheet_by_index(0)
        for row in range(1, sheet.nrows):
            hgnc_symbol = sheet.cell_value(row, 5)
            ensembl_id = sheet.cell_value(row, 4)
            if ensembl_id not in human_genes:
                human_genes[ensembl_id] = hgnc_symbol  # Ensembl ID
                homologous_genes[hgnc_symbol] = []
            # 0 - WormBase ID, 6 - number of programs out of 6
            homologous_genes[hgnc_symbol].append(tuple([sheet.cell_value(row, 0), int(sheet.cell_value(row, 6))]))
        return human_genes, homologous_genes

    # receives a table of human genes and returns a dict that converts the human gene name to its id
    def genes_id_and_names(self):
        print("def: genes_id_and_names")
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        f.readline()  # headlines
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            genes[line[1]] = line[0]  # geneName = geneId
        f.close()
        return genes

    # returns a dict with gene names as keys and their length in bp as values
    def get_genes_length(self, gene_name_index, start_index, end_index):
        print("def: get_genes_length")
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        f.readline()  # headlines
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            # line[3] - gene end, line[2] = gene start
            genes[line[gene_name_index]] = int(line[end_index])-int(line[start_index])
        f.close()
        return genes

    # returns a list of all c.elegans genes
    def get_genes_list(self, column: int = 0, delete_first: bool = False):
        print("def: get_genes_list")
        genes = []
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first:
            f.readline()
        for row in f:
            genes.append(row.strip("\n").split(FILE_TYPES_DELIMETER[self.type])[column])
        f.close()
        return genes

    def get_c_elegans_genes_from_output(self, column):
        print("def: get_c_elegans_genes_from_output")
        genes = []
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            genes.extend(FileReader.fromStringToGenesList(line[column], 0))
        f.close()
        return genes

    # receives (1) a key index, (2) value index, (3) boolean value determining whether to disregard the first sentence,
    # and returns a dictionary of all genes and the longest coding sequence length
    def get_genes_cd_length(self, key_index, value_index, delete_first_line: bool = False):
        print("def: get_genes_cd_length")
        lengths = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            try:
                gene = line[key_index]
                str_length = line[value_index]
                if str_length:
                    length = int(str_length)
                if gene in lengths:
                    if length > lengths[gene]:
                        lengths[gene] = length
                else:
                    lengths[gene] = length
            except Exception as e:
                print("There has been a problem with extracting the following gene's length:", gene, e)
        f.close()
        return lengths

    def from_file_to_dict(self, key_index: int = 0, value_index: int = 1, delete_first_line: bool = False):
        print("def: from_file_to_dict")
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            try:
                genes[line[key_index]] = line[value_index]
            except Exception as e:
                print("An error has occured with reading the file", self.name, "in row", row, ":", e)
        f.close()
        return genes

    def from_file_to_dict_with_tuple_key(self, first_key_index: int = 0, second_key_index: int = 1, value_index: int = 2,
                                         delete_first_line: bool = False):
        print("def: from_file_to_dict_with_tuple_key")
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            try:
                genes[(line[first_key_index], line[second_key_index])] = line[value_index]
            except Exception as e:
                print("An error has occured while reading the file", self.name, "in row", row, ":", e)
        f.close()
        return genes

    def from_file_to_list(self, delete_first_line: bool = False):
        print("def: from_file_to_list")
        keys = []
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            keys.append(row.rstrip('\n'))
        return keys

    def from_file_to_dict_with_plural_values(self, key_index, value_index, delete_first: bool = False):
        print("def: from_file_to_dict_with_plural_values")
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first:
            f.readline()  # headline
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            try:
                if line[key_index] not in genes:
                    genes[line[key_index]] = {line[value_index]}
                else:
                    genes[line[key_index]].add(line[value_index])
            except Exception as e:
                print("An error has occured while reading the file", self.name, "in row", row, ":", e)
        f.close()
        return genes

    def from_MMP_file_to_dict_with_listed_values(self, key_index, list_of_value_indexes: list, delete_first: bool = False):
        print("def: from_MMP_file_to_dict_with_listed_values")
        dic = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first:
            f.readline()  # headline
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            if line[7] == "intergenic":
                # not in a gene
                continue
            values = []
            for value_index in list_of_value_indexes:
                values.append(line[value_index])
            dic[line[key_index].strip(" ")] = values
        f.close()
        return dic

    def from_file_to_dict_with_lists_as_values(self, key_index, value_index, delete_first: bool = False):
        print("def: from_file_to_dict_with_lists_as_values")
        dic = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first:
            f.readline()  # headline
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            key = line[key_index]
            values = line[value_index].strip("\n").split(", ")
            dic[key] = values
        f.close()
        return dic

    def read_results_file(self):
        print("def: read_results_file")
        try:
            f = open(self.path + self.name, FileMode.READ.value)
        except Exception as e:
            print("file cannot be opened, maybe it doesn't exist:", e)
            exit()
        return f

    def get_celegans_id_to_human_name_dict(self, value_index: int = 1, key_index: int = 2):
        print("def: get_celegans_id_to_human_name_dict")
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            cElegansGenesTuples = FileReader.fromStringToTuplesList(line[key_index])
            for t in cElegansGenesTuples:
                genes[t[0]] = line[value_index]
        f.close()
        return genes

    def sizeOfFile(self):
        print("def: sizeOfFile")
        count = 0
        f = open(self.name, FileMode.READ.value)
        for _ in f:
            count += 1
        f.close()
        self.size = count
        return count

    @staticmethod
    def fromStringToTuplesList(s: str):
        print("def: fromStringToTuplesList")
        tuples_list = []
        str_list = s.replace("(", "").replace("'", "").replace(" ", "").split("),")
        for item in str_list:
            item = item.replace(")", "")
            tuples_list.append(tuple(item.split(",")))
        return tuples_list

    @staticmethod
    def fromStringToGenesList(s: str, c_elegans_column: int):
        print("def: fromStringToGenesList")
        genes = []
        tuples_list = FileReader.fromStringToTuplesList(s)
        for t in tuples_list:
            genes.append(t[c_elegans_column])
        return genes

    def readRedundantFile(self, key_index, value_index, delete_first_line: bool = False):
        print("def: readRedundantFile")
        genes = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            if line[key_index] not in genes:
                genes[line[key_index]] = line[value_index]
        f.close()
        return genes

    def fromPluralFileToOneFile(self, num_of_files, target_file_path):
        print("def: fromPluralFileToOneFile")
        target_file = open(target_file_path, FileMode.WRITE.value)
        for i in range(num_of_files):
            f = open(self.path+self.name+str(i))
            for row in f:
                target_file.write(row)
            print("File number " + str(i) + " is now completely downloaded")
            f.close()
        target_file.close()

    def extractCelegansFromData(self, orthologs_column, gene_column):
        print("def: extractCelegansFromData")
        genes = []
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            c_elegans_part = line[orthologs_column]
            genes.extend(FileReader.fromStringToGenesList(c_elegans_part, gene_column))
        f.close()
        return genes

    def makeDictFromSummary(self, key_index, value_index, delete_first_line: bool = False):
        print("def: makeDictFromSummary")
        genesDescription = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            try:
                genesDescription[line[key_index].strip("\"")] = line[value_index]
            except Exception as e:
                print("Exception in makeDictFromSummary:", e)
                f.close()
                exit()
        f.close()
        return genesDescription

    def addSummaryToGenesAndSeparateTuples(self, summary_path, new_data_plus_summary, key_index, summary_index,
                                           summary_delimeter, delete_first_line: bool = False):
        print("def: addSummaryToGenesAndSeparateTuples")
        summary = open(summary_path)
        summary_dic = FileReader(summary_path, "", summary_delimeter).makeDictFromSummary(key_index, summary_index,
                                                                                          delete_first_line)
        summary.close()

        data = open(self.path + self.name, FileMode.READ.value)
        new_data = open(new_data_plus_summary, FileMode.WRITE.value)

        for row in data:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            human_gene_id = line[0]
            human_gene_name = line[1]
            orthologs = FileReader.fromStringToTuplesList(line[2])
            phenotype = line[3]

            for t in orthologs:
                ortholog_id = t[0]
                if ortholog_id in summary_dic and len(summary_dic[ortholog_id]) > 1:
                    summary_value = summary_dic[ortholog_id]
                else:
                    summary_value = "Not mentioned"
                new_data.write(human_gene_id + "\t" + human_gene_name + "\t" + ", ".join(t) + "\t" +
                                phenotype + "\t" + summary_value + "\n")
        data.close()
        new_data.close()

    @staticmethod
    def lineFixer(line):
        print("def: lineFixer")
        new_line = ''
        line = unicodedata.normalize("NFKD", line)
        for ch in line:
            if ('A' <= ch <= 'z') or ch == " " or '1' <= ch <= '9':
                new_line += ch
        return new_line

    def fromHGMDtoDict(self):
        print("def: fromHGMDtoDict")
        genes_and_variants = {}
        f = open(self.path + self.name, FileMode.READ.value)
        line = f.readline().rstrip("\n")
        while line:
            if line.endswith(":"):
                gene = line.strip("\n")[:line.find(":")]
                genes_and_variants[gene] = []
                line = f.readline()
                while not line.startswith("\n") and len(line) > 3:
                    variant = line.strip("\n").split("\t")[1]
                    genes_and_variants[gene].append(variant)
                    line = f.readline()
                line = f.readline().rstrip("\n")
        f.close()
        return genes_and_variants

    def readData(self, key_index, delete_first_line: bool = False):
        print("def: readData")
        genes_data = {}
        f = open(self.path + self.name, FileMode.READ.value)
        if delete_first_line:
            f.readline()
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            key = line[key_index]
            genes_data[key] = row.strip("\n")
        f.close()
        return genes_data

    # # raises errors in A2HOSTING
    # def get_list_from_excel_using_pandas(self, column_name='WormBase Gene ID', sheet_name='kinase'):
    #     xls = pd.ExcelFile(self.path + self.name)
    #     df = pd.read_excel(xls, sheet_name)
    #     #df = df.sort_values(["Public Name"])
    #     return [gene_id for gene_id in df[column_name]]

    # raises errors in A2HOSTING
    # def get_dictionary_from_excel_using_pandas(self, key_column_name, value_column_name, sheet_name='kinase'):
    #     print(key_column_name, value_column_name)
    #     xls = pd.ExcelFile(self.path + self.name)
    #     df = pd.read_excel(xls, sheet_name)
    #     two_list = [(key, value) for key, value in zip(df[key_column_name], df[value_column_name])]
    #     return dict(two_list)

    # reads parsed enzymes and returns list of objects
    def get_parsed_restriction_enzymes_list(self, name_column=0, site_column=1, parsed_sites_column=2,
                                            full_site_column=3):
        print("def: get_parsed_restriction_enzymes_list")
        restriction_enzymes = []
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            name = line[name_column]
            filtered_site = line[site_column]
            full_site = line[full_site_column]
            parsed_sites = tuple(line[parsed_sites_column].split(","))
            restriction_enzymes.append(RestrictionEnzyme(name,
                                                         filtered_site,
                                                         parsed_sites,
                                                         full_site))
        f.close()
        return restriction_enzymes

    # calculates parsed enzymes and returns list of objects
    def get_restriction_enzymes_list(self, name_column=0, site_column=1):
        print("def: get_restriction_enzymes_list")
        restriction_enzymes = []
        f = open(self.path + self.name, FileMode.READ.value)
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            name = line[name_column]
            site = line[site_column]
            filtered_site = ''
            for ch in site:
                if 'A' <= ch <= 'Z':
                    filtered_site += ch

            restriction_enzymes.append(RestrictionEnzyme(name,
                                                         filtered_site,
                                                         self.find_enzymes_derivatives(filtered_site),
                                                         site))
        f.close()
        return restriction_enzymes

    @staticmethod
    def is_derivative(word, sample, i=0):
        print("def: is_derivative")
        variables = {'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
                     'H': ('A', 'C', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                     'N': ('A', 'C', 'G', 'T'), 'R': ('A', 'G'), 'S': ('C', 'G'),
                     'V': ('A', 'C', 'G'), 'W': ('A', 'T'), 'Y': ('C', 'T')}
        if i < len(word):
            if word[i] in variables:
                terminals = variables[word[i]]
                for terminal in terminals:
                    new_word = list(word)
                    new_word[i] = terminal
                    if FileReader.is_derivative("".join(new_word), sample, i+1):
                        return True
                return False
            else:
                if FileReader.is_derivative(word, sample, i+1):
                    return True
        return word.replace("/", "") == sample

    @staticmethod
    def find_enzymes_derivatives(word, rec=True):
        print("def: find_enzymes_derivatives")
        variables = {'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
                     'H': ('A', 'C', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                     'N': ('A', 'C', 'G', 'T'), 'R': ('A', 'G'), 'S': ('C', 'G'),
                     'V': ('A', 'C', 'G'), 'W': ('A', 'T'), 'Y': ('C', 'T')}
        if rec:
            s = set()
            return tuple(FileReader.find_derivatives_rec(word, 0, s, variables))
        else:
            return tuple(FileReader.find_derivatives_iter(word, variables))

    # the recursive function that parses word with non terminals to terminals
    @staticmethod
    def find_derivatives_rec(word, i, s, variables):
        print("def: find_derivatives_rec")
        if i == len(word) - 1:
            if word[i] in variables:
                new_word = list(word)
                terminals = variables[word[i]]
                for terminal in terminals:
                    new_word[i] = terminal
                    s.add("".join(new_word))
            else:
                s.add(word)
        else:
            if word[i] in variables:
                new_word = list(word)
                terminals = variables[word[i]]
                for terminal in terminals:
                    new_word[i] = terminal
                    FileReader.find_derivatives_rec("".join(new_word), i + 1, s, variables)
            else:
                FileReader.find_derivatives_rec(word, i + 1, s, variables)
        return s

    # the iterative function that parses word with non terminals to terminals
    @staticmethod
    def find_derivatives_iter(word, variables):
        print("def: find_derivatives_iter")
        s = {word}
        for i in range(len(word)):
            if word[i] in variables:
                terminals = variables[word[i]]
                for derivative in s.copy():
                    for terminal in terminals:
                        new_derivative = list(derivative)
                        new_derivative[i] = terminal
                        s.add("".join(new_derivative))
                    s.remove(derivative)
        return s

    @staticmethod
    def get_plain_site(site):
        print("def: get_plain_site")
        plain_site = ''
        for ch in site:
            if 'A' <= ch <= 'Z':
                plain_site += ch
        return plain_site

    @staticmethod
    def get_human_genes_ids():
        print("def: get_human_genes_ids")
        return FileReader(FileReader.research_path + r"\Data",
                          r"\human_gene_id-gene_name.txt").get_genes_list(0, True)

    @staticmethod
    def get_c_elegans_genes_ids():
        print("def: get_c_elegans_genes_ids")
        return FileReader(FileReader.research_path + r"\Data",
                          r"\c-elegans-gene-id_gene-name_ncbi-id.txt").get_genes_list(0, True)

    def from_redundant_key_to_redundant_value(self, key, value, delete_first = True):
        print("def: from_redundant_key_to_redundant_value")
        dic = {}
        f = open(self.path + self.name)
        if delete_first:
            f.readline()  # headline
        for row in f:
            line = row.rstrip('\n').split(FILE_TYPES_DELIMETER[self.type])
            if not line[value]:
                continue
            if line[key] in dic:
                if line[value] not in dic[line[key]]:
                    dic[line[key]].append(line[value])
            else:
                dic[line[key]] = [line[value]]
        return dic




