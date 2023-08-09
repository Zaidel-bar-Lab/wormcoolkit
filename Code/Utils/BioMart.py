# from pybiomart import Dataset

from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.Ensembl import Ensembl
from sys import platform

class BioMart:
    def __init__(self):
        # self.dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
        # self.df = self.dataset.query(attributes=['ensembl_gene_id', 'uniprotswissprot']).dropna()

        if platform == "darwin":
            self.human_id_to_uniprot_id = FileReader(FileReader.research_path + r"/Data",
                                                     r"/human_gene_id_gene_name_uniprot_id.txt").from_redundant_key_to_redundant_value(
                0, 2)
            self.worm_id_to_uniprot_id = FileReader(FileReader.research_path + r"/Data",
                                                    r"/c_elegans_gene_id_gene_name_uniprot_id.txt").from_redundant_key_to_redundant_value(
                0, 2)
        else:
            self.human_id_to_uniprot_id = FileReader(FileReader.research_path + r"\Data",
                                r"\human_gene_id_gene_name_uniprot_id.txt").from_redundant_key_to_redundant_value(0, 2)
            self.worm_id_to_uniprot_id = FileReader(FileReader.research_path + r"\Data",
                                r"\c_elegans_gene_id_gene_name_uniprot_id.txt").from_redundant_key_to_redundant_value(0, 2)

    # formerly - from the biomart package, now - from the biomart sheet
    def get_human_swissprot_sequence_from_biomart(self, gene_id):
        # uni_ids = self.df.loc[self.df['Gene stable ID'] == gene_id]['UniProtKB/Swiss-Prot ID']
        if gene_id not in self.human_id_to_uniprot_id:
            return None
        uni_ids = self.human_id_to_uniprot_id[gene_id]
        if len(uni_ids) == 1:
            return HttpRequester.get_protein_seq_by_uniprot_swissprot_id(uni_ids[0])
        elif len(uni_ids) == 0:
            return None
        else:  # many ids
            # check if uniprot reviewed all ids
            reviewed_html = HttpRequester.get_human_uniprot_html(gene_id)
            for uni_id in uni_ids[:]:
                if uni_id not in reviewed_html:
                    uni_ids.remove(uni_id)
            # choose by length
            chosen_seq = ''
            for uni_id in uni_ids:
                optional_seq = HttpRequester.get_protein_seq_by_uniprot_swissprot_id(uni_id)
                if len(optional_seq) > len(chosen_seq):
                    chosen_seq = optional_seq
            return chosen_seq

    def get_c_elegans_swissprot_sequence_from_biomart(self, gene_id):
        # uni_ids = self.df.loc[self.df['Gene stable ID'] == gene_id]['UniProtKB/Swiss-Prot ID']
        print("def: get_c_elegans_swissprot_sequence_from_biomart")
        print(gene_id)
        if gene_id not in self.worm_id_to_uniprot_id:
            print("1")
            return None
        uni_ids = self.worm_id_to_uniprot_id[gene_id]
        print(uni_ids)
        if len(uni_ids) == 1:
            print("2")
            return HttpRequester.get_protein_seq_by_uniprot_swissprot_id(uni_ids[0])
        elif len(uni_ids) == 0:
            print("3")
            return None
        else:  # many ids
            # check if uniprot reviewed all ids
            reviewed_html = HttpRequester.get_worm_uniprot_html(gene_id)
            for uni_id in uni_ids[:]:
                if uni_id not in reviewed_html:
                    uni_ids.remove(uni_id)
            # choose by length
            chosen_seq = ''
            for uni_id in uni_ids:
                optional_seq = HttpRequester.get_protein_seq_by_uniprot_swissprot_id(uni_id)
                if len(optional_seq) > len(chosen_seq):
                    chosen_seq = optional_seq
            return chosen_seq

    def get_human_protein_from_uniprot_by_gene_id(self, gene_id):
        # first try canonical sequence
        canonical_protein_sequence = self.get_human_swissprot_sequence_from_biomart(gene_id)
        if canonical_protein_sequence:
            return canonical_protein_sequence
        else:
            # get longest sequence
            longest_seq = HttpRequester().get_longest_human_protein_sequence_from_uniprot(gene_id)
            if not longest_seq:
                return None
            return longest_seq

    def get_human_protein_from_uniprot_by_gene_name(self, human_gene_name):
        human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
        if not human_gene_id:
            print("Human gene id for", human_gene_name, "cannot be found")
            return None
        else:
            return self.get_human_protein_from_uniprot_by_gene_id(human_gene_id)


