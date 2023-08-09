from Bio import Entrez
from Bio import pairwise2
from Bio.Blast import NCBIWWW
import os;
os.environ["OMP_NUM_THREADS"] = "1"
from Bio.Blast import NCBIXML
from Bio.Blast import Record
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SubsMat import MatrixInfo as matlist
#from Bio.Align import substitution_matrices as matlist
from Bio.pairwise2 import format_alignment
from sys import platform
from Code.Enum.BlastType import BlastType
from Code.Enum.FileType import FileType
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.Ensembl import Ensembl
from Code.Utils.Strings import Strings
from Code.Utils.BioMart import BioMart

E_VALUE_THRESH = 0.04


class BioPython:
    def __init__(self):
        if platform == "darwin":
            self.c_elegans_id_multiple_accessions = FileReader(FileReader.research_path + r"/Data",
                                                               r"/all-genes-ids-and-accession-numbers.txt",
                                                               FileType.TSV).from_file_to_dict_with_plural_values(0, 1)
        else:
            self.c_elegans_id_multiple_accessions = FileReader(FileReader.research_path + r"\Data",
                                                           r"\all-genes-ids-and-accession-numbers.txt",
                                                           FileType.TSV).from_file_to_dict_with_plural_values(0, 1)
        self.blast_results = {}
        self.bm = BioMart()

    @staticmethod
    def get_human_gene_sequence(geneName: str, converter: dict):
        print("def: get_human_gene_sequence")
        #print("input: geneName: str, converter: dict")
        if geneName in converter:
            geneId = converter[geneName]
        else:
            print("There is no record for this gene in the dictionary")
            return None
        Entrez.email = "liranavda@gmail.com"
        handle = Entrez.esearch(db="gene", retmax=10, term=geneId)
        record = Entrez.read(handle)
        handle.close()
        recordId = record["IdList"][0]
        gene = Entrez.efetch(db="gene", id=recordId, rettype="gb", retmode="text")
        data = gene.read()
        return data

    # receives a list of genes' ids(WB) and returns a dict of a gene id(WB) and its id(number)
    @staticmethod
    def get_genes_id(c_elegans_genes: list):
        print("def: get_genes_id")
        print("we have " + str(len(c_elegans_genes)) + " genes to extract id to. let's begin!\n")
        genes = {}
        index = 0
        for gene in c_elegans_genes:
            Entrez.email = "liranavda@gmail.com"
            gene_id = BioPython.get_gene_id_by_entrez(gene)
            if gene_id:
                genes[gene] = gene_id
                index += 1
                if not index % 100:
                    print("got to " + str(index) + "!\n")
            else:
                print("Gene id can't be found using Entrez: " + gene)
                continue
        return genes

    # retrieves the id (number)
    @staticmethod
    def get_gene_id_by_entrez(c_elegans_gene_id: str):
        print("def: get_gene_id_by_entrez")
        print("the input is: ", c_elegans_gene_id)
        Entrez.email = "liranavda@gmail.com"
        try:
            handle = Entrez.esearch(db="gene", retmax=1, term=c_elegans_gene_id)
            record = Entrez.read(handle)
            handle.close()
            if len(record["IdList"]) > 0:
                gene_ncbi_id = record["IdList"][0]
                return gene_ncbi_id
        except Exception as e:
            print("Exception in get_gene_id_by_entrez:", e)
            return None

    def get_human_gene_id_from_ENTREZ_by_name(gene_name):
        print("def: get_human_gene_id_from_ENTREZ_by_name")
        print(gene_name)
        Entrez.email = "liranavda@gmail.com"
        length = len(gene_name)
        search_string = gene_name + "[gene] AND Homo sapiens"
        # print(search_string)
        try:
            handle = Entrez.esearch(db="gene", term=search_string)
            record = Entrez.read(handle)
            handle.close()
            print(record)
            if len(record["IdList"]) > 0:
                ids = record['IdList']
                for i in range(len(ids)):
                    seq_id = ids[i]
                    handle = Entrez.efetch(db="gene", id=seq_id, retmode="text")
                    record = handle.read()
                    # print(record[4:4+length])
                    if (record[4:4 + length] == gene_name):     ## gene name start from chr 4
                        print(seq_id)
                        return (seq_id)
        except Exception as e:
            print("Exception in get_gene_id_by_entrez:", e)
            return None

    @staticmethod
    def blast(program, db, accession_number):
        print("def: blast")
        print("Accession Number: " + str(accession_number))
        try:
            result_handle = NCBIWWW.qblast(program=program, database=db, sequence=accession_number)
        except ValueError:
            result_handle = NCBIWWW.qblast(BlastType.fromNtoP(program), db, accession_number)
        record = NCBIXML.read(result_handle)
        result_handle.close()
        for alignment in record.alignments:
            if alignment.hit_def.startswith("Homo sapiens"):
                return True
        return False
        # print("****Alignment****")
        # print("sequence: ", alignment.title)
        # print("length: ", alignment.length)
        # for hsp in alignment.hsps:
        # if hsp.expect < E_VALUE_THRESH:
        #     print("****HSP****")
        #     print("e value:", hsp.expect)
        #     print(hsp.query[0:75] + "...")
        #     print(hsp.match[0:75] + "...")
        #     print(hsp.sbjct[0:75] + "...")

    @staticmethod
    def blast_with_seq(program, db, accession_number, sequence):
        print("def: blast_with_seq")
        state = False
        print("Accession Number: " + str(accession_number))
        result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequence,
                                       entrez_query="Homo sapiens[organism]")
        record: Record = NCBIXML.read(result_handle)
        result_handle.close()
        for alignment in record.alignments:
            if BioPython.get_species(alignment.hit_def).startswith("Homo sapiens"):
                state = True
                print("****Alignment****")
                print(alignment.hit_id)
                print(alignment.hit_def)
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print("****HSP****")
                        print("e value:", hsp.expect)
                        print("hsp score:", hsp.score)
        return state

    def blast_with_sequences(self, genes_names_sequences_dict, program="blastp", db="nr", query_organism: str = "Homo sapiens[organism]",
                                org_name: str = "Homo sapiens"):
        print("def: blast_with_sequences")
        sequences_string = ''
        for gene_name in genes_names_sequences_dict:
            sequences_string += ">" + gene_name + "\n" + genes_names_sequences_dict[gene_name] + "\n"
        print("sequences:", sequences_string)
        records_list = []
        try:
            print(program)
            print(db)
            print(sequences_string)
            print(query_organism)
            #save_file = open("C:/Users/sapir/Desktop/new_project/Data/blast_test.xml", 'w')
            #result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequences_string,
            #                               entrez_query=query_organism , hitlist_size = 3 , alignments=3 )
            result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequences_string,
                                           entrez_query=query_organism)

            print("check blast resultsssss: ")
            print(result_handle)
            #blast_results = result_handle.read()
            #save_file.write(blast_results)
            #save_file.close()

            print("blaaaaaaaa")
            records_list = list(NCBIXML.parse(result_handle))
            print(records_list)
            result_handle.close()
            print("done qblast")
        except Exception as e:
            print(e)

        for record in records_list:
            hit_ids = []
            alignments = record.alignments
            for alignment in alignments:
                try:
                    if BioPython.get_species(alignment.hit_def).startswith(org_name):
                        hit_id_record = alignment.hit_id
                        hit_id = hit_id_record[
                                 hit_id_record.find("|") + 1:hit_id_record.find("|", hit_id_record.find("|") + 1)]
                        extra_letter = hit_id_record[hit_id_record.find(hit_id) +
                                                     len(hit_id) + 1:hit_id_record.find(hit_id) + len(hit_id) + 2]
                        if extra_letter:
                            hit_id = hit_id + "_" + extra_letter
                        hit_ids.append(hit_id)
                except Exception as e:
                    print("Exception in blast_with_sequences:", e)
                    continue
            self.blast_results[record.query] = hit_ids
        print("blast results:", self.blast_results)

    def build_genes_hit_ids_dictionary(self, genes_names_to_blast, key_species="C.elegans"):
        print("def: build_genes_hit_ids_dictionary")
        if key_species == "C.elegans":
            genes_sequences_dict = self.get_aa_seqs_dict_by_c_elegans_genes_names(genes_names_to_blast)
        else:  # human genes
            genes_sequences_dict = self.get_aa_seqs_dict_by_human_genes_names(genes_names_to_blast)

        query_organism = "Homo sapiens[organism]" if key_species == "C.elegans" else "Caenorhabditis elegans[organism]"
        org_name = "Homo sapiens" if key_species == "C.elegans" else "Caenorhabditis elegans"
        self.blast_with_sequences(genes_sequences_dict, query_organism=query_organism, org_name=org_name)


    # receives (1) an accession numbers and (2) the relevant sequencec, runs BLAST and returns a list of all
    # id of homo sapiens reported (first 50 by default)
    @staticmethod
    def pipeline_blast_with_seq(program, db, sequence, query_organism: str = "Homo sapiens[organism]",
                                org_name: str = "Homo sapiens"):
        print("def: pipeline_blast_with_seq")
        hit_ids = []
        try:
            result_handle = NCBIWWW.qblast(program=program, database=db, sequence=sequence,
                                           entrez_query=query_organism)
            record: Record = NCBIXML.read(result_handle)
            result_handle.close()
        except Exception as e:
            print("Exception in pipeline_blast_with_seq:", e)
            return None
        for alignment in record.alignments:
            try:
                if BioPython.get_species(alignment.hit_def).startswith(org_name):
                    hit_id_record = alignment.hit_id
                    hit_id = hit_id_record[hit_id_record.find("|") + 1:hit_id_record.find("|", hit_id_record.find("|") + 1)]
                    extra_letter = hit_id_record[hit_id_record.find(hit_id) +
                                                len(hit_id) + 1:hit_id_record.find(hit_id) + len(hit_id) + 2]
                    if extra_letter:
                        hit_id = hit_id + "_" + extra_letter
                    hit_ids.append(hit_id)
            except Exception as e:
                print("Exception in pipeline_blast_with_seq:", e)

                continue
        return hit_ids

    @staticmethod
    def blast_for_all(blast_type, db, accessions_dict: dict):
        print("def blast_for_all")
        f = open("blast-for-all", "w+")
        index = 0
        keys = accessions_dict.keys()
        for key in keys:
            accession = accessions_dict[key]
            res = BioPython.blast(blast_type, db, accession)
            index += 1
            if not index % 100:
                print("got to " + str(index) + "!")
            f.write(key + "\t" + accession + "\t" + str(res))
            print(res)

    # receives a dictionary of accession numbers and sequences
    @staticmethod
    def blastp_by_accessions(blast_program, db, accessions_dict: dict):
        print("def: blastp_by_accessions")
        accessions = accessions_dict.keys()
        for accession in accessions:
            seq = accessions_dict[accession]
            print("seq: " + seq)
            res = BioPython.blast_with_seq(blast_program, db, accession, seq)
            print(res)

    # receives (1) a sequence with white spaces, counts only the nucleotide letters, and returns the count
    @staticmethod
    def get_sequence_length(seq: str):
        print("def: get_sequence_length")
        return seq.count('t') + seq.count('a') + seq.count('c') + seq.count('g')

    # receives (1) accession numbers for specific gene id, enact the function that chooses the longest isoform
    # of all accession number, and returns the chosen longest accession number
    @staticmethod
    def get_longest_accession(accessions: set(), index_word ="ORIGIN", database: str = "nucleotide"):
        print("def: get_longest_accession")
        Entrez.email = "liranavda@gmail.com"
        longest_isoform = None
        longest_isoform_length = 0
        if len(accessions) == 1:
            return accessions.pop()
        for accession_number in accessions:
            info = Entrez.efetch(db=database, id=accession_number, rettype="gb", retmode="text")
            record = info.read()
            info.close()
            if index_word in record:
                index = record.find(index_word)
                nt_seq = record[index + len(index_word):]
                seq_length = BioPython.get_sequence_length(nt_seq)
                if seq_length > longest_isoform_length:
                    longest_isoform = accession_number
                    longest_isoform_length = seq_length
            else:
                print(index_word, "does not exist in record for", accession_number)
        if not longest_isoform:  # longest isoform stayed the empty string
            print("Error: no isoform was found longer than the empty string")
        return longest_isoform

    @staticmethod
    def get_aa_seq_from_ensembl(gene_id , isoform_name=None):
        print("def: get_aa_seq_from_ensembl")
        print(isoform_name)
        if isoform_name == None:
            chosen_seq = ''
            try:
                res = HttpRequester(url="https://rest.ensembl.org/sequence/id/").get_protein_sequence_from_ensembl(gene_id=gene_id)
                print("res1")
                print(res)
                if not res:
                    return None
            except Exception as e:
                print("Exception in get_aa_seq_from_ensembl:", e)
                return None
            sequences = res.split(">")
            for seq in sequences:
                if seq == "":
                    continue
                else:
                    print(seq)
                    parsed_seq = Strings.from_fasta_seq_to_seq(seq)
                    #print(parsed_seq)
                    if len(chosen_seq) < len(parsed_seq):
                        chosen_seq = parsed_seq
            return chosen_seq if chosen_seq else None
        else:
            print(isoform_name)
            chosen_seq = ''
            try:
                res = HttpRequester(url="https://rest.ensembl.org/sequence/id/").get_protein_sequence_from_ensembl(gene_id=gene_id)
                print("res2")
                #print(res)
                if not res:
                    return None
            except Exception as e:
                print("Exception in get_aa_seq_from_ensembl:", e)
                return None
            sequences = res.split(">")
            for seq in sequences:
                if seq == "":
                    continue
                else:
                    #print(seq)
                    split = seq.split("\n")
                    #print ("dd" , split[0] , isoform_name)
                    if split[0] == str(isoform_name+"."):
                        #print("seqqqqqqqqqqqqqqqq: " ,seq)
                        chosen_seq = Strings.from_fasta_seq_to_seq(seq)
                    #print("chosen_seq : ", chosen_seq)
                    #if len(chosen_seq) < len(parsed_seq):
                        #chosen_seq = parsed_seq
            return chosen_seq if chosen_seq else None

    @staticmethod
    def get_species(alignment_hit: str):
        print("def: get_species")
        start_index = alignment_hit.find("[")
        end_index = alignment_hit.find("]")
        species = alignment_hit[start_index + 1:end_index]
        if species == "imported":
                return BioPython.get_species(alignment_hit[end_index + 1:])
        return species

    @staticmethod
    def get_gene_name_from_protein_accession(hit_ids: [], search_term: str = "gene=\"", end_search_term: str = "\"\n"):
        print("def: get_gene_name_from_protein_accession")
        accessions_and_names = {}
        Entrez.email = "liranavda@gmail.com"
        for hit_id in hit_ids:
            try:
                info = Entrez.efetch(db="protein", id=hit_id, rettype="gb", retmode="text")
                record = info.read()
                info.close()
            except Exception as e:
                print("Exception in get_gene_name_from_protein_accession:", e)
                print("extraction was not successful for: " + hit_id)
                continue
            # extracting aa sequence out of chosen isoform
            start_index = record.find(search_term)
            if start_index < 0:  # no gene name in record
                continue
            end_index = record.find(end_search_term, start_index)
            gene_name = record[start_index + len(search_term):end_index]
            accessions_and_names[hit_id] = gene_name
        return accessions_and_names

    @staticmethod
    def get_alignments(human_seq, c_elegans_seq):
        print("def: get_alignments")
        classic_alignments = pairwise2.align.globalxx(human_seq, c_elegans_seq)
        pam_alignments = pairwise2.align.globaldx(human_seq, c_elegans_seq, matlist.pam250)
        print(matlist.pam250)
        with_penalty_alignments = pairwise2.align.globalms(human_seq, c_elegans_seq, 2, -1, -.5, -.1)
        blosum_alignments = pairwise2.align.globaldx(human_seq, c_elegans_seq, matlist.blosum62)
        print(matlist.blosum62)
        alignments = blosum_alignments + pam_alignments + classic_alignments + with_penalty_alignments
        return alignments

    # receives (1) the human sequence and (2) the C.elegans sequence, align them and check the maximum conservation out
    # of all the alignment, returns the score
    @staticmethod
    def get_conservation_score(human_seq, c_elegans_seq):
        print("def: get_conservation_score")
        alignments = BioPython.get_alignments(human_seq, c_elegans_seq)
        max_score = 0
        for alignment in alignments:
            score = BioPython.alignment_window_conservation_checker(alignment, 0, 0, whole_alignment=True)
            if score > max_score:
                max_score = score
        return max_score

    @staticmethod
    def pairwise_alignment_inspector(human_seq, c_elegans_seq, original_amino_acid, variant_index, window_size: int = 30):
        print("def: pairwise_alignment_inspector")
        answers_count = [0, 0, 0]
        ultimate_result = ''
        similar_answer_result = ''
        different_answer_result = ''
        wrong_answer_result = ''
        count = 0
        c_elegans_aa_location = '-'
        sum_window_conservation_scores = 0
        alignments = BioPython.get_alignments(human_seq, c_elegans_seq)
        for alignment in alignments:
            human_alignment = alignment[0]
            c_elegans_alignment = alignment[1]
            alignment_index = Strings.get_amino_acid_in_location_in_alignment(variant_index, human_alignment)
            status = human_alignment[alignment_index - 1], c_elegans_alignment[alignment_index - 1]

            if original_amino_acid == '-':  # termination
                return "Stop Codon", count, c_elegans_aa_location, 0
            if human_alignment[alignment_index - 1] != original_amino_acid:
                answers_count[2] += 1
                wrong_answer_result = "not in the human sequence, wanted: " + original_amino_acid + \
                         " and found: " + human_alignment[alignment_index - 1] + " alignment index " + str(alignment_index)
            elif human_alignment[alignment_index - 1] == c_elegans_alignment[
                        alignment_index - 1] == original_amino_acid:
                start = max(0, alignment_index - window_size)
                end = min(len(alignment[0]) - 1, alignment_index + window_size)
                tmp_window_conservation_score = BioPython.alignment_window_conservation_checker(alignment, start, end)
                sum_window_conservation_scores += tmp_window_conservation_score
                count += 1
                ultimate_result = "conserved: " + original_amino_acid
                c_elegans_aa_location = BioPython.get_place_in_sequence(c_elegans_alignment, alignment_index)

            elif Strings.are_amino_acids_similar(human_alignment[alignment_index - 1],
                                                 c_elegans_alignment[alignment_index - 1]):
                answers_count[0] += 1
                similar_answer_result = "not conserved, but similar: " + ", ".join(status)
            else:
                answers_count[1] += 1
                different_answer_result = "not conserved: " + ", ".join(status)

        if answers_count[0] >= answers_count[1] and answers_count[0] >= answers_count[2]:
            answer_result = similar_answer_result
        elif answers_count[1] >= answers_count[0] and answers_count[1] >= answers_count[2]:
            answer_result = different_answer_result
        else:
            answer_result = wrong_answer_result
        print("counts:", answers_count)
        return (ultimate_result, count, c_elegans_aa_location, sum_window_conservation_scores/count) if ultimate_result else \
            (answer_result, count, c_elegans_aa_location, 0)

    @staticmethod
    def get_place_in_sequence(c_elegans_alignment, alignment_index):
        print("def: get_place_in_sequence")
        sequence_index = 0
        for char_index in range(alignment_index):
            if 'A' <= c_elegans_alignment[char_index] <= 'z':
                sequence_index += 1
        return sequence_index

    @staticmethod
    def alignment_window_conservation_checker(alignment, start, end, whole_alignment=False):
        print("def: alignment_window_conservation_checker")
        lst = format_alignment(*alignment).split("\n")
        equals = 0
        if whole_alignment:
            # print(format_alignment(*alignment))
            segment_one = lst[0]
            segment_two = lst[2]
        else:
            segment_one = lst[0][start:end]
            # print(segment_one)
            segment_two = lst[2][start:end]
        for i in range(len(segment_one)):
            # print(segment_one[i], segment_two[i])
            if segment_one[i] == segment_two[i]:
                equals += 1
        return (equals / len(segment_one))*100

    # receives C.elegans ncbi gene id and return the longest isoform (accession number) of the gene
    def get_c_elegans_accession_number(self, c_elegans_gene_id):
        print("def: get_c_elegans_accession_number")
        if c_elegans_gene_id in self.c_elegans_id_multiple_accessions:
            c_elegans_accession_numbers = set(self.c_elegans_id_multiple_accessions[c_elegans_gene_id])
            c_elegans_accession_number = BioPython.get_longest_accession(accessions=c_elegans_accession_numbers)
            return c_elegans_accession_number
        else:
            return None

    # receives a list of C.elegans' genes and returns a dictionary with genes names as keys and sequences as values
    def get_aa_seqs_dict_by_c_elegans_genes_names(self, c_elegans_genes_names):
        print("def: get_aa_seqs_dict_by_c_elegans_genes_names")
        genes_names_and_sequences = {}
        for gene_name in c_elegans_genes_names:
            seq = self.get_aa_seq_by_c_elegans_gene_name(gene_name)
            if seq:
                genes_names_and_sequences[gene_name] = seq
        print(len(c_elegans_genes_names), "genes names, and gene name-sequence dictionary of",
              len(genes_names_and_sequences), "items")
        return genes_names_and_sequences

    def get_aa_seq_by_c_elegans_gene_name(self, c_elegans_gene_name , isoform_name=None):
        print("def: get_aa_seq_by_c_elegans_gene_name")
        gene_id_WB = Ensembl.get_c_elegans_gene_id_by_gene_name(c_elegans_gene_name)
        if not gene_id_WB:
            print("Couldn't find gene id (WB) for:", c_elegans_gene_name)
            return None
        return self.get_c_elegans_aa_seq_by_id(gene_id_WB , isoform_name)

    def get_aa_seqs_dict_by_human_genes_names(self, human_genes_names):
        print("def: get_aa_seqs_dict_by_human_genes_names")
        genes_names_and_sequences = {}
        for gene_name in human_genes_names:
            gene_id = Ensembl.get_human_gene_id_by_gene_name(gene_name)
            if gene_id:
                human_seq = self.get_human_aa_seq_by_id(gene_id)
                if human_seq:
                    genes_names_and_sequences[gene_name] = human_seq
        return genes_names_and_sequences

    # this function unites all the ways to extract a human sequence from human gene id
    def get_human_aa_seq_by_id(self, human_gene_id):
        print("def: get_human_aa_seq_by_id")
        seq = self.bm.get_human_protein_from_uniprot_by_gene_id(human_gene_id)  # try from uniprot
        if not seq:
            seq = BioPython.get_aa_seq_from_ensembl(human_gene_id)  # try from ensembl
            if not seq:
                print("Couldn't find gene sequence for gene: " + human_gene_id)
                return None
        return seq

    def get_human_aa_seq_by_name(self, human_gene_name):
        print("def: get_human_aa_seq_by_name")
        human_gene_id = Ensembl.get_human_gene_id_by_gene_name(human_gene_name)
        if not human_gene_id:
            print("Human gene id for", human_gene_name, "cannot be found")
            return None
        else:
            return self.get_human_aa_seq_by_id(human_gene_id)

    # receives human gene id (ENSG) and returns a set of all related accession numbers
    @staticmethod
    def get_human_accession_number(human_gene_id):
        print("def: get_human_accession_number")
        accessions_dic = FileReader(FileReader.research_path + r"\Data",
                                    r"\human_genes_ids_to_accession_numbers.txt").from_file_to_dict_with_plural_values(0, 1, True)
        if human_gene_id not in accessions_dic:
            print("Human gene id", human_gene_id, "didn't have an accession number")
            return None
        accessions = set(accessions_dic[human_gene_id])
        accession = BioPython.get_longest_accession(accessions=accessions)
        return accession

    # receives gene id (WB/ENSG) and returns gene ncbi id (number)
    @staticmethod
    def get_ncbi_id(gene_id):
        print("def: get_ncbi_id")
        c_elegans_id_number = Ensembl.get_ncbi_id_by_gene_id(gene_id)
        if c_elegans_id_number:
            return c_elegans_id_number
        else:
            c_elegans_id_number = BioPython.get_gene_id_by_entrez(gene_id)
            if c_elegans_id_number:
                return c_elegans_id_number
        print("Couldn't find C.elegans's ncbi id for gene id:", gene_id)
        return None

    # receives the C.elegans gene id (WB) and returns the gene's sequence
    def get_c_elegans_aa_seq_by_id(self, c_elegans_id_WB , isoform_name=None):
        print("def: get_c_elegans_aa_seq_by_id")
        print("c_elegans_seq = self.get_aa_seq_from_ensembl(c_elegans_id_WB, isoform_name)")
        c_elegans_seq = self.get_aa_seq_from_ensembl(c_elegans_id_WB, isoform_name)
        #print("c_elegans_seq: " , c_elegans_seq)
        if c_elegans_seq != None:
            return c_elegans_seq
        else:
            print("canonical_protein_sequence = self.bm.get_c_elegans_swissprot_sequence_from_biomart(c_elegans_id_WB)")
            canonical_protein_sequence = self.bm.get_c_elegans_swissprot_sequence_from_biomart(c_elegans_id_WB)
            if canonical_protein_sequence:
                return canonical_protein_sequence
        return None

    def extract_Human_gene_description_from_ENTREZ(gene_id):
        """ This function is extract Human gene description from ENTREZ."""
        Entrez.email = "sapir2@mail.tau.ac.il"
        print("extract_Human_gene_description_from_ENTREZ")
        print(gene_id)
        if not gene_id or gene_id == [None]:
            return "No description was found"
        request = Entrez.epost("gene", id=gene_id)
        result = Entrez.read(request)
        webEnv = result["WebEnv"]
        queryKey = result["QueryKey"]
        data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
        annotations = Entrez.read(data)
        description = annotations['DocumentSummarySet']['DocumentSummary'][0]['Summary']
        print(description)
        if description == None or description == "":
            return "no description was found"
        return (annotations['DocumentSummarySet']['DocumentSummary'][0]['Summary'])

# human_gene_name = "CAPZA1"
# c_elegans_gene_name = "cap-1"
# human_seq = BioPython.get_human_protein_from_uniprot_by_gene_name(human_gene_name)
# print("human seq:", human_seq)
# c_elegans_seq = BioPython().get_aa_seq_by_c_elegans_gene_name(c_elegans_gene_name)
# print("c-elegans seq", c_elegans_seq)

# print(BioPython.get_aa_seq_by_c_elegans_gene_name("repo-1"))
#
# print(BioPython.get_human_protein_from_uniprot_by_gene_name('TCP1'))
# print(BioPython.get_human_protein_from_uniprot_by_gene_name('DDR1'))
# print(BioPython.get_ncbi_id(Ensembl.get_gene_id_by_gene_name('TCP1', "Human")))
# print(BioPython.get_ncbi_id(Ensembl.get_gene_id_by_gene_name('DDR1', "Human")))
