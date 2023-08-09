import random
import re
from statistics import mean
from itertools import islice

from Code.CRISPR.NamedTuples.CodonData import CodonData
from Code.CRISPR.NamedTuples.CodonMutation import CodonMutation
from Code.CRISPR.NamedTuples.PointMutation import PointMutation
from Code.CRISPR.NamedTuples.MutateSection import MutateSection
from Code.CRISPR.NamedTuples.RestrictionMutation import RestrictionMutation
from Code.CRISPR.NamedTuples.RestrictionSite import RestrictionSite
from Code.CRISPR.NamedTuples.Result import Result
from Code.CRISPR.NamedTuples.SequenceSites import SequenceSites
from Code.CRISPR.Enum.AminoAcid import AminoAcid
from Code.CRISPR.Enum.DNASection import DNASection
from Code.CRISPR.Enum.MutationDirection import MutationDirection
from Code.CRISPR.Enum.RestrictionSiteType import RestrictionSiteType
from Code.Files.FileReader import FileReader
from Code.Http.HttpRequester import HttpRequester
from Code.Utils.BioPython import BioPython
from Code.Utils.Ensembl import Ensembl


class CrisprPlanner:
    amino_acid_dic = {'A': {'GCT', 'GCC', 'GCA', 'GCG'}, 'C': {'TGT', 'TGC'}, 'D': {'GAT', 'GAC'}, 'E': {'GAA', 'GAG'},
                      'F': {'TTT', 'TTC'}, 'G': {'GGT', 'GGC', 'GGA', 'GGG'}, 'H': {'CAT', 'CAC'},
                      'I': {'ATT', 'ATC', 'ATA'}, 'K': {'AAA', 'AAG'}, 'L': {'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'},
                      'M': {'ATG'}, 'N': {'AAT', 'AAC'}, 'P': {'CCT', 'CCC', 'CCA', 'CCG'}, 'Q': {'CAA', 'CAG'},
                      'R': {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'}, 'S': {'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'},
                      'T': {'ACT', 'ACC', 'ACA', 'ACG'}, 'V': {'GTT', 'GTC', 'GTA', 'GTG'}, 'W': {'TGG'},
                      'Y': {'TAT', 'TAC'}, 'STOP': {'TAA', 'TAG', 'TGA'}}

    amino_acid_dic_with_names = {AminoAcid.ALANINE: {'GCT', 'GCC', 'GCA', 'GCG'}, AminoAcid.CYSTEINE: {'TGT', 'TGC'}, AminoAcid.ASPARTIC_ACID: {'GAT', 'GAC'}, AminoAcid.GLUTAMIC_ACID: {'GAA', 'GAG'},
                      AminoAcid.PHENYLALANINE: {'TTT', 'TTC'}, AminoAcid.GLYCINE: {'GGT', 'GGC', 'GGA', 'GGG'}, AminoAcid.HISTIDINE: {'CAT', 'CAC'},
                      AminoAcid.ISOLEUCINE: {'ATT', 'ATC', 'ATA'}, AminoAcid.LYSINE: {'AAA', 'AAG'}, AminoAcid.LEUCINE: {'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'},
                      AminoAcid.METHIONINE: {'ATG'}, AminoAcid.ASPARAGINE: {'AAT', 'AAC'}, AminoAcid.PROLINE: {'CCT', 'CCC', 'CCA', 'CCG'}, AminoAcid.GLUTAMINE: {'CAA', 'CAG'},
                      AminoAcid.ARGININE: {'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'}, AminoAcid.SERINE: {'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'},
                      AminoAcid.THREONINE: {'ACT', 'ACC', 'ACA', 'ACG'}, AminoAcid.VALINE: {'GTT', 'GTC', 'GTA', 'GTG'}, AminoAcid.TRYPTOPHAN: {'TGG'},
                      AminoAcid.TYROSINE: {'TAT', 'TAC'}}

    codon_dic = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TGT': 'C', 'TGC': 'C', 'GAT': 'D', 'GAC': 'D',
                 'GAA': 'E', 'GAG': 'E', 'TTT': 'F', 'TTC': 'F', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
                 'CAT': 'H', 'CAC': 'H', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K', 'TTA': 'L',
                 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATG': 'M', 'AAT': 'N', 'AAC': 'N',
                 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
                 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                 'AGT': 'S', 'AGC': 'S', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GTT': 'V', 'GTC': 'V',
                 'GTA': 'V', 'GTG': 'V', 'TGG': 'W', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
                 'TGA': 'STOP'}

    variables = {'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
                 'H': ('A', 'C', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                 'N': ('A', 'C', 'G', 'T'), 'R': ('A', 'G'), 'S': ('C', 'G'),
                 'V': ('A', 'C', 'G'), 'W': ('A', 'T'), 'Y': ('C', 'T')}

    codon_usage = {'UUU': 23.3, 'UCU': 16.7, 'UAU': 17.5, 'UGU': 11.2, 'UUC': 23.9, 'UCC': 10.6, 'UAC': 13.7,
                   'UGC': 9.1,
                   'UUA': 9.8, 'UCA': 20.6, 'UAA': 1.6, 'UGA': 1.4, 'UUG': 20.0, 'UCG': 12.2, 'UAG': 0.6, 'UGG': 11.1,
                   'CUU': 21.2, 'CCU': 8.8, 'CAU': 14.1, 'CGU': 11.2, 'CUC': 14.8, 'CCC': 4.4, 'CAC': 9.2, 'CGC': 5.1,
                   'CUA': 7.9, 'CCA': 26.1, 'CAA': 27.4, 'CGA': 12.1, 'CUG': 12.1, 'CCG': 9.7, 'CAG': 14.4, 'CGG': 4.7,
                   'AUU': 32.2, 'ACU': 18.9, 'AAU': 30.2, 'AGU': 12.1, 'AUC': 18.9, 'ACC': 10.4, 'AAC': 18.3,
                   'AGC': 8.4,
                   'AUA': 9.5, 'ACA': 20.0, 'AAA': 37.5, 'AGA': 15.4, 'AUG': 26.1, 'ACG': 8.9, 'AAG': 25.8, 'AGG': 4.0,
                   'GUU': 24.1, 'GCU': 22.4, 'GAU': 35.8, 'GGU': 10.9, 'GUC': 13.6, 'GCC': 12.6, 'GAC': 17.1,
                   'GGC': 6.7,
                   'GUA': 9.8, 'GCA': 19.8, 'GAA': 40.8, 'GGA': 31.7, 'GUG': 14.3, 'GCG': 8.2, 'GAG': 24.5, 'GGG': 4.4}

    def __init__(self, gene_name, aa_mutation_site , end ,given_isoform, type, sense_strand: str = "", amino_acid_sequence: str = "",
                favourite_enzymes_names=None, max_results=1):
        print("init")
        print(gene_name, aa_mutation_site,given_isoform, end, type, sense_strand, amino_acid_sequence,
                favourite_enzymes_names, max_results)
        print("gene_name")
        self.gene_name = gene_name
        print(self.gene_name)
        print("type")
        self.type = type
        print(self.type)
        print("gene_id")
        self.gene_id = Ensembl.get_c_elegans_gene_id_by_gene_name(gene_name)
        print(self.gene_id)
        self.result = Result(None,None, None, None, False, None, [], [], [], [], None, None, 0, None, 0, [], [], [], [])
        print("given_isoform")
        self.given_isoform= given_isoform if given_isoform else None
        print(self.given_isoform)
        self.from_aa = None
        self.to_aa = None
        print("from_aa: ",self.from_aa, "to_aa: ",self.to_aa)
        print("sense_strand")
        transcript_result = HttpRequester.get_transcript(self.gene_id, self.given_isoform,self.type, self.result)
        self.sense_strand = sense_strand if sense_strand else transcript_result[0]
        print(self.sense_strand)
        print("help_isoform_id")
        self.help_isoform_id= given_isoform if given_isoform else transcript_result[1]
        print(self.help_isoform_id)
        self.anti_sense_strand = CrisprPlanner.get_complementary_sequence(self.sense_strand)
        print("amino_acid_sequence")
        self.amino_acid_sequence = amino_acid_sequence if amino_acid_sequence else \
            BioPython().get_aa_seq_by_c_elegans_gene_name(gene_name , self.help_isoform_id)
        print(self.amino_acid_sequence)
        print("lennnnn: ", len(self.amino_acid_sequence))
        #self.end = False
        self.end = end
        print("self.amino_acid_site")
        self.amino_acid_site = len(self.amino_acid_sequence) if self.end else aa_mutation_site
        print("site   + + end   " ,self.end , self.amino_acid_site)
        #self.amino_acid_site = aa_mutation_site
        self.codon = None
        self.sense_mutation_site = -1
        self.anti_sense_mutation_site = -1
        self.ssODN_mutation_codon_start = -1
        self.mutated_strand = None
        self.mutation_direction = None
        self.dsb_merged = False
        self.ssODN_vs_crRNA_strand = 1
        self.DSB = None
        self.pam_sites = None
        self.mutation_zone = None
        self.ssODN_direction = None
        self.mutated_sites = []
        self.codon_mutations = []
        self.reattachment_mutations = []
        self.restriction_site_mutations = []
        self.restriction_enzymes = self.get_favourite_restriction_enzymes(favourite_enzymes_names)
        self.max_results = max_results
        self.maybe_problematic=(-1,"")
        self.problematic=False
        self.reverse_linker = None
        self.forward_linker = None
        self.primer_reverse = None
        self.primer_forward = None
        self.insertion = None
        self.intronic_codon = False
        self.type_of_intronic_codon = False
        self.tag_for_intronic = False         ## 1 - the guides locations are OK, 2 - cant do this, throw an error , 3 - change the guides locations
        self.length_of_intron = 0

    def get_relevant_strand(self, direction):
        print("def get_relevant_strand")
        if direction > 0:
            return self.sense_strand
        else:
            return self.anti_sense_strand

    def plan_my_crispr(self,
                       from_aa: AminoAcid,
                       to_aa: AminoAcid,
                       crrna: str = None,
                       crrna_strand: int = 1,
                       check_consistency: bool = False,
                       window_size: int = 30,
                       PAM_size: int = 3):

        print("def: plan_my_crispr")
        self.type = "point"

        if from_aa == to_aa:
            e = "Your request to change " + from_aa.name + " to " + to_aa.name + " is invalid."
            return self.result, e
        self.from_aa = from_aa
        self.to_aa = to_aa
        self.result.query_data = (self.gene_name, self.from_aa, self.amino_acid_site, self.to_aa)


        possible_error = self.initiate_crispr(check_consistency)
        if possible_error:
            print(possible_error)
            return self.result, possible_error

        sense_crrnas, anti_sense_crrnas = self.find_crrnas(window_size, PAM_size)
        print("line to debug")
        print("sense crrnas:", sense_crrnas, sep="\n")
        print("anti-sense crrnas:", anti_sense_crrnas, sep="\n")
        self.result.sense_crrnas = sense_crrnas
        self.result.anti_sense_crrnas = anti_sense_crrnas
        self.result.isoform =self.help_isoform_id

        if not crrna:  # find crrnas and return them
            if not sense_crrnas and not anti_sense_crrnas:
                return self.result, "Couldn't find crRNA sequences"
            sense_crrna_values_dic, sense_error = HttpRequester.check_crRNA([crrna[0] for crrna in sense_crrnas])
            self.result.sense_crrnas_results = list(sense_crrna_values_dic.items())
            anti_sense_crrna_values_dic, anti_sense_error = HttpRequester.check_crRNA([crrna[0] for crrna in
                                                                                       anti_sense_crrnas])
            self.result.anti_sense_crrnas_results = list(anti_sense_crrna_values_dic.items())
            print("crrans results:", self.result.anti_sense_crrnas_results, self.result.sense_crrnas_results)
            self.result.crrnas_results_status = (sense_error, anti_sense_error)
            return self.result, None
        else:
            chosen_crrna = self.get_chosen_crrna(crrna, crrna_strand)
            if not chosen_crrna:
                return self.result, "crRNA sequence was not found among optional crRNAs, please check if you chose " \
                                    "the right strand.\nAlso Please note that the system is case-sensitive, thus in " \
                                    "copying down the relevant sequence, keep introns in lower case and exons in " \
                                    "upper case."
            self.result.crRNA = chosen_crrna
            self.result.crRNA_strand = crrna_strand

        # now we have our cr_rna
        strand_str = "anti-" if crrna_strand < 0 else ""
        print("chosen crRNA:", self.result.crRNA, "The relevant strand is", strand_str + "sense")

        self.complete_fields(crrna_strand, to_aa, self.result.crRNA)

        # first - when at last we try to find restriction site
        first_strand = self.mutated_strand
        first_mutated_sites = self.mutated_sites

        for restriction_site_type in RestrictionSiteType:
            print("Restriction site type:", restriction_site_type)
            self.mutated_strand = first_strand
            self.mutated_sites = first_mutated_sites[:]
            is_successful = self.insert_mutations(self.ssODN_mutation_codon_start, restriction_site_type)
            print("is_successful " , is_successful ,"170")
            if is_successful:
                print("successful!")
                return self.result, None
            # else
            print("sapir5")
            print("problematic=" ,self.problematic)
            if self.problematic==True:
                e = "Reattachment section could not be mutated, try a different crRNA."
                print(e)
                return self.result, e
            print("reattachment_mutations:" , self.reattachment_mutations)
            print("dsb_merged:", self.dsb_merged)
            if not self.reattachment_mutations:
                if not self.dsb_merged:
                    # couldn't find a way to mutate reattachment site
                    e = "Reattachment section could not be mutated, try a different crRNA."
                    print(e)
                    return self.result, e
                else:
                    print("PAM site could not be mutated. Moving on to trying to mutate crRNA")
                    return self.dsb_merged_scheme(first_strand, crrna_strand, chosen_crrna)
        e = "We could not find any restriction sites to insert or to remove."
        return self.result, e


    def plan_my_crispr_insertion(self,
                       from_aa: AminoAcid,
                       to_aa: AminoAcid,
                       crrna: str = None,
                       crrna_strand: int = 1,
                       check_consistency: bool = False,
                       window_size: int = 30,
                       PAM_size: int = 3,
                       reverse_linker: str = None,
                       forward_linker: str = None,
                       primer_reverse: str = None,
                       primer_forward: str = None,
                       insertion: str = None,
                       #end = False
    ):

        aa_names_dict = {'A': "ALANINE", 'R': "ARGININE", 'N': "ASPARAGINE", 'D': "ASPARTIC_ACID", 'C': "CYSTEINE",
                         'E': "GLUTAMIC_ACID", 'Q': "GLUTAMINE", 'G': "GLYCINE", 'H': "HISTIDINE", 'I': "ISOLEUCINE",
                         'L': "LEUCINE",
                         'K': "LYSINE", 'M': "METHIONINE", 'F': "PHENYLALANINE", 'P': "PROLINE", 'S': "SERINE",
                         'T': "THREONINE",
                         'W': "TRYPTOPHAN", 'Y': "TYROSINE", 'V': "VALINE"}

        #self.from_aa = from_aa
        self.from_aa = AminoAcid[(aa_names_dict[BioPython().get_aa_seq_by_c_elegans_gene_name(self.gene_name,self.help_isoform_id)[self.amino_acid_site - 1]])]
        self.to_aa = to_aa
        self.result.query_data = (self.gene_name, self.from_aa, self.amino_acid_site, self.to_aa)
        self.reverse_linker = reverse_linker
        self.forward_linker = forward_linker
        self.primer_reverse = primer_reverse
        self.primer_forward = primer_forward
        self.insertion = insertion
        self.type = "insertion"
        #self.end = end

        print("func: plan_my_crispr_insertion")
        print(self.end)

        possible_error = self.initiate_crispr(check_consistency)
        if possible_error:
            print(possible_error)
            return self.result, possible_error

        char_list = ["a","c","g","t","A","C","G","T",""]
        string = reverse_linker+primer_reverse+primer_forward+forward_linker
        matched_list = [characters in char_list for characters in string]
        string_contains_chars = all(matched_list)
        if string_contains_chars == False:
            return self.result, "Check that you entered valid sequences for linkers and primers."
        if ((primer_forward == "") and (primer_reverse!="")) or ((primer_forward != "") and (primer_reverse=="")):
            return self.result, "Check that you entered valid sequences for primers."
        if primer_reverse=="" and primer_forward=="":
            self.result.without_primers = True
        sense_crrnas, anti_sense_crrnas = self.find_crrnas(window_size, PAM_size)
        print("line to debug 2")
        print("sense crrnas:", sense_crrnas, sep="\n")
        print("anti-sense crrnas:", anti_sense_crrnas, sep="\n")
        self.result.sense_crrnas = sense_crrnas
        self.result.anti_sense_crrnas = anti_sense_crrnas
        self.result.isoform =self.help_isoform_id

        if not crrna:  # find crrnas and return them
            if not sense_crrnas and not anti_sense_crrnas:
                return self.result, "Couldn't find crRNA sequences"
            sense_crrna_values_dic, sense_error = HttpRequester.check_crRNA([crrna[0] for crrna in sense_crrnas])
            self.result.sense_crrnas_results = list(sense_crrna_values_dic.items())
            anti_sense_crrna_values_dic, anti_sense_error = HttpRequester.check_crRNA([crrna[0] for crrna in
                                                                                       anti_sense_crrnas])
            self.result.anti_sense_crrnas_results = list(anti_sense_crrna_values_dic.items())
            print("crrans results:", self.result.anti_sense_crrnas_results, self.result.sense_crrnas_results)
            self.result.crrnas_results_status = (sense_error, anti_sense_error)
            return self.result, None



        else:
            chosen_crrna = self.get_chosen_crrna(crrna, crrna_strand)
            if not chosen_crrna:
                return self.result, "crRNA sequence was not found among optional crRNAs, please check if you chose " \
                                    "the right strand.\nAlso Please note that the system is case-sensitive, thus in " \
                                    "copying down the relevant sequence, keep introns in lower case and exons in " \
                                    "upper case."
            self.result.crRNA = chosen_crrna
            self.result.crRNA_strand = crrna_strand

        # now we have our cr_rna
        strand_str = "anti-" if crrna_strand < 0 else ""
        print("chosen crRNA:", self.result.crRNA, "The relevant strand is", strand_str + "sense")

        self.complete_fields(crrna_strand, to_aa, self.result.crRNA)


        print("pam_sites:",self.pam_sites)
        print("DSB:",self.DSB)
        print("sense_mutation_site:",self.sense_mutation_site)
        print("anti_sense_mutation_site:",self.anti_sense_mutation_site)
        print("ssODN_direction:",self.ssODN_direction)
        print("ssODN_mutation_codon_start:",self.ssODN_mutation_codon_start)
        print("try")
        print(self.end)
        if self.end == True:
            print("end, true")
            self.sense_mutation_site = self.sense_mutation_site +3
            self.anti_sense_mutation_site = self.anti_sense_mutation_site -3
            if self.ssODN_mutation_codon_start == self.sense_mutation_site:
                self.ssODN_mutation_codon_start = self.ssODN_mutation_codon_start+3
            else:
                self.ssODN_mutation_codon_start = self.ssODN_mutation_codon_start - 3
        print("pam_sites:",self.pam_sites)
        print("DSB:",self.DSB)
        print("sense_mutation_site:",self.sense_mutation_site)
        print("anti_sense_mutation_site:",self.anti_sense_mutation_site)
        print("ssODN_direction:",self.ssODN_direction)
        print("ssODN_mutation_codon_start:",self.ssODN_mutation_codon_start)
        print("try")

        print ("the crna strand we are looking at is:", crrna_strand)
        self.result.insertion = self.insertion
        
        #if self.ssODN_mutation_codon_start == self.sense_mutation_site:   ## looking the sense strand
            ## if crrna_strand == 1 and if crrna_strand == -1 we want to
        if crrna_strand == 1:
            dsb = self.result.crRNA[1].end - 3
            print("dsb: " ,dsb)
            #if self.DSB >= self.sense_mutation_site:                     ##  the cut is right to the insertion
            if self.sense_mutation_site < self.result.crRNA[1].end -3:
                print("case 4")
                print(dsb)
                cant_destroy_pam = False
                # print(self.sense_mutation_site,self.result.crRNA[1].start+4)
                pam_start = dsb + 4
                pam_end = dsb + 6
                if self.sense_mutation_site < self.result.crRNA[1].start + 4:
                    print("need to insert new mutations!")
                    print(self.sense_strand[pam_start:pam_end + 1])
                    if "GG" in self.sense_strand[pam_start:pam_end + 1]:
                        print("we have an exon to deal with")
                        self.check_speical_states_to_destroy_pam()
                        print(self.maybe_problematic)
                        if (self.maybe_problematic[0] == 0):
                            if (self.maybe_problematic[1] == "CGG" or self.maybe_problematic[1] == "AGG" or
                                    self.maybe_problematic[1] == "GGG"):
                                self.sense_strand = self.sense_strand[:pam_start + 2] + "A" + self.sense_strand[pam_end + 1:]
                                # self.sense_strand[self.pam_sites.end + 1] = "A"    ## replace the last G with A
                            elif (self.maybe_problematic[1] == "TGG" or self.maybe_problematic[1] == "GGG"):
                                cant_destroy_pam = True     # check
                        if (self.maybe_problematic[0] == 1 and self.maybe_problematic[1][3:5] == "GG"):
                            cant_destroy_pam = True
                        if (self.maybe_problematic[0] == 2):
                            if (self.maybe_problematic[1][0:4] == "ATGG" or self.maybe_problematic[1][0:4] == "TGGG"):
                                cant_destroy_pam = True
                            else:          #check
                                self.sense_strand = self.sense_strand[:pam_start + 1] + "A" + self.sense_strand[pam_end:]
                                # self.sense_strand[self.pam_sites.start+1] = "A"   ## replace the middle G of pam to A
                        print("cant_destroy_pam ", cant_destroy_pam)
                    else:  ## the gg are in intronic region
                        self.sense_strand = self.sense_strand[:pam_start] + "aa" + self.sense_strand[pam_start + 2:]
                        # self.sense_strand[self.pam_sites.start] = "a"    ## the gg are in intronic region

                    if cant_destroy_pam == True:
                        print("create_new_guide_with_mutations")
                        self.create_new_guide_with_mutations(pam_start,1)

                #print(self.result.crRNA[1])
                #print(self.result.crRNA[1].start)
                HA1 = self.sense_strand[self.sense_mutation_site-35+3:self.sense_mutation_site+3]
                print('HA1',HA1)
                complete = self.sense_strand[self.sense_mutation_site+3:dsb+1]
                print("complete" ,complete)
                HA2 = self.sense_strand[dsb+1:dsb+1+35]
                print("HA2",HA2)
                forward_result = HA1+self.forward_linker+self.primer_forward
                reverse_result = self.create_complementary_and_reverse_string(HA2)+self.create_complementary_and_reverse_string(complete)+self.reverse_linker+self.primer_reverse
                print("forward_result",forward_result)
                print("reverse_result",reverse_result)
                self.result.side = True

            if self.sense_mutation_site >= self.result.crRNA[1].end -3:                    ##  the cut is right to the insertion
                print("case 3")
                cant_destroy_pam = False
                # print(self.sense_mutation_site,self.result.crRNA[1].start+4)
                pam_start = dsb + 4
                pam_end = dsb + 7
                print("pam: ",pam_start,pam_end)
                print(self.sense_strand[pam_start])
                print(self.sense_strand[pam_end])
                if pam_end < self.sense_mutation_site:
                    print("need to insert new mutations!")
                    print(self.sense_strand[pam_start:pam_end])
                    if "GG" in self.sense_strand[pam_start:pam_end]:
                        print("we have an exon to deal with")
                        self.check_speical_states_to_destroy_pam()
                        print(self.maybe_problematic)
                        if (self.maybe_problematic[0] == 0):
                            if (self.maybe_problematic[1] == "CGG" or self.maybe_problematic[1] == "AGG" or
                                    self.maybe_problematic[1] == "GGG"):
                                self.sense_strand = self.sense_strand[:pam_start + 2] + "A" + self.sense_strand[pam_end:]
                                # self.sense_strand[self.pam_sites.end + 1] = "A"    ## replace the last G with A
                            elif (self.maybe_problematic[1] == "TGG" or self.maybe_problematic[1] == "GGG"):
                                cant_destroy_pam = True
                        if (self.maybe_problematic[0] == 1 and self.maybe_problematic[1][3:5] == "GG"):
                            cant_destroy_pam = True
                        if (self.maybe_problematic[0] == 2):
                            if (self.maybe_problematic[1][0:4] == "ATGG" or self.maybe_problematic[1][0:4] == "TGGG"):
                                cant_destroy_pam = True
                            else:
                                self.sense_strand = self.sense_strand[:pam_start + 1] + "A" + self.sense_strand[pam_end-1:]
                                # self.sense_strand[self.pam_sites.start+1] = "A"   ## replace the middle G of pam to A
                        print("cant_destroy_pam ", cant_destroy_pam)
                    else:  ## the gg are in intronic region
                        self.sense_strand = self.sense_strand[:pam_start+1] + "aa" + self.sense_strand[pam_start + 3:]
                        # self.sense_strand[self.pam_sites.start] = "a"    ## the gg are in intronic region

                    if cant_destroy_pam == True:
                        print("create_new_guide_with_mutations")
                        self.create_new_guide_with_mutations(pam_start,1)

                HA1 = self.sense_strand[dsb-35+1:dsb+1]
                print('HA1',HA1)
                complete = self.sense_strand[dsb+1:self.sense_mutation_site+3]
                print("complete", complete)
                HA2 = self.sense_strand[self.sense_mutation_site+3:self.sense_mutation_site+35+3]
                print("HA2",HA2)
                self.result.side = False

            self.result.forward_seq = HA1+self.forward_linker+self.primer_forward
            self.result.reverse_seq = self.create_complementary_and_reverse_string(HA2)+self.create_complementary_and_reverse_string(complete)+self.reverse_linker+self.primer_reverse
            self.result.insertion = self.insertion
            self.result.HA1_forward = HA1
            self.result.HA2_forward = HA2
            self.result.linker_forward = self.forward_linker
            self.result.linker_reverse_forward = self.create_complementary_and_reverse_string(reverse_linker)
            self.result.primer_forward = self.primer_forward
            self.result.primer_reverse_forward = self.create_complementary_and_reverse_string(primer_reverse)
            self.result.complete = complete
            #self.result.side = False


        #if self.ssODN_mutation_codon_start == self.anti_sense_mutation_site:   ## looking the anti-sense strand
        if crrna_strand == -1:
            dsb = self.result.crRNA[1].end - 3
            print(dsb)
            if self.anti_sense_mutation_site < self.result.crRNA[1].end - 3:
            #if self.DSB >= self.anti_sense_mutation_site+3:                     ##  the cut is right to the insertion
                print("case 1")
                print ("dsb" , dsb)
                cant_destroy_pam = False
                # print(self.sense_mutation_site,self.result.crRNA[1].start+4)
                pam_start = dsb + 4
                pam_end = dsb + 7
                print("pam: ", pam_start, pam_end)
                print(self.anti_sense_strand[pam_start])
                print(self.anti_sense_strand[pam_end])
                print(self.anti_sense_mutation_site)
                print(self.result.crRNA[1].start)
                if self.anti_sense_mutation_site < self.result.crRNA[1].start + 4:
                    print("need to insert new mutations!")
                    print(self.anti_sense_strand[pam_start:pam_end])
                    if "GG" in self.anti_sense_strand[pam_start:pam_end]:
                        print("we have an exon to deal with")
                        self.check_speical_states_to_destroy_pam()
                        print(self.maybe_problematic)
                        #self.maybe_problematic = (self.maybe_problematic[0],self.create_complementary_and_reverse_string(self.maybe_problematic[1]))
                        print(self.maybe_problematic)
                        if (self.maybe_problematic[0] == 0):
                            if (self.maybe_problematic[1] == "CCC" or self.maybe_problematic[1] == "CCG" or
                                self.maybe_problematic[1] == "CCT" or self.maybe_problematic[1] == "CCA"):
                                cant_destroy_pam = True
                        if (self.maybe_problematic[0] == 1):
                            print("self.maybe_problematic[0] == 1")     ## check
                            print(self.anti_sense_strand[pam_start + 1-10:pam_start + 2])
                            print(self.anti_sense_strand[pam_end:pam_end+10])
                            self.anti_sense_strand = self.anti_sense_strand[
                                                     :pam_start + 2] + "A" + self.anti_sense_strand[
                                                                             pam_end:]
                        if (self.maybe_problematic[0] == 2):
                            print("self.maybe_problematic[0] == 2)")
                            self.anti_sense_strand = self.anti_sense_strand[:pam_start + 1] + "A" + self.anti_sense_strand[
                                                                                              pam_end - 1:]
                                # self.sense_strand[self.pam_sites.start+1] = "A"   ## replace the middle G of pam to A
                        print("cant_destroy_pam ", cant_destroy_pam)
                    else:  ## the gg are in intronic region
                        self.anti_sense_strand = self.anti_sense_strand[:pam_start + 1] + "tt" + self.anti_sense_strand[pam_start + 3:]
                        # self.sense_strand[self.pam_sites.start] = "a"    ## the gg are in intronic region

                    if cant_destroy_pam == True:
                        print("create_new_guide_with_mutations")
                        self.create_new_guide_with_mutations(pam_start,-1)


                HA1 = self.anti_sense_strand[dsb+1:dsb+1+35]
                print('HA1',HA1)
                complete = self.anti_sense_strand[self.anti_sense_mutation_site:dsb+1]
                print("complete:",complete)
                HA2 = self.anti_sense_strand[self.anti_sense_mutation_site-35:self.anti_sense_mutation_site]
                print("HA2",HA2)
                self.result.side = False

            if self.anti_sense_mutation_site >= self.result.crRNA[1].end - 3:                     ##  the cut is right to the insertion
                print("case 2")
                print ("dsb" , dsb)
                cant_destroy_pam = False
                # print(self.sense_mutation_site,self.result.crRNA[1].start+4)
                pam_start = dsb + 4
                pam_end = dsb + 7
                print("pam: ", pam_start, pam_end)
                print(self.anti_sense_strand[pam_start])
                print(self.anti_sense_strand[pam_end])
                print(self.anti_sense_mutation_site)
                print(self.result.crRNA[1].start)
                if pam_end < self.anti_sense_mutation_site:
                    print("need to insert new mutations!")
                    print(self.anti_sense_strand[pam_start:pam_end])
                    if "GG" in self.anti_sense_strand[pam_start:pam_end]:
                        print("we have an exon to deal with")
                        self.check_speical_states_to_destroy_pam()
                        print(self.maybe_problematic)
                        #self.maybe_problematic = (self.maybe_problematic[0],self.create_complementary_and_reverse_string(self.maybe_problematic[1]))
                        print(self.maybe_problematic)
                        if (self.maybe_problematic[0] == 0):
                            if (self.maybe_problematic[1] == "CCC" or self.maybe_problematic[1] == "CCG" or
                                self.maybe_problematic[1] == "CCT" or self.maybe_problematic[1] == "CCA"):
                                cant_destroy_pam = True
                        if (self.maybe_problematic[0] == 1):
                            print("self.maybe_problematic[0] == 1")     ## check
                            print(self.anti_sense_strand[pam_start + 1-10:pam_start + 2])
                            print(self.anti_sense_strand[pam_end:pam_end+10])
                            self.anti_sense_strand = self.anti_sense_strand[
                                                     :pam_start + 2] + "A" + self.anti_sense_strand[
                                                                             pam_end:]
                        if (self.maybe_problematic[0] == 2):
                            print("self.maybe_problematic[0] == 2)")
                            self.anti_sense_strand = self.anti_sense_strand[:pam_start + 1] + "A" + self.anti_sense_strand[
                                                                                              pam_end - 1:]
                                # self.sense_strand[self.pam_sites.start+1] = "A"   ## replace the middle G of pam to A
                        print("cant_destroy_pam ", cant_destroy_pam)
                    else:  ## the gg are in intronic region
                        self.anti_sense_strand = self.anti_sense_strand[:pam_start + 1] + "tt" + self.anti_sense_strand[pam_start + 3:]
                        # self.sense_strand[self.pam_sites.start] = "a"    ## the gg are in intronic region

                    if cant_destroy_pam == True:
                        print("create_new_guide_with_mutations")
                        self.create_new_guide_with_mutations(pam_start,-1)

                HA1 = self.anti_sense_strand[self.anti_sense_mutation_site:self.anti_sense_mutation_site+35]
                print('HA1',HA1)
                complete = self.anti_sense_strand[dsb+1:self.anti_sense_mutation_site]
                print("complete:",complete )
                HA2 = self.anti_sense_strand[dsb-35+1:dsb+1]
                print("HA2",HA2)
                self.result.side = True

            self.result.forward_seq = self.create_complementary_and_reverse_string(HA1)+self.create_complementary_and_reverse_string(complete)+self.forward_linker+self.primer_forward
            self.result.reverse_seq = HA2+self.reverse_linker+self.primer_reverse
            self.result.insertion = self.insertion
            self.result.HA1_forward = self.create_complementary_and_reverse_string(HA1)
            self.result.HA2_forward = self.create_complementary_and_reverse_string(HA2)
            self.result.linker_forward = self.forward_linker
            self.result.linker_reverse_forward = self.create_complementary_and_reverse_string(reverse_linker)
            self.result.primer_forward = self.primer_forward
            self.result.primer_reverse_forward = self.primer_reverse
            self.result.complete = self.create_complementary_and_reverse_string(complete)
            #self.result.side = True

        return self.result ,""


    # Takes care of the possibility where the mutation codon is merged with dsb and PAM site cannot be mutated
    def dsb_merged_scheme(self, first_strand, crrna_strand, chosen_crrna):
        print("def dsb_merged_scheme")
        first_strand = self.get_complementary_sequence(first_strand)
        self.ssODN_vs_crRNA_strand = -1 * self.ssODN_vs_crRNA_strand
        self.complete_fields(crrna_strand, self.to_aa, chosen_crrna)
        for r_site_type in RestrictionSiteType:
            self.mutated_strand = first_strand
            self.mutated_sites = []
            is_successful = self.insert_mutations(self.ssODN_mutation_codon_start, r_site_type)
            print("is_successful ", is_successful, "305")
            if is_successful:
                return self.result, None
            if not self.reattachment_mutations:
                # couldn't find a way to mutate reattachment site
                e = "Reattachment section could not be mutated"
                print(e)
                return self.result, e
        e = "We could not find any restriction sites to insert or to remove."
        return self.result, e

    def create_new_guide_with_mutations(self,pam_start,strand):
        print("def: create_new_guide_with_mutations")
        guide = self.result.crRNA[0]
        case = self.maybe_problematic[0]
        if strand == 1:
            if case == 0:
                lst = [pam_start-21,pam_start-1]
            if case == 1:
                lst = [pam_start-20,pam_start]
            if case == 2:
                lst = [pam_start-19,pam_start-2]
            print(self.sense_strand[lst[0]:lst[1]])
            start = lst[0]
            end = lst[1]
            text = self.sense_strand
            print(lst)
            print(guide)
            print(self.sense_strand[lst[0]:lst[1]])
            print(self.anti_sense_strand[lst[0]+1:lst[1]+1])
            print(case)
            print(lst[0],lst[1])
            print(strand)
            for i in range(start,end,3):
                codon = text[i:i+3]
                print(codon)
                if codon != "ATG" and codon != "TGA" and codon != "TGG" and codon != "ATA":
                    print("if 1")
                    #print(self.sense_strand[:i + 1])
                    #print(self.sense_strand[i + 2:])
                    if codon.endswith("G"):
                        text = text[:i + 2] + "A" + text[i + 3:]
                    if codon.endswith("A"):
                        text = text[:i + 2] + "G" + text[i + 3:]
                    if codon.endswith("C"):
                        text = text[:i + 2] + "T" + text[i + 3:]
                    if codon.endswith("T"):
                        text = text[:i + 2] + "C" + text[i + 3:]
                if codon == "ATA":
                    print("if 2")
                    text = text[:i + 1] + "C" + text[i + 2:]
            self.sense_strand = text
        if strand == -1:
            print(strand)
            print(case)
            if case == 0:
                lst = [pam_start-22,pam_start-1]
            if case == 1:
                lst = [pam_start-20,pam_start]
            if case == 2:
                lst = [pam_start-19,pam_start-2]
            print(lst)
            print(self.anti_sense_strand[lst[0]+1:lst[1]+1])
            start = lst[0]+1
            end = lst[1]+1
            text = self.anti_sense_strand
            for i in range(start,end,3):
                codon = text[i:i+3]
                print(codon)
                if codon != "CAT" and codon != "CCA" and codon != "TCA" and codon != "TAT":
                    print("if 1")
                    #print(self.sense_strand[:i + 1])
                    #print(self.sense_strand[i + 2:])
                    if codon.startswith("G"):
                        text = text[:i] + "A" + text[i + 1:]
                    if codon.startswith("A"):
                        text = text[:i] + "G" + text[i + 1:]
                    if codon.startswith("C"):
                        text = text[:i] + "T" + text[i + 1:]
                    if codon.startswith("T"):
                        text = text[:i] + "C" + text[i + 1:]
                if codon == "TAT":
                    print("if 2")
                    text = text[:i] + "G" + text[i + 1:]
            self.anti_sense_strand = text
        return


    @staticmethod
    def starts_with_ATG(strand):
        print("def: starts_with_ATG")
        for i in range(len(strand)):
            if strand[i].islower():  # intron
                continue
            else:
                if strand[i:i + 3] == 'ATG':
                    return True
                else:
                    return False

    def choose_aa_codon(self):
        print("def: starts_with_ATG")
        codon = self.sense_strand[self.sense_mutation_site:self.sense_mutation_site + 3]
        print(codon)
        #print(self.sense_strand)
        #print(self.sense_mutation_site,self.sense_mutation_site + 3)
        lowercase_letters = ''.join(c for c in codon if c.islower())
        real_codon = ''.join(c for c in codon if c.isupper())
        #print(lowercase_letters)
        #print(real_codon)
        if len(lowercase_letters) > 0:
            print("there is a intron!")
            #print(codon[1])
            intron_len = 0
            new_exon_start = self.sense_mutation_site+3
            print("new_exon_start",new_exon_start)
            if codon[1].isupper() == True:
            #    new_exon_start = self.sense_mutation_site + 3
            #    print(new_exon_start)
                intron_len = 1
            else:
            #    new_exon_start = self.sense_mutation_site + 4
                intron_len = 2
            while True:
                if self.sense_strand[new_exon_start].islower() == True:
                    new_exon_start +=1
                    intron_len +=1
                else:
                    self.length_of_intron = intron_len
                    break
            print(new_exon_start)
            print("gg" , real_codon)
            print(intron_len)
            if len(real_codon) == 1:      ## the codon i splitten by 1/2
                self.type_of_intronic_codon = 1
                real_codon = real_codon+(self.sense_strand[new_exon_start:new_exon_start+2])
                #if_intronic_check_which_aa_need_to_change_case1(self)

            if len(real_codon) == 2:        ## the codon i splitten by 2/1
                #print(new_exon_start)
                real_codon = real_codon+(self.sense_strand[new_exon_start])
                self.type_of_intronic_codon = 2
            print(real_codon)
            self.codon = real_codon
            self.intronic_codon = True
        else:
            self.codon = codon


    def for_intronic_codon_check_new_codon_if_posibble(self):
        print("func = for_intronic_codon_check_new_codon_if_posibble")
        codon = self.codon
        print(codon)
        print("self.type_of_intronic_codon" , self.type_of_intronic_codon)
        if self.type_of_intronic_codon == 1:
            for nuc in ["A","C","T","G"]:
                one_change1 = nuc + codon[1:3]
                one_change2 = codon[0] + nuc + codon[2]
                one_change3 = codon[0:2] + nuc
                print(one_change1 ,one_change2 ,one_change3)
                print(self.amino_acid_dic_with_names[self.to_aa])
                if one_change1 in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 1
                    print(self.tag_for_intronic)
                    return (one_change1)
                if one_change2 in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 3
                    print(self.tag_for_intronic)
                    return (one_change2)
                if one_change3 in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 3
                    print(self.tag_for_intronic)
                    return (one_change3)
            for double in ["AA","CC","GG","TT","AG","AC","AT","GA","GT","GC","CA","CG","CT","TA","TC","TG"]:
                double_change = codon[0]+double
                if double_change in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 3
                    return (double_change)
        if self.type_of_intronic_codon == 2:
            for nuc in ["A","C","T","G"]:
                one_change1 = nuc + codon[1:2]
                one_change2 = codon[0] + nuc + codon[2]
                one_change3 = codon[0:1] + nuc
                if one_change1 in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 1
                    return (one_change1)
                if one_change2 in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 1
                    return (one_change2)
                if one_change3 in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 3
                    return (one_change3)
            for double in ["AA","CC","GG","TT","AG","AC","AT","GA","GT","GC","CA","CG","CT","TA","TC","TG"]:
                double_change = double+codon[2]
                if double_change in self.amino_acid_dic_with_names[self.to_aa]:
                    self.tag_for_intronic = 1
                    return (double_change)
        self.tag_for_intronic = 2
        return (None)


    # initiates the CRISPR planner. fills in the amino_acid_sequence if not already filled in, and finds the nt site in
    # sense strand and anti-sense strand
    def initiate_crispr(self, check_consistency: bool):
        print("def: initiate_crispr")
        if not self.sense_strand:
            # sense strand was not given and extraction has failed
            e = "Exception in initiate_crispr: no sense strand sequence was delivered or could be extracted, " \
                "please re-send the request with the needed sequence."
            return e
        # sense strand exists, let's check if it's the right direction:
        if not self.starts_with_ATG(self.sense_strand):
            if self.starts_with_ATG(self.anti_sense_strand):  # switched
                self.sense_strand, self.anti_sense_strand = self.anti_sense_strand, self.sense_strand
            else:
                e = "Exception in initiate_crispr: valid sense strand could not be found, " \
                    "please re-send the request with the needed sequence."
                return e
        #print("Gene's sense strand:", self.sense_strand, "\nGene's anti-sense strand:", self.anti_sense_strand)
        self.result.sense_strand = self.sense_strand

        if not self.amino_acid_sequence:
            # amino acid couldn't be extracted
            e = "Exception in initiate_crispr: no amino acid sequence could be extracted for your gene. something went worng connected ensembel"
            return e

        if self.amino_acid_site > len(self.amino_acid_sequence):
            e = "Exception in initiate_crispr: The amino acid site you requested: " + str(self.amino_acid_site) + \
                " does not exist, the protein sequence is " + str(len(self.amino_acid_sequence)) + " amino acid long."
            return e

        print("Gene's amino acid sequence: ", self.amino_acid_sequence)
        if self.amino_acid_sequence[self.amino_acid_site - 1] != self.from_aa.value:
            print("Exception in initiate_crispr: Amino acid in site ")
            e = "Amino acid in site " + str(self.amino_acid_site) + " is not " + \
                str(self.from_aa.name) + " but " + \
                str(AminoAcid(self.amino_acid_sequence[self.amino_acid_site - 1]).name) + "."
            return e


        self.sense_mutation_site = self.find_site_in_nt_seq(amino_acid_site=self.amino_acid_site,
                                                            check_sequence_consistency=check_consistency)
        if (self.sense_mutation_site==None):
            e = "In isoform ", str(self.given_isoform), "no amino acid in position " + str(self.amino_acid_site) + " was found." ,\
                "try to change to another isoform from WormBase."

            return e
            p
            rint("no amino acid in this position was found")

        self.choose_aa_codon()
        print(self.intronic_codon , self.type)
        if self.intronic_codon == True and self.type == "point":
            the_new_codon = self.for_intronic_codon_check_new_codon_if_posibble()
            print(the_new_codon)
            print(self.tag_for_intronic)
            if self.tag_for_intronic == 2:
                e = "The system cannot generate a Point Mutation CRISPR plan. The codon for the amino acid of interest is" \
                    " divided by an intron in the gene sequence. In such an instance generating point mutations is tricky " \
                    "and we recommend altering the amino acid codon by insertion."
                return e

        #self.codon = self.sense_strand[self.sense_mutation_site:self.sense_mutation_site + 3]


        print("mutation site in sense strand is:", self.sense_mutation_site, "and the codon is:", self.codon)
        #if self.type=="insertion" and self.codon.islower()==True:
         #   self.codon = self.codon.upper()
        if self.type=="insertion":
            self.anti_sense_mutation_site = self.get_mutation_site_for_anti_sense(self.sense_mutation_site)
            return None
        if self.codon_dic[self.codon] != self.from_aa.value:
            e = "Exception in initiate_crispr: The nucleotide site corresponding to " + \
                " amino acid site " + str(self.amino_acid_site) + " is " + str(self.sense_mutation_site) + \
                " with codon " + self.codon + ", which does not code for " + str(self.from_aa) + " but for " + \
                str(AminoAcid(self.codon_dic[self.codon])) + "."
            return e

        self.anti_sense_mutation_site = self.get_mutation_site_for_anti_sense(self.sense_mutation_site)
        print("mutation site in anti-sense strand:", self.anti_sense_mutation_site)
        return None

    def get_chosen_crrna(self, crrna, strand):
        print("crRNA we are looking for:", crrna, "in strand:", strand)
        if strand == 1:
            for sense_crrna in self.result.sense_crrnas:
                if sense_crrna[0] == crrna:
                    print("chosen:", sense_crrna)
                    return sense_crrna
        else:
            for anti_sense_crrna in self.result.anti_sense_crrnas:
                if anti_sense_crrna[0] == crrna:
                    print("chosen:", anti_sense_crrna)
                    return anti_sense_crrna
        return None

    def find_crrnas(self, window_size, PAM_size):
        print("find_crrnas sense")
        print("window_size, PAM_size",window_size, PAM_size)
        print("len of intron is - " , self.length_of_intron)
        print(self.sense_strand,self.sense_mutation_site)
        print(self.anti_sense_strand, self.anti_sense_mutation_site)
        print("self.tag_for_intronic" , self.tag_for_intronic)
        if self.tag_for_intronic == 3:
            self.sense_mutation_site = self.sense_mutation_site + self.length_of_intron + 1
            print(self.sense_strand, self.sense_mutation_site , window_size, PAM_size)
            sense_crrnas = self.get_list_of_potential_sequences(self.sense_strand,
                                                                self.sense_mutation_site,
                                                                window_size,
                                                                PAM_size)
        else:
            sense_crrnas = self.get_list_of_potential_sequences(self.sense_strand,
                                                                self.sense_mutation_site,
                                                                window_size,
                                                                PAM_size)
        print("find_crrnas anti sense")
        if self.tag_for_intronic == 3:
            self.anti_sense_mutation_site = self.anti_sense_mutation_site - self.length_of_intron
            print(self.anti_sense_strand,self.anti_sense_mutation_site,window_size,PAM_size)
            print(self.anti_sense_strand, self.anti_sense_mutation_site, window_size, PAM_size)
            anti_sense_crrnas = self.get_list_of_potential_sequences(self.anti_sense_strand,
                                                                 self.anti_sense_mutation_site,
                                                                 window_size,
                                                                 PAM_size)
        else:
            print(self.anti_sense_strand,self.anti_sense_mutation_site,window_size,PAM_size)
            anti_sense_crrnas = self.get_list_of_potential_sequences(self.anti_sense_strand,
                                                                 self.anti_sense_mutation_site,
                                                                 window_size,
                                                                 PAM_size)

        print("this is self right now:")
        print(self.gene_name
,self.type
,self.gene_id
,self.given_isoform
,self.from_aa
,self.to_aa
,self.help_isoform_id
,self.amino_acid_sequence
,self.amino_acid_site
,self.codon
,self.sense_mutation_site
,self.anti_sense_mutation_site
,self.ssODN_mutation_codon_start
,self.mutated_strand
,self.mutation_direction
,self.dsb_merged
,self.ssODN_vs_crRNA_strand
,self.DSB
,self.pam_sites
,self.mutation_zone
,self.ssODN_direction
,self.mutated_sites
,self.codon_mutations
,self.reattachment_mutations
,self.intronic_codon
,self.type_of_intronic_codon
,self.tag_for_intronic
,self.length_of_intron
   )
        return sense_crrnas, anti_sense_crrnas

    # start_crrna_index = for testing purposes
    def get_crrna(self, window_size, PAM_size, start_crrna_index=None):
        # get optional sequences for crRNA
        sense_crrnas, anti_sense_crrnas = self.find_crrnas(window_size, PAM_size)
        print("Sense crRNAs: ")
        for sequence in sense_crrnas:
            print(sequence)
            pam_site_start = sequence[1][1] + 1
            pam_site_end = sequence[1][1] + 3
            print("PAM sequence: " + self.sense_strand[pam_site_start:pam_site_end + 1])
        print("Anti Sense crRNAs: ")
        for sequence in anti_sense_crrnas:
            print(sequence)
            pam_site_start = sequence[1][1] + 1
            pam_site_end = sequence[1][1] + 3
            print("PAM sequence: " + self.anti_sense_strand[pam_site_start:pam_site_end + 1])

        return self.choose_best_crrna(sense_crrnas, anti_sense_crrnas, start_crrna_index)

    def complete_fields(self, crrna_strand_direction, to_aa, chosen_crrna):
        self.pam_sites = SequenceSites(start=chosen_crrna[1][1] + 1, end=chosen_crrna[1][1] + 3)

        self.ssODN_direction, self.mutation_direction, self.dsb_merged = \
            self.choose_ssODN_strand(crrna_strand_direction, to_aa)
        self.result.ssODN_strand = self.ssODN_direction
        print("Mutation direction is:", str(self.mutation_direction), "thus ssODN direction is:",
              str(self.ssODN_direction))
        if self.ssODN_direction != crrna_strand_direction:
            self.ssODN_vs_crRNA_strand = -1 * self.ssODN_vs_crRNA_strand
            print("The ssODN is on the opposite strand!")
            # ssODN is not on the strand from which we got the crRNA, thus we need to find its pam sites again
            self.pam_sites = SequenceSites(self.find_index_in_parallel_strand(self.pam_sites.end),
                                           self.find_index_in_parallel_strand(self.pam_sites.start))
            self.result.pam_sites = self.pam_sites
            self.result.crRNA_strand = crrna_strand_direction
            print("pam sites on the ssODN strand:", self.pam_sites)
        self.DSB = self.get_dsb(self.pam_sites, self.mutation_direction)
        self.result.pam_sites = self.pam_sites

        # to be mutated
        self.mutated_strand = self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand
        self.ssODN_mutation_codon_start = self.sense_mutation_site if self.ssODN_direction > 0 else \
            self.anti_sense_mutation_site

        # getting the zone where all mutations will be located: after to before DSB
        # sanity check
        print("mutation codon is at", self.ssODN_mutation_codon_start)

    def create_complementary_and_reverse_string(self,string):
        new=[]
        for i in range(len(string)-1,-1,-1):
            if string[i]=='A':
                new.append('T')
            if string[i]=='T':
                new.append('A')
            if string[i]=='C':
                new.append('G')
            if string[i]=='G':
                new.append('C')
            if string[i]=='a':
                new.append('t')
            if string[i]=='t':
                new.append('a')
            if string[i]=='c':
                new.append('g')
            if string[i]=='g':
                new.append('c')
        new_str="".join(new)
        print(new_str)
        return(new_str)


    #checking if we are in one of the cases where we cant destroy the PAM site (we have 4 for sense and 1 for anti sense)
    def check_speical_states_to_destroy_pam(self):
        num=0
        string=""
        print("check speical states")
        if (self.ssODN_direction == 1) and self.tag_for_intronic == 3 and self.type_of_intronic_codon ==1:
            mut_site = self.sense_mutation_site - 1
        elif (self.ssODN_direction == 1) and self.tag_for_intronic == 3 and self.type_of_intronic_codon ==2:
            mut_site = self.sense_mutation_site - 2
        else:
            mut_site = self.sense_mutation_site
        print(mut_site)

        #need to check if it work both strands
        if (self.ssODN_direction == 1):
            if(mut_site<=self.pam_sites.start):                       #case 1
                print("case1")
                num=(self.pam_sites.start-1)-(mut_site+2)
                num=num%3
                print("num= ",num)
                if (num==0):
                    print(self.pam_sites.start,self.pam_sites.end)
                    string=str(self.sense_strand[self.pam_sites.start:self.pam_sites.end+1])
                    print("string" , string)
                    self.maybe_problematic=(0,string)
                if (num==1):
                    string=str(self.sense_strand[self.pam_sites.start-1:self.pam_sites.end+3])
                    self.maybe_problematic = (2, string)
                if (num==2):
                    string = str(self.sense_strand[self.pam_sites.start - 2:self.pam_sites.end + 2])
                    self.maybe_problematic = (1, string)
            else:                                                                     #case 2
                print("case2")
                num = (mut_site - 1) - (self.pam_sites.start + 2)
                num = num % 3
                if (num == 0):
                    print(self.pam_sites.start, self.pam_sites.end)
                    string = str(self.sense_strand[self.pam_sites.start:self.pam_sites.end + 1])
                    self.maybe_problematic = (0, string)
                if (num == 1):
                    string = str(self.sense_strand[self.pam_sites.start - 1:self.pam_sites.end + 3])
                    self.maybe_problematic = (2, string)
                if (num == 2):
                    string = str(self.sense_strand[self.pam_sites.start - 2:self.pam_sites.end + 2])
                    self.maybe_problematic = (1, string)


        else:                                                               #self.ssODN_direction==-1
            if (self.anti_sense_mutation_site>=self.pam_sites.start):                                     #case 3
                print("case3")
                num=(self.anti_sense_mutation_site)-(self.pam_sites.end)-1
                num=num%3
                print("num= ",num)
                if (num==0):
                    #print(self.pam_sites.start,self.pam_sites.end)
                    string=str(self.anti_sense_strand[self.pam_sites.start:self.pam_sites.end+1])
                    string=self.create_complementary_and_reverse_string(string)
                    print("string" , string)
                    self.maybe_problematic=(0,string)
                if (num==1):
                    string=str(self.anti_sense_strand[self.pam_sites.start-2:self.pam_sites.end+2])
                    print(string)
                    string = self.create_complementary_and_reverse_string(string)
                    print(string)
                    self.maybe_problematic = (1, string)
                if (num==2):
                    string = str(self.anti_sense_strand[self.pam_sites.start - 1:self.pam_sites.end + 3])
                    string = self.create_complementary_and_reverse_string(string)
                    self.maybe_problematic = (2, string)


            if (self.anti_sense_mutation_site<self.pam_sites.start):                             #case 4
                print("case4")
                num=(self.pam_sites.start)-(self.anti_sense_mutation_site+2)-1
                num=num%3
                print("num= ",num)
                if (num==0):
                    #print(self.pam_sites.start,self.pam_sites.end)
                    string=str(self.anti_sense_strand[self.pam_sites.start:self.pam_sites.end+1])
                    string=self.create_complementary_and_reverse_string(string)
                    print("string" , string)
                    self.maybe_problematic=(0,string)
                if (num==1):
                    string=str(self.anti_sense_strand[self.pam_sites.start-1:self.pam_sites.end+3])
                    print(string)
                    string = self.create_complementary_and_reverse_string(string)
                    print(string)
                    self.maybe_problematic = (1, string)
                if (num==2):
                    string = str(self.anti_sense_strand[self.pam_sites.start - 2:self.pam_sites.end + 2])
                    string = self.create_complementary_and_reverse_string(string)
                    self.maybe_problematic = (2, string)
        return ()


    # this function is responsible for inserting all needed mutations to the string: to change amino acid, to prevent
    # reattachment and to add\remove restriction site
    def insert_mutations(self, ssODN_codon_mutation_site, restriction_site_type):
        # 1. mutation to change codon to change the amino acid
        possible_codon_mutations = self.get_possible_mutations_demands(self.from_aa,
                                                                       self.to_aa, self.sense_mutation_site)
        second_strand = self.mutated_strand
        second_mutated_sites = self.mutated_sites[:]

        # print("original mutated sites:", second_mutated_sites)
        for possible_codon_mutation in possible_codon_mutations:
            print("Now working on codon mutation:", possible_codon_mutation)
            self.mutated_strand = second_strand
            self.mutated_sites = second_mutated_sites[:]
            print("apply_codon_mutation" , ssODN_codon_mutation_site , possible_codon_mutation.dict_of_mutations)
            self.apply_codon_mutation(ssODN_codon_mutation_site,
                                      possible_codon_mutation.dict_of_mutations)
            self.codon_mutations = self.mutated_sites[:]

            print("check 14:14")
            print(self.mutated_strand)
            print(self.mutated_sites)
            print(self.codon_mutations)

            # making sure that the strand is good so far
            print("mutated strand after adding codon mutation:", self.mutated_strand)
            print("sites changed:", self.mutated_sites)

            # 2. mutations to change nt to prevent re-attachments - in PAM or crRNA sequence
            # if ssODN vs crRNA is -1, the section to be changed is PAM site, if 1 - crRNA.

            if self.ssODN_vs_crRNA_strand < 0:
                print("need to insert mutation in the pam site!")
                section_to_mutate = MutateSection(1,
                                                  DNASection.PAM_SITE,
                                                  self.pam_sites)
            else:
                print("need to insert mutation in the cRNA!")
                section_to_mutate = MutateSection(4,
                                                  DNASection.CR_RNA,
                                                  SequenceSites(self.pam_sites.start - 20,
                                                                self.pam_sites.start - 1 - 3))

            print("section to mutate to prevent re-attachment:", section_to_mutate.section_sites, "(close segment) in",
                  section_to_mutate.section_type)

            print("pam start= " ,self.pam_sites.start , "pam end= " , self.pam_sites.end)
            print("dsb =", self.DSB)
            #print("pam site= " , self.anti_sense_strand[self.pam_sites.start:self.pam_sites.end+1])
            if section_to_mutate.section_type == DNASection.PAM_SITE:
                self.check_speical_states_to_destroy_pam()
                print("the pam area and the case classification:" , self.maybe_problematic)
                lst=["CCT","CCC","CCG","CCA","TGG"]
                if (self.maybe_problematic[0]==0) and (self.maybe_problematic[1] in lst):
                    print("we are in case 1 or 4" , self.maybe_problematic)
                    self.problematic=True
                    return False
                if (self.maybe_problematic[0]==1) and (self.maybe_problematic[1][3:5]=="GG"):
                    print("we are in case 2" , self.maybe_problematic)
                    self.problematic = True
                    return False
                if (self.maybe_problematic[0]==2) and (self.maybe_problematic[1][0:4]=="ATGG" or self.maybe_problematic[1][0:4]=="TGGG"):
                    print("we are in case 3" , self.maybe_problematic)
                    self.problematic = True
                    return False
            possible_mutations, number_of_mutants = self.get_possible_codon_mutations(section_to_mutate)
            print(possible_mutations)
            if number_of_mutants == 0:  # PAM site has been mutated
                print("PAM site has been already mutated!")
                self.reattachment_mutations = self.get_mutations_in_pam()
                # 3. mutations to add/remove restriction sites
                print("Add or Remove restriction sites:")
                result = self.add_remove_restriction_sites(restriction_site_type)
                if result:
                    return True
            else:
                # subsets with indexes that represent mutations that sum up to enough mutations overall
                valid_index_subsets = CrisprPlanner.get_valid_subsets(number_of_mutants, possible_mutations, self.pam_sites,
                                                                      self.mutation_direction)
                third_strand = self.mutated_strand
                third_mutated_sites = self.mutated_sites[:]
                print("all valid index subsets to mutate reattachment section:", valid_index_subsets)
                for index_subset in valid_index_subsets:
                    self.mutated_strand = third_strand
                    self.mutated_sites = third_mutated_sites[:]
                    point_mutations = self.get_point_mutations(possible_mutations,
                                                               index_subset)
                    self.mutated_strand = self.change_chars_in_string(self.mutated_strand, point_mutations)
                    print("For reattachment prevention now working on index subset:", index_subset, "and point mutations:",
                          point_mutations)
                    print("mutated strand after adding anti-reattachment mutations:", self.mutated_strand)
                    self.mutated_sites.extend(point_mutations)
                    self.reattachment_mutations = point_mutations
                    self.change_chars_in_string(self.mutated_strand, self.reattachment_mutations)
                    print("all mutated sites:", self.mutated_sites)

                    # 3. mutations to add/remove restriction sites
                    print("Add or Remove restriction sites:")
                    result = self.add_remove_restriction_sites(restriction_site_type)
                    if result:
                        return True
        return False

    # receives (1) number of mutations needed, (2) the list of all possible anti-reattachment section mutations, and
    # (3) power set of all mutations indexes, and leave out only the ones with sufficient amount of mutations. if there
    # are not any, leaves only subsets with maximum number of mutations


    @staticmethod
    def get_valid_subsets(number_of_mutations, possible_mutations, pam_sites, mutation_direction):
        s = list(range(len(possible_mutations)))
        x = len(s)
        valid_subsets = []
        underscored_subsets = []
        max_underscored = 0
        for i in range(1 << x):
            subset = [s[j] for j in range(x) if (i & (1 << j))]
            #print("subset" , subset)
            mutations_in_subset = 0
            sum_usage = 0
            for index in subset:
                if mutations_in_subset >= number_of_mutations:
                    break
                mutations_in_subset += possible_mutations[index][1].number_of_mutations
                sum_usage += possible_mutations[index][1].usage
            else:
                if mutations_in_subset < number_of_mutations:  # not enough mutations
                    if mutations_in_subset >= max_underscored:
                        if mutations_in_subset > max_underscored:
                            underscored_subsets = []
                            max_underscored = mutations_in_subset
                        if not CrisprPlanner.are_codons_overriden(subset, possible_mutations):
                            underscored_subsets.append((subset, mutations_in_subset, sum_usage))
                else:
                    # check if mutations don't override each other by checking codon sites
                    if not CrisprPlanner.are_codons_overriden(subset, possible_mutations):
                        valid_subsets.append((subset, mutations_in_subset, sum_usage))
        print("valid_subsets ",valid_subsets)
        if valid_subsets != [([], 0, 0)]:
            # sort subsets according to number of mutations max usage
            print("valid subsets:", valid_subsets)
            # sort according to minimal distance from the DSB, minimal number of mutations and maximal usage percentage
            valid_subsets.sort(key=lambda mutation: (CrisprPlanner.get_group_distance_from_dsb(mutation[0],
                                                                                               possible_mutations,
                                                                                               pam_sites,
                                                                                               mutation_direction),
                                                     mutation[1],
                                                     -mutation[2]))
            print("valid subsets after sort:", valid_subsets)
            return list(map(lambda item: item[0], valid_subsets))
        else:
            print("WARNING: the maximum amount of mutations is", max_underscored, ", there are no subsets with",
                  number_of_mutations, "mutations")
            # sort according to number of mutations first, and codons usage second
            underscored_subsets.sort(key=lambda mutation: (mutation[1], -mutation[2]))
            return list(map(lambda item: item[0], underscored_subsets))




    """
    @staticmethod
    def get_valid_subsets(number_of_mutations, possible_mutations, pam_sites, mutation_direction):
        index_power_set = CrisprPlanner.get_power_set(lst=list(range(len(possible_mutations))),
                                                      max_num=len(possible_mutations))
        valid_subsets = []
        underscored_subsets = []
        max_underscored = 0
        for subset in index_power_set:
            print("the subset ",subset)
            mutations_in_subset = 0
            sum_usage = 0
            for index in subset:
                if mutations_in_subset >= number_of_mutations:
                    break
                mutations_in_subset += possible_mutations[index][1].number_of_mutations
                sum_usage += possible_mutations[index][1].usage
            else:
                if mutations_in_subset < number_of_mutations:  # not enough mutations
                    if mutations_in_subset >= max_underscored:
                        if mutations_in_subset > max_underscored:
                            underscored_subsets = []
                            max_underscored = mutations_in_subset
                        if not CrisprPlanner.are_codons_overriden(subset, possible_mutations):
                            underscored_subsets.append((subset, mutations_in_subset, sum_usage))
                else:
                    # check if mutations don't override each other by checking codon sites
                    if not CrisprPlanner.are_codons_overriden(subset, possible_mutations):
                        valid_subsets.append((subset, mutations_in_subset, sum_usage))
        print("valid_subsets ",valid_subsets)
        if valid_subsets != [([], 0, 0)]:
            # sort subsets according to number of mutations max usage
            print("valid subsets:", valid_subsets)
            # sort according to minimal distance from the DSB, minimal number of mutations and maximal usage percentage
            valid_subsets.sort(key=lambda mutation: (CrisprPlanner.get_group_distance_from_dsb(mutation[0],
                                                                                               possible_mutations,
                                                                                               pam_sites,
                                                                                               mutation_direction),
                                                     mutation[1],
                                                     -mutation[2]))
            print("valid subsets after sort:", valid_subsets)
            return list(map(lambda item: item[0], valid_subsets))
        else:
            print("WARNING: the maximum amount of mutations is", max_underscored, ", there are no subsets with",
                  number_of_mutations, "mutations")
            # sort according to number of mutations first, and codons usage second
            underscored_subsets.sort(key=lambda mutation: (mutation[1], -mutation[2]))
            return list(map(lambda item: item[0], underscored_subsets))

"""

    # returns the distance of the furthest point from the DSB
    @staticmethod
    def get_group_distance_from_dsb(subset, possible_mutations, pam_sites, mutation_direction):
        DSB = CrisprPlanner.get_dsb(pam_sites, mutation_direction)
        max_distance = 0
        for index in subset:
            start = possible_mutations[index][0].codon_sites.start
            average_mutated_index = mean([start + item for item in possible_mutations[index][1].dict_of_mutations])
            if abs(average_mutated_index - DSB) > max_distance:
                max_distance = abs(average_mutated_index - DSB)
        return max_distance

    @staticmethod
    def get_dsb(pam_sites, mutation_direction):
        if mutation_direction == MutationDirection.DOWNSTREAM:
            DSB = pam_sites.end + 4
        else:
            DSB = pam_sites.start - 3
        return DSB

    # checks if the subset points to mutations that try to mutate the same codons
    @staticmethod
    def are_codons_overriden(subset, possible_mutations):
        codons = set()
        for index in subset:
            mutation = possible_mutations[index]
            codon_data = mutation[0]
            if codon_data in codons:
                return True
            codons.add(codon_data)
        return False

    # receives (1) the section that needs to be mutated to prevent reattachment, (2) the indexes of the amino acid's
    # mutation codon on the ssODN strand,(3) the mutation direction, and (4) the sites that have
    #  changed while mutating the amino acid, computes which indexes can be changed in order to insert the section with
    # silent mutations and returns the possibilities for silent mutations
    def get_possible_codon_mutations(self, section_to_mutate: MutateSection):
        # format of codons_data: ('GCG', (1, 3))
        if (self.ssODN_direction == 1) and self.tag_for_intronic == 3 and self.type_of_intronic_codon ==1:
            mut_site = self.sense_mutation_site - 1
        elif (self.ssODN_direction == 1) and self.tag_for_intronic == 3 and self.type_of_intronic_codon ==2:
            mut_site = self.sense_mutation_site - 2
        else:
            mut_site = self.sense_mutation_site
        if (self.ssODN_direction == -1):
            mut_site = self.ssODN_mutation_codon_start

        codons_data = CrisprPlanner.get_relevant_codons(section_to_mutate, mut_site,
                                                        self.mutated_strand)
        possible_mutations = []
        for codon_data in codons_data:
            # all codons_data that translate to the same amino acid as the relevant codon does
            same_aa_codons = self.get_similar_codons(codon_data.codon)
            if same_aa_codons is None:
                continue  # intron
            print("for codon:", codon_data, "same aa codons:", same_aa_codons)
            if section_to_mutate.section_type == DNASection.PAM_SITE:
                # removes PAM mutations of sequence NGA and NGG
                CrisprPlanner.check_PAM_mutation(codon_data, same_aa_codons, section_to_mutate)
            # removes codons that run over the nucleotide change to change amino acid
            CrisprPlanner.check_already_mutated_sites(codon_data, same_aa_codons, self.mutated_sites)
            # removes codons in which the difference is in part that in outside of the mutation zone
            CrisprPlanner.check_outside_codon_mutations(codon_data, same_aa_codons, section_to_mutate)

            optional_codon_mutations = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(codon_data.codon,
                                                                                            same_aa_codons)
            for optional_mutation in optional_codon_mutations:
                possible_mutations.append((codon_data, optional_mutation))
        # sort mutations according to which is closest the the mean of former point mutations
        possible_mutations.sort(
            key=lambda x: abs(x[0].codon_sites.start + mean(x[1].dict_of_mutations.keys()) - self.DSB))
        print("possible mutations:", *possible_mutations, sep="\n")
        number_of_mutants = CrisprPlanner.get_number_of_mutants(section_to_mutate, self.mutated_sites,
                                                                self.mutation_direction)
        return possible_mutations, number_of_mutants

    # receives (1) codon data (2) list of codons that are translated into same amino acid as the codon, and (3) the
    # section to mutate, and removes all codons from the same_aa_codons_ list in which the difference in nucleotides
    # is outside of the section to mutate. Also if the section is crRNA we leave out all codons that write on the DSB.
    @staticmethod
    def check_outside_codon_mutations(codon_data, same_aa_codons, section_to_mutate: MutateSection):
        codon_start = codon_data.codon_sites.start
        for same_aa_codon in same_aa_codons[:]:
            _, demands = CrisprPlanner.codon_distance(codon_data.codon, same_aa_codon)
            print(same_aa_codon, ":", demands)
            valid_demands = False
            dsb_safe = True
            for index in demands:
                if section_to_mutate.section_sites.start <= codon_start + index <= section_to_mutate.section_sites.end:
                    valid_demands = True
                if section_to_mutate.section_type == DNASection.CR_RNA:
                    if section_to_mutate.section_sites.end < codon_start + index:
                        dsb_safe = False
            if not valid_demands:
                print("codon:", same_aa_codon, "is removed because of outside mutations")
                same_aa_codons.remove(same_aa_codon)
            if same_aa_codon in same_aa_codons and not dsb_safe:
                print("codon:", same_aa_codon, "is removed because of writing over DSB")
                same_aa_codons.remove(same_aa_codon)

    # receives (1) all possible mutations to prevent reattachment and (2) list of indexes corresponding to the mutations
    # list and returns the point mutations in format of namedTuples: PointMutations(index, new nucleotide)
    def get_point_mutations(self, possible_mutations, indexes_subset):
        # possible mutation format: (CodonData(codon='TTT', codon_sites=SequenceSites(start=1174, end=1176)),
        # OptionalCodonChange('TTC', 1, {2: 'C'}, 23.9))
        point_mutations = []
        for index in indexes_subset:
            possible_mutation = possible_mutations[index]
            dict_of_changes = possible_mutation[1].dict_of_mutations
            for key in dict_of_changes:
                point_mutations.append(
                    PointMutation(possible_mutation[0].codon_sites.start + key,
                                  possible_mutation[0].codon[key],
                                  dict_of_changes[key]))

        return point_mutations

    @staticmethod
    def get_power_set(lst, max_num, power_set=[[]]):
        if lst:
            new_power_set = []
            new_item = lst[0]
            for item in power_set:
                if len(item) < max_num:
                    new_power_set.append(item + [new_item])
                new_power_set.append(item)
            return CrisprPlanner.get_power_set(lst[1:], max_num, new_power_set)
        else:
            print("power_set")
            print(power_set)
            return power_set

    # Takes into consideration whether mutations to change amino acid have also modified the reattachment sequence, and
    # thus we need fewer mutations to introduce to the reattachment sequence
    @staticmethod
    def get_number_of_mutants(section_to_mutate, mutated_sites, mutation_direction):
        number_of_mutants = section_to_mutate.number_of_mutations
        plus = 0 if section_to_mutate.section_type == DNASection.PAM_SITE else 2
        for mutated_site in mutated_sites:
            if section_to_mutate.section_sites.start <= mutated_site.index <= section_to_mutate.section_sites[1] + plus:
                if section_to_mutate.section_type == DNASection.PAM_SITE:
                    if CrisprPlanner.disregard_mutation_in_pam(section_to_mutate, mutated_site, mutation_direction):
                        continue
                number_of_mutants -= 1
                print("already a mutant in nucleotide:", mutated_site)
        if number_of_mutants < 0:
            number_of_mutants = 0
        return number_of_mutants

    # picks the point mutations in the codon point mutations that mutate the PAM site
    def get_mutations_in_pam(self):
        pam_mutations = []
        for codon_mutation in self.codon_mutations:
            if self.pam_sites.start <= codon_mutation.index <= self.pam_sites.end:
                pam_mutations.append(codon_mutation)
        print("pam mutations:", pam_mutations)
        return pam_mutations

    @staticmethod
    def disregard_mutation_in_pam(section_to_mutate, mutated_site, mutation_direction):
        if mutation_direction == MutationDirection.UPSTREAM and \
                mutated_site.index == section_to_mutate.section_sites.start:  # N to N
            return True
        if mutation_direction == MutationDirection.UPSTREAM and \
                mutated_site.index == section_to_mutate.section_sites.end and mutated_site.new_nucleotide == 'A':  # NGA
            return True
        if mutation_direction == MutationDirection.DOWNSTREAM and \
                mutated_site.index == section_to_mutate.section_sites.end:  # N to N
            return True
        if mutation_direction == MutationDirection.DOWNSTREAM and \
                mutated_site.index == section_to_mutate.section_sites.start and mutated_site.new_nucleotide == 'T':  # TCN
            return True
        return False

    # receives (1) codon_data the contains the codon and its sites, (2) possible codons that the former codons could be
    # mutated into but keep the amino acid, and (3) sites that have changed and the new mutations shouldn't run over.
    # The function goes through each of the mutated site (input 3) and checks if it is inside the range of the new
    # codon sites. If so, it checks whether in the possible codons list there are codons that run over that site (change
    #  nucleotide), and if there are, removes them from the possible codons list.
    @staticmethod
    def check_already_mutated_sites(codon_data, possible_codons, mutated_sites):
        print("checking if some codons need to be removed to prevent them from changing mutations")
        for mutated_site in mutated_sites:
            if codon_data.codon_sites.start <= mutated_site.index <= codon_data.codon_sites.end:
                for possible_codon in possible_codons[:]:
                    # mutated_codon_details format = (1, {2: 'G'}, 14.3)
                    mutated_codon_details = CrisprPlanner.how_to_get_b_codons_from_a(codon_data.codon, [possible_codon])
                    for key in mutated_codon_details.dict_of_mutations:
                        # for each nucleotide mutation, key is the nucleotide place in the codon: 0,1,2
                        if codon_data.codon_sites[0] + key == mutated_site.index:
                            print(mutated_codon_details, codon_data.codon_sites, "removed because", mutated_site.index,
                                  "has already been mutated")
                            possible_codons.remove(possible_codon)
                            break

    # receives (1) list of mutations to check, (2) list of possible codons in that site, (3) section to mutate (PAM site
    #  or crRNA and the details), and (4) mutation direction (upstream or downstream) and checks if one of the PAM
    # mutations is substituting NGG site with an NGA or NGG site, and if so, removes it from the valid list of codons
    # returned.
    @staticmethod
    def check_PAM_mutation(codon_data: CodonData, possible_codons, section_to_mutate):
        if codon_data.codon_sites.start < section_to_mutate.section_sites.start:
            print("case 1 - need to check!!!!!!")
            # before beginning of PAM site
            if codon_data.codon_sites.end == section_to_mutate.section_sites.start:
                # only contains first nucleotide
                for possible_codon in possible_codons[:]:
                    if possible_codon[2] == 'T':
                        print("possible codon removed:", possible_codon, "in", codon_data)
                        possible_codons.remove(possible_codon)
            elif codon_data.codon_sites.end > section_to_mutate.section_sites.start:
                # contains first two amino acids
                for possible_codon in possible_codons[:]:
                    if possible_codon[1:] == 'TC':
                        print("possible codon removed:", possible_codon, "in", codon_data)
                        possible_codons.remove(possible_codon)
        elif codon_data.codon_sites == section_to_mutate.section_sites:
            print("case 2 - checked")
            for possible_codon in possible_codons[:]:
                # the site is CCN
                print(possible_codon)
                print(possible_codon[:2])
                if (possible_codon[:2] == 'TC' or  possible_codon[:2] =='CC'):
                    print("possible codon removed:", possible_codon, "in", codon_data)
                    possible_codons.remove(possible_codon)
        else:
            print("case 3 - need to check!!!!!!")
            if codon_data.codon_sites.start < section_to_mutate.section_sites.end:
                # contains two last amino acids
                for possible_codon in possible_codons[:]:
                    if possible_codon[0] == 'C':
                        print("possible codon removed:", possible_codon, "in", codon_data)
                        possible_codons.remove(possible_codon)
            elif codon_data.codon_sites.start == section_to_mutate.section_sites.end:
                # contains only last PAM nucleotide
                print(
                    "all mutations of this sort must be erased since mutation in that codon will leave out the NGG")
                possible_codons.clear()

    @staticmethod
    def is_codon_in_intron(codon):
        for nuc in codon:
            if 'a' <= nuc <= 'z':  # small letter, intron
                return True
        return False

    # receives (1) a codon and returns all other codons that translated to the same amino acid
    def get_similar_codons(self, codon):
        if self.is_codon_in_intron(codon):
            return None  # intron
        if self.ssODN_direction > 0:
            amino_acid = CrisprPlanner.codon_dic[codon]
            relevant_codons = self.amino_acid_dic[amino_acid].copy()
            relevant_codons.remove(codon)
            return list(relevant_codons)
        else:  # ssODN is on anti-sense, PAM is changed
            complementary_codon = self.get_complementary_sequence(codon)
            amino_acid = CrisprPlanner.codon_dic[complementary_codon]
            relevant_codons = self.amino_acid_dic[amino_acid].copy()
            relevant_codons.remove(complementary_codon)
            complementary_similar_codons = set()
            for relevant_codon in relevant_codons:
                complementary_similar_codons.add(self.get_complementary_sequence(relevant_codon))
            return list(complementary_similar_codons)

    # receives (1) the section to be mutated (PAM site or crRNA indexes), (2) the aa mutation codon indexes of the ssODN
    # strand, and (3) the mutated strand (sense or anti-sense) and returns a list of all codons and their indexes that
    # are part of the section to be mutated indexes
    # format of the results: ('GCG', (1, 3))
    @staticmethod
    def get_relevant_codons(section_to_mutate: MutateSection, ssODN_mutation_codon_start, mutated_strand):
        print("func: get_relevant_codons")
        print( section_to_mutate, ssODN_mutation_codon_start, mutated_strand)
        codons_data = []
        remainder = abs(section_to_mutate.section_sites.start - ssODN_mutation_codon_start) % 3
        print("section_to_mutate.section_sites.start" , section_to_mutate.section_sites.start)
        print("ssODN_mutation_codon_start" , ssODN_mutation_codon_start)

        print(remainder)
        if section_to_mutate.section_sites.start >= ssODN_mutation_codon_start:
            start_point = section_to_mutate.section_sites.start - remainder
        else:
            if remainder > 0:
                start_point = section_to_mutate.section_sites.start - (3 - remainder)
            else:
                start_point = section_to_mutate.section_sites.start
        if section_to_mutate.section_type == DNASection.PAM_SITE:
            end_point = section_to_mutate.section_sites.end  # codons can cross the end of section
        else:
            # if section is crRNA or mutation zone, codons shouldn't cross the DSB
            end_point = section_to_mutate.section_sites.end - 2
        while start_point <= end_point:  # far enough so a codon won't cross the end:
            codon = mutated_strand[start_point:start_point + 3]
            codons_data.append(CodonData(codon, SequenceSites(start_point, start_point + 2)))
            start_point += 3
        print(str(len(codons_data)) + " relevant codons that we can mutate:", *codons_data, sep="\n")
        return codons_data

    # This function receives (1) an amino acid sequence site and (2) value that indicates whether to check the
    # correlation between amino acid sequence and nucleotide sequence (codon to amino acid) or not, and returns the
    # site of the mutation in the nucleotide sequence, meaning the site or the corresponding codon
    def find_site_in_nt_seq(self, amino_acid_site, check_sequence_consistency: bool = True):
        print("b"," ",amino_acid_site," ", check_sequence_consistency)
        no_utr_strand = self.remove_utr()
        print("sense without utr:", no_utr_strand)
        if check_sequence_consistency:
            # consistency check
            nt_index = 0
            aa_index = 0
            while aa_index < len(self.amino_acid_sequence):
                CrisprPlanner.check_sequence_consistency(no_utr_strand[nt_index:nt_index + 3],
                                                         self.amino_acid_sequence[aa_index])
                nt_index += 3
                aa_index += 1
        nt_site_in_sense_without_utr = (amino_acid_site - 1) * 3
        print("nt site " , nt_site_in_sense_without_utr)
        utrs = 0
        seq_index = 0
        exon_index = 0
        print("len of seq ", len(self.sense_strand))
        while seq_index < len(self.sense_strand):
            #print(self.sense_strand[seq_index])
            if self.sense_strand[seq_index].islower():
                # region is intron or utr
                utrs += 1
            else:
                # region is exon
                if nt_site_in_sense_without_utr == 0:
                    return exon_index + utrs
                exon_index += 1
                #print(exon_index)
                if exon_index == nt_site_in_sense_without_utr:
                    print("hereeeeeeeeeeee ",exon_index + utrs)
                    return exon_index + utrs
            seq_index += 1


    # receives (1) a codon string and (2) the relevant amino and checks if the codon translates to this amino acid
    @staticmethod
    def check_sequence_consistency(input_codon, amino_acid):
        codon_table = {'GCT': 'A', 'TGT': 'C', 'GAT': 'D', 'GAA': 'E', 'TTT': 'F', 'GGT': 'G', 'CAT': 'H',
                       'ATT': 'I', 'AAA': 'K', 'TTA': 'L', 'ATG': 'M', 'AAT': 'N', 'CCT': 'P', 'CAA': 'Q',
                       'CGT': 'R', 'TCT': 'S', 'ACT': 'T', 'GTT': 'V', 'TGG': 'W', 'TAT': 'Y', 'TTC': 'F',
                       'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATC': 'I', 'ATA': 'I',
                       'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCC': 'P',
                       'CCA': 'P', 'CCG': 'P', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCC': 'A', 'GCA': 'A',
                       'GCG': 'A', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', 'CAC': 'H', 'CAG': 'Q', 'AAC': 'N',
                       'AAG': 'K', 'GAC': 'D', 'GAG': 'E', 'TGC': 'C', 'TGA': 'STOP', 'CGC': 'R', 'CGA': 'R',
                       'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGC': 'G', 'GGA': 'G',
                       'GGG': 'G', '': 'END'}
        if codon_table[input_codon] != amino_acid:
            print("Instead of " + codon_table[input_codon], ", the amino acid in the sequence is " + amino_acid)

    # receives a (1) strand of nt, (2) the mutation site, (3) the window size around the mutation site where DSB point
    # can be located and (4) size of PAM sequence, checks in the strand for 20 bp sequences around the mutation size
    # that end in PAM sequence ('NGG') and returns a list of those sequences, as a tuple of the sequences and its start
    # and end indexes.
    @staticmethod
    def get_list_of_potential_sequences(strand, mutation_site, window_size, PAM_size, cr_rna_size: int = 20):
        print("get_list_of_potential_sequences")
        print(strand, mutation_site, window_size, PAM_size, cr_rna_size)
        potential_cr_rnas = []

        start_point = max(0, mutation_site - window_size + PAM_size)
        end_point = min(mutation_site + window_size + PAM_size, len(strand) - 1)

        print("range to find crRNA:", start_point, end_point)
        for i in range(start_point, end_point):
            if i + 3 < len(strand):
                if strand[i + 1:i + 3] == "gg" or strand[i + 1:i + 3] == "GG":
                    print("in crna finding!")
                    print(strand[i + 1:i + 3] , i)
                    print(cr_rna_size)
                    try:
                        if i - cr_rna_size >= 0:
                            potential_cr_rna = strand[i - cr_rna_size:i]
                            potential_cr_rnas.append((potential_cr_rna, SequenceSites(i - cr_rna_size, i - 1)))
                    except Exception as e:
                        print("Not enough space for a whole crRNA, i = " + str(i))
                        print("Exception in get_list_of_potential_sequences:", e)
                        continue
        return potential_cr_rnas

    # receives (1) the mutation position in the sense strand, and returns the position of this mutation in the
    # anti-sense
    def get_mutation_site_for_anti_sense(self, sense_mutation_site):
        return len(self.anti_sense_strand) - sense_mutation_site - 3

    # receives (1) a sequence and returns its complementary, from 5' to 3'
    @staticmethod
    def get_complementary_sequence(given_sense_strand):
        if not given_sense_strand:
            return None
        anti_sense_strand = ''
        pair_nucleotides = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
        for letter in given_sense_strand:
            anti_sense_strand += pair_nucleotides[letter]
        return anti_sense_strand[::-1]

    # receives (1) a list of sense strands and returns their complementary strands, from 5' to 3'
    @staticmethod
    def get_anti_sense_strands(given_sense_strands):
        complementary_strands = []
        for sense_strand in given_sense_strands:
            anti_sense_strand = CrisprPlanner.get_complementary_sequence(sense_strand)
            complementary_strands.append(anti_sense_strand)
        return complementary_strands

    # this function gets a nucleotide sequence and adds to a dictionary all restriction sites in the sequence
    def find_restriction_sites(self, nucleotide_seq, start_point, r_type):
        restriction_sites_in_seq = set()
        for restriction_enzyme in self.restriction_enzymes:
            for restriction_site in restriction_enzyme.derivatives:
                complementary_restriction_site = self.get_complementary_sequence(restriction_site)
                sites = [restriction_site, complementary_restriction_site]
                for site in sites:
                    if site in nucleotide_seq:
                        indexes = [m.start() for m in re.finditer('(?=' + site + ')', nucleotide_seq)]
                        for index in indexes:
                            if restriction_site == "GGTGA":
                                if index + 8 >= len(nucleotide_seq):
                                    continue
                            seq_sites = SequenceSites(start_point + index,
                                                      start_point + index + len(site) - 1)  # closed section
                            if seq_sites.end < self.mutation_zone.start or seq_sites.start > self.mutation_zone.end:
                                continue
                            restriction_sites_in_seq.add(RestrictionSite(seq_sites, restriction_enzyme, r_type))
        # print("Restriction Sites:", *restriction_sites_in_seq, sep='\n')
        return list(restriction_sites_in_seq)

    @staticmethod
    def change_chars_in_string(seq, point_mutations: list):
        lst = list(seq)
        for point_mutation in point_mutations:
            lst[point_mutation.index] = point_mutation.new_nucleotide
        return "".join(lst)

    @staticmethod
    def change_char_in_string(seq, position, new_char):
        lst = list(seq)
        lst[position] = new_char
        return "".join(lst)

    @staticmethod
    def choose_best_crrna(sense_options: list, anti_sense_options: list, start_crrna_index):
        if not start_crrna_index:
            num = random.uniform(0, 1)
            if num > 0.5:
                return random.choice(sense_options), 1
            else:
                return random.choice(anti_sense_options), -1
        else:  # TESTING
            for crrna in sense_options:
                if crrna[1].start == start_crrna_index:
                    return crrna, 1
            for crrna in anti_sense_options:
                if crrna[1].start == start_crrna_index:
                    return crrna, -1
        return None, None

    def choose_ssODN_strand(self, crrna_strand_direction, to_aa: AminoAcid):
        if self.dsb_merged:
            return crrna_strand_direction, MutationDirection.UPSTREAM, False
        codon_dsb_merged = False

        def check_inside_dsb(sense_strand, sense_mutation_site):
            inner_codon_dsb_merged = False
            downstream = 0
            upstream = 0
            if crrna_strand_direction < 0:
                # anti-sense, to_aa needs to get complementary
                comp_codon_mutation = \
                    CrisprPlanner.how_to_get_b_from_a(sense_strand[sense_mutation_site:sense_mutation_site + 3], to_aa)
                dict_of_mutations = {}
                for key in comp_codon_mutation.dict_of_mutations:
                    dict_of_mutations[2 - key] = CrisprPlanner.get_complementary_sequence(
                        comp_codon_mutation.dict_of_mutations[key])
                codon = CrisprPlanner.get_complementary_sequence(comp_codon_mutation.codon)
                codon_mutation = CodonMutation(codon=codon,
                                               number_of_mutations=comp_codon_mutation.number_of_mutations,
                                               dict_of_mutations=dict_of_mutations,
                                               usage=CrisprPlanner.codon_usage[CrisprPlanner.change_t_to_u(codon)])
            else:
                codon_mutation = CrisprPlanner.how_to_get_b_from_a(sense_strand[mutation_site:mutation_site + 3], to_aa)
            print(codon_mutation)
            for key in codon_mutation.dict_of_mutations:
                if mutation_site + key >= DSB:
                    downstream += 1
                else:
                    upstream += 1
            print("downstream:", downstream, "upstream:", upstream)
            if upstream > downstream:
                return MutationDirection.UPSTREAM, inner_codon_dsb_merged
            if upstream == downstream:
                inner_codon_dsb_merged = True
            return MutationDirection.DOWNSTREAM, inner_codon_dsb_merged

        mutation_site = self.sense_mutation_site if crrna_strand_direction > 0 else self.anti_sense_mutation_site

        DSB = self.pam_sites.start - 3
        # default assignment
        mutation_direction = MutationDirection.UPSTREAM
        print(
            "mutation site: " + str(mutation_site) + ", pam start: " + str(self.pam_sites.start) + ", DSB: " + str(DSB))
        if mutation_site <= DSB <= mutation_site + 2:
            # DSB is inside the codon we need to mutate, further investigation is needed
            mutation_direction, codon_dsb_merged = check_inside_dsb(self.sense_strand, self.sense_mutation_site)
        elif mutation_site > DSB:
            mutation_direction = MutationDirection.DOWNSTREAM
        elif mutation_site < DSB:
            mutation_direction = MutationDirection.UPSTREAM
        if mutation_direction == MutationDirection.DOWNSTREAM:
            return -1 * crrna_strand_direction, mutation_direction, codon_dsb_merged
        else:
            return crrna_strand_direction, mutation_direction, codon_dsb_merged

    def find_index_in_parallel_strand(self, index):
        return len(self.sense_strand) - index - 1

    # this function will return the minimal nucleotides you need to change (and the range of options for change) in
    # order to get from codon a to amino acid b
    @staticmethod
    def how_to_get_b_from_a(a_codon, b_amino_acid: AminoAcid):
        b_codons = CrisprPlanner.amino_acid_dic[b_amino_acid.value]
        return CrisprPlanner.how_to_get_b_codons_from_a(a_codon, b_codons)

    # this function will return the minimal nucleotides you need to change (and the range of options for change) in
    # order to get from codon a to possible codons b. result in format: ('GCG', (1, {0: 'G'}, 8.2))
    @staticmethod
    def how_to_get_b_codons_from_a(a_codon, b_codons):
        min_distance = 4
        options = {}
        for b_codon in b_codons:
            distance, demands = CrisprPlanner.codon_distance(a_codon, b_codon)
            if distance < min_distance:
                options = {}
                min_distance = distance
                options[b_codon] = (distance, demands, CrisprPlanner.codon_usage[CrisprPlanner.change_t_to_u(b_codon)])
            elif distance == min_distance:
                options[b_codon] = (distance, demands, CrisprPlanner.codon_usage[CrisprPlanner.change_t_to_u(b_codon)])
        if options:
            # result is a tuple of codon, details
            return CrisprPlanner.get_codon_with_higher_usage(options)
        return None

    # this function will return a list of all possible codon changes, sorted by the number of point mutations needed and
    # by the usage value of the codon, in order to get from codon a to possible codons b.
    # result in format of list of OptionlCodonChange
    @staticmethod
    def get_list_of_how_to_get_b_codons_from_a(a_codon, b_codons):
        options = []
        for b_codon in b_codons:
            distance, demands = CrisprPlanner.codon_distance(a_codon, b_codon)
            options.append(CodonMutation(b_codon,
                                         distance,
                                         demands,
                                         CrisprPlanner.codon_usage[CrisprPlanner.change_t_to_u(b_codon)]))
        if options:
            # result is a list of sorted OptionalCodonChange
            options.sort(key=lambda option: (option.number_of_mutations, -option.usage))
        return options

    @staticmethod
    def codon_distance(a_codon, b_codon):
        distance = 0
        demands = {}
        for i in range(3):
            if a_codon[i] != b_codon[i]:
                distance += 1
                demands[i] = b_codon[i]
        return distance, demands

    @staticmethod
    def change_t_to_u(codon):
        new_codon = ""
        for nt in codon:
            if nt == 'T':
                new_codon += 'U'
            else:
                new_codon += nt
        return new_codon

    # gets a dictionary of options, such that the keys are tuples of the original codon (doesn't change) and new, and
    # the values are 3-tuples which consists of the number of mutations, the instructions (dictionary of codon indexes and
    # wanted nucleotides, and the usage percentage, and returns a MutatedCodonDetails that holds the said info + the key
    # codon with the highest codon usage
    @staticmethod
    def get_codon_with_higher_usage(codon_options):
        max_codon = ""
        max_usage = 0
        for codon in codon_options:
            value = codon_options[codon]
            usage = value[2]
            if usage > max_usage:
                max_usage = usage
                max_codon = codon
        codon_mutation = CodonMutation(max_codon, codon_options[max_codon][0], codon_options[max_codon][1], max_usage)
        # print(codon_mutation)
        return codon_mutation

    # calculated the demands to go through the codon from the mutation site to one of the codons translating to the
    # wanted amino acid
    def get_possible_mutations_demands(self, from_amino_acid: AminoAcid, to_amino_acid: AminoAcid, mutation_site):
        #sense_codon = self.sense_strand[mutation_site:mutation_site + 3]
        sense_codon = self.codon
        codon_options = self.amino_acid_dic[from_amino_acid.value]



        if sense_codon not in codon_options:
            print(f'Codon in mutation site {mutation_site}: {sense_codon} does not translate into {from_amino_acid}')
            return []
        a_codon = sense_codon if self.ssODN_direction > 0 else CrisprPlanner.get_complementary_sequence(sense_codon)
        b_codons = CrisprPlanner.amino_acid_dic[to_amino_acid.value]
        if self.ssODN_direction > 0:
            b_codons_info = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(a_codon, list(b_codons))
        else:
            b_codons_info = CrisprPlanner.get_list_of_how_to_get_b_codons_from_a(a_codon,
                                                                                 CrisprPlanner.get_anti_sense_strands(
                                                                                     b_codons))
        print("\n" + str(len(b_codons_info)), "possible codons change:", b_codons_info)
        print("The best option is to change", a_codon, "to", b_codons_info[0].codon, "with demands:",
              str(b_codons_info[0].dict_of_mutations))

        return b_codons_info

    # receives (1) the mutation site in the strand (internal indexing), and (2) the demands to make the change with
    # indexes in respect to the codon, not the strand (so, 0, 1, 2...) and mutate the
    # mutated strand and adds to the list of mutated sites in format of PointMutation tuples: (index, new nucleotide)
    def apply_codon_mutation(self, mutation_site, demands):
        print("mutation_site" , mutation_site)
        print("demands" , demands)
        point_mutations = []
        for index in demands:
            nt = demands[index]
            if self.tag_for_intronic == 3 and self.type_of_intronic_codon == 1:
                print("apply_codon_mutation first if")
                position = mutation_site + int(index) - 1
            elif self.tag_for_intronic == 3 and self.type_of_intronic_codon == 2:
                print("apply_codon_mutation second if")
                position = mutation_site + int(index) - 2
            else:
                print("apply_codon_mutation else")
                position = mutation_site + int(index)
            point_mutations.append(PointMutation(position,
                                                 self.mutated_strand[position],
                                                 nt))
            print(nt , position , point_mutations)
        self.mutated_strand = CrisprPlanner.change_chars_in_string(self.mutated_strand, point_mutations)
        self.mutated_sites.extend(point_mutations)

    @staticmethod
    def check_if_pam_site_mutated(mutated_sites: list, pam_sites: SequenceSites):
        pam_site_mutated = False
        for mutated_site in mutated_sites:
            if pam_sites.start <= mutated_site <= pam_sites.end:
                pam_site_mutated = True
                break

        if pam_site_mutated:
            print("pam site mutated")
        else:
            print("pam site not mutated")
        return pam_site_mutated

    # receives (1) the mutation direction (UPSTREAM or DOWNSTREAM) and (2) the pam sites and according to the scheme
    # that shows that for UPSRTREAM mutations, the homology arm that can undergo changes is from the DSB and for
    # DOWNSTREAM mutations, the homology arm that can undergo changes is up until the DSB, returns the right zone of
    # 13 nt. The range is open in the end (doesn't include the end nucleotide)
    def get_mutations_zone(self):
        lowest_site = min(list(map(lambda item: item.index, self.mutated_sites)))

        DSB = CrisprPlanner.get_dsb(self.pam_sites, self.mutation_direction)
        return SequenceSites(min(lowest_site, max(0, DSB - 20)), DSB - 1)  # closed segment

    @staticmethod
    def is_within_mutation_zone(mutation_zone, place):
        if mutation_zone[0] < place < mutation_zone[1]:
            return True
        else:
            return False

    # extracts the part of original and mutated sequences and checks to see whether one of them has a restriction site
    # the other one doesn't. returns the distinctive restriction sites if the lists don't match, and a None value
    # if they do
    def get_distinctive_restriction_sites(self, rest_site_type: RestrictionSiteType):
        # first we'll find the boundaries in which to search
        lowest_site = min(list(map(lambda item: item.index, self.mutated_sites)))
        highest_site = max(list(map(lambda item: item.index, self.mutated_sites)))
        # max_restriction_site = max(list(map(lambda item: len(item.site), self.restriction_enzymes)))
        max_restriction_site = 30  # prior analysis

        original_sequence = self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand
        original_section = original_sequence[max(0, lowest_site - max_restriction_site):
                                             min(len(original_sequence) - 1, highest_site + max_restriction_site)]
        mutated_section = self.mutated_strand[max(0, lowest_site - max_restriction_site):
                                              min(len(self.mutated_strand) - 1, highest_site + max_restriction_site)]

        restrictions_sites_in_original_seq = self.find_restriction_sites(original_section,
                                                                         max(0, lowest_site - max_restriction_site),
                                                                         RestrictionSiteType.REMOVED)
        restrictions_sites_in_mutated_seq = self.find_restriction_sites(mutated_section,
                                                                        max(0, lowest_site - max_restriction_site),
                                                                        RestrictionSiteType.INSERTED)

        if self.do_restriction_sites_lists_match(restrictions_sites_in_original_seq,
                                                 restrictions_sites_in_mutated_seq):
            return None
        else:
            return self.find_extra_restriction_sites(restrictions_sites_in_original_seq,
                                                     restrictions_sites_in_mutated_seq,
                                                     rest_site_type)

    def get_restriction_sites(self):
        max_len_site = max(list(map(lambda item: len(item.site), self.restriction_enzymes)))

        start = max(0, self.mutation_zone.start - max_len_site)
        end = min(self.mutation_zone.end + 1 + max_len_site, len(self.mutated_strand) - 1)
        mutated_section = self.mutated_strand[start:end]
        return self.find_restriction_sites(mutated_section,
                                           start,
                                           RestrictionSiteType.REMOVED)

    # receives a RestrictionSite list and returns a RestrictionMutation list
    def from_site_to_mutation(self, restriction_sites):
        restriction_mutations = []
        for r_site in restriction_sites:
            mutated_strand = self.change_chars_in_string(self.mutated_strand, self.mutated_sites)
            r_mutation = RestrictionMutation(restriction_site=r_site,
                                             number_of_mutations=len(self.mutated_sites),
                                             mutated_sites=[],
                                             mutated_strand=mutated_strand,
                                             codon_mutations=self.codon_mutations,
                                             reattachment_mutations=self.reattachment_mutations,
                                             repair_template=self.get_repair_template(self.mutated_sites,
                                                                                      mutated_strand))
            restriction_mutations.append(r_mutation)
        return restriction_mutations

    # checks if the mutated strand so far (with codon and reattachment mutations) has different inserted
    # restriction sites than the original one
    def find_added_restriction_sites(self):
        print("INSERTED WITHOUT: Checking for already distinctive inserted rest sites...")
        restriction_sites_list = self.get_distinctive_restriction_sites(RestrictionSiteType.INSERTED)
        if restriction_sites_list:  # distinctive restriction site found
            print("Lists don't match!")  # different restriction sites
            # sorts and filters restriction sites
            valid_rest_sites_list = self.check_mutated_restriction_sites(restriction_sites_list)
            if valid_rest_sites_list:
                print("Inserted restriction sites found without extra mutations!")
                self.result.no_extra_inserted_mutations += self.from_site_to_mutation(valid_rest_sites_list)
                if len(self.result.no_extra_inserted_mutations) >= self.max_results:
                    print("Inserted", self.max_results, "restriction sites!")
                    self.result.success = True
        return self.result.success

    def add_restriction_sites(self):
        print("INSERT: Let's modify some sites!\nfirst, trying to add for mutation zone:", self.mutation_zone)
        inserted_restriction_mutations = self.get_new_restriction_site()
        if inserted_restriction_mutations:
            print("Inserted restriction mutations!", inserted_restriction_mutations)
            self.result.inserted_mutations += inserted_restriction_mutations
            if len(self.result.no_extra_inserted_mutations) + len(
                    self.result.inserted_mutations) >= self.max_results:
                print("Inserted", self.max_results, "restriction sites!")
                self.result.success = True
        return self.result.success

    def find_removed_restriction_sites(self):
        print("REMOVE WITHOUT: Checking for already distinctive removed rest sites...")
        restriction_sites_list = self.get_distinctive_restriction_sites(RestrictionSiteType.REMOVED)
        if restriction_sites_list:  # distinctive restriction site found
            print("lists don't match!")  # different restriction sites
            # sorts and filters restriction sites
            valid_rest_sites_list = self.check_mutated_restriction_sites(restriction_sites_list)
            if valid_rest_sites_list:
                print("Removed restriction sites found without extra mutations!")
                self.result.no_extra_removed_mutations += self.from_site_to_mutation(valid_rest_sites_list)
                if len(self.result.no_extra_removed_mutations) >= self.max_results:
                    print("Removed", self.max_results, "restriction sites!")
                    self.result.success = True
        return self.result.success

    def remove_existing_restriction_sites(self):
        print("REMOVE: Now let's try to actively remove...")
        removed_restriction_mutations = self.get_removed_restriction_sites()
        if removed_restriction_mutations:
            print("Removed restriction mutation found!", removed_restriction_mutations)
            self.result.removed_mutations += removed_restriction_mutations
            if len(self.result.no_extra_removed_mutations) + len(
                    self.result.removed_mutations) >= self.max_results:
                print("Removed", self.max_results, "restriction sites!")
                self.result.success = True
        return self.result.success

    def add_remove_restriction_sites(self, restriction_site_type):
        self.mutation_zone = self.get_mutations_zone()
        print("Mutation zone:", self.mutation_zone)
        if restriction_site_type == RestrictionSiteType.INSERTED_WITHOUT:
            # trying to passively find an inserted restriction site
            return self.find_added_restriction_sites()
        elif restriction_site_type == RestrictionSiteType.INSERTED:
            # trying to actively insert a restriction site
            return self.add_restriction_sites()
        elif restriction_site_type == RestrictionSiteType.REMOVED_WITHOUT:
            # trying to passively find a removed restriction site
            return self.find_removed_restriction_sites()
        else:  # trying to actively remove a restriction site
            return self.remove_existing_restriction_sites()

    # outdated
    def is_in_favourites(self, restriction_sites_found):
        for r_site in restriction_sites_found:
            if type(r_site) == RestrictionMutation:
                r_site = r_site.restriction_site
            if r_site.enzyme in self.restriction_enzymes:
                return True
        return False

    # receives a lists of RestrictionSite that have been mutated into insertion\deletion in prior steps.
    # filters out irrelevant sites (too close to same sites) and sorts
    # them.
    def check_mutated_restriction_sites(self, restriction_sites_list: list):
        for rest_site in restriction_sites_list[:]:
            if self.check_distance(rest_site) < 150:
                restriction_sites_list.remove(rest_site)
        # get read of duplicates
        restriction_sites_list = list(set(restriction_sites_list))
        restriction_sites_list.sort(key=lambda r_site: (-(abs(r_site.index.get_mean() - self.DSB)),
                                                        self.check_distance(r_site, vicinity=250),
                                                        - self.check_rareness(r_site, 250)),
                                    reverse=True)
        print("Found restriction sites after sort:", *restriction_sites_list, sep="\n")
        return restriction_sites_list

    # uses the mutation zone and extracts all restriction sites there. after filtering and sorting, tries to remove
    # all extracted restriction sites and returns a list of RestrictionMutation of sites removed
    def get_removed_restriction_sites(self):
        rest_sites = self.get_restriction_sites()
        if not rest_sites:
            print("No restriction sites in the relevant section")
            return None
        self.filter_out_restriction_sites_with_no_space(rest_sites)
        if rest_sites:
            # if there are still restriction sites left, we  will now sort them
            self.sort_restriction_sites(rest_sites)
            print("rest sites after sorting:", rest_sites)
            restriction_sites_removed = self.remove_restriction_sites(rest_sites)
            if restriction_sites_removed:  # succeeded in finding and removing a restriction site
                print("Restriction site can be removed! restriction mutations:", restriction_sites_removed)
                return restriction_sites_removed
            print("No restriction site to remove was found")
            return None
        else:
            print("No restriction sites left after filtering out")
            return None

    # receives (1)list of RestrictionSite and (2) the start of the ssODN mutation codon tries to insert silent
    # mutations to remove them, and returns a list of RestrictionMutation that tell what changes can be made to each
    # restriction site to remove it
    def remove_restriction_sites(self, rest_sites: list):
        print("### remove_restriction_sites ###")
        restriction_mutations = []
        # fourth_strand = self.mutated_strand
        fourth_mutated_sites = self.mutated_sites[:]
        for rest_site in rest_sites:
            print("for rest site:", rest_site)
            # it's not a reattachment section, but I need this namedTuple for the get_relevant_codons method
            section_to_mutate = MutateSection(1, DNASection.MUTATION_ZONE, rest_site.index)
            possible_mutations, _ = self.get_possible_codon_mutations(section_to_mutate)

            for index in range(len(possible_mutations)):
                print("for mutation: ", possible_mutations[index])
                # self.mutated_strand = fourth_strand
                self.mutated_sites = fourth_mutated_sites[:]
                point_mutations = self.get_point_mutations(possible_mutations, [index])
                self.mutated_sites.extend(point_mutations)
                self.restriction_site_mutations = point_mutations
                if self.is_restriction_site_well_removed(point_mutations, rest_site):
                    mutated_strand = self.change_chars_in_string(self.mutated_strand, self.restriction_site_mutations)
                    restriction_mutations.append(
                        RestrictionMutation(restriction_site=rest_site,
                                            number_of_mutations=len(point_mutations),
                                            mutated_sites=self.restriction_site_mutations,
                                            mutated_strand=mutated_strand,
                                            codon_mutations=self.codon_mutations,
                                            reattachment_mutations=self.reattachment_mutations,
                                            repair_template=self.get_repair_template(self.mutated_sites,
                                                                                     mutated_strand)))
                    print("Restriction site removed with:", point_mutations)
        return restriction_mutations

    # receives (1) point mutations and the restriction site we wish to remove, mutate a copy of the relevant strand
    # according to the point mutations, take the part in which the restriction site is, and check if after the mutations
    # the restriction site is removed and can't function as a derivative of another restriction site
    def is_restriction_site_well_removed(self, point_mutations, rest_site: RestrictionSite):
        mutated_strand = self.change_chars_in_string(self.mutated_strand, point_mutations)
        print("strand after restriction site removal mutation:", mutated_strand)
        start_index, end_index = rest_site.index.start, rest_site.index.end
        new_rest_site = mutated_strand[start_index:end_index + 1]
        if FileReader.is_derivative(rest_site.enzyme.site, new_rest_site):
            print("New rest site", new_rest_site, "is a derivative for", rest_site.enzyme.name, rest_site.enzyme.site)
            return False
        sequence = mutated_strand[max(0, start_index - len(new_rest_site)):
                                  min(end_index + len(new_rest_site) + 1, len(mutated_strand) - 1)]
        for derivative in rest_site.enzyme.derivatives:
            if derivative in sequence or self.get_complementary_sequence(derivative) in sequence:
                print("New rest site has been created accidentally or has always been there, rest site is cancelled")
                return False
        print("New rest site", new_rest_site, "is not a derivative for", rest_site.enzyme.name)
        return True

    # receives a list of restriction sites (RestrictionSite) and sorts them according to what restriction sites the lab
    # has, its distance from same type restriction sites and frequency along the sequence
    def sort_restriction_sites(self, rest_sites):
        # first, sort according to existing restriction sites in the lab, second, sort according to the distance to next
        # same type restriction site, and third
        rest_sites.sort(key=lambda r_site: (-(abs(r_site.index.get_mean() - self.DSB)),
                                            self.check_distance(r_site, vicinity=250),
                                            -self.check_rareness(r_site, 250)),
                        reverse=True)

    def check_rareness(self, restriction_site, distance=250):
        occurrences = []
        site_length = len(restriction_site.enzyme.site)
        mutated_section = self.mutated_strand[max(0, restriction_site.index.start - distance - site_length):
                                              min(len(self.mutated_strand) - 1,
                                                  restriction_site.index.end + distance + site_length)]
        for derivative in restriction_site.enzyme.derivatives:
            occurrences += [m.start() for m in re.finditer('(?=' + derivative + ')', mutated_section.upper())]
            occurrences += [m.start() for m in re.finditer('(?=' + self.get_complementary_sequence(derivative) + ')',
                                                           mutated_section.upper())]
        return len(occurrences)

    # receives a list of restriction sites and the strand that they are taken from, and looks for other same restriction
    # sites in the vicinity of certain nucleotides. if the vicinity is too close, remove the restriction site from the
    # list
    def filter_out_restriction_sites_with_no_space(self, rest_sites, vicinity: int = 150):
        for restriction_site in rest_sites[:]:
            if self.check_distance(restriction_site, vicinity) < vicinity:
                rest_sites.remove(restriction_site)

    # receives a restriction site and a strand, and checks in what distance there is another restriction site from the
    # same type
    def check_distance(self, restriction_site: RestrictionSite, vicinity: int = 150):
        occurrences = []
        distance = vicinity + 1
        site_length = len(restriction_site.enzyme.site)
        start_point = max(0, restriction_site.index.start - vicinity - site_length)
        end_point = min(len(self.mutated_strand) - 1, restriction_site.index.end + vicinity + site_length)
        mutated_section = self.mutated_strand[start_point:end_point]
        # print("current site being tested:", restriction_site)
        for derivative in restriction_site.enzyme.derivatives:
            occurrences += [m.start() for m in re.finditer('(?=' + derivative + ')', mutated_section.upper())]
            occurrences += [m.start() for m in
                            re.finditer('(?=' + self.get_complementary_sequence(derivative) + ')',
                                        mutated_section.upper())]
        # print("occurrences:", occurrences, ", start:", start_point)
        for occurrence in occurrences:
            if 0 < abs(start_point + occurrence - restriction_site.index.start) < distance:
                distance = abs(start_point + occurrence - restriction_site.index.start)
                # print("distance from other same rest site:", distance, "in index:", start_point + occurrence)
        return distance

    def get_mean_of_subtracted_lists(self, item: RestrictionMutation):
        def subtract_two_lists(lst1, lst2):
            return list(set(lst1) - set(lst2))

        only_new_mutated_sites = subtract_two_lists(item.mutated_sites, self.mutated_sites)
        if not only_new_mutated_sites:
            # if empty
            only_new_mutated_sites = self.mutated_sites
        return mean(list(map(lambda site: site.index, only_new_mutated_sites)))

    # gets a list of RestrictionMutation objects, sorts them and returns the list
    def get_new_restriction_site(self):
        possible_restriction_mutations = self.get_possible_restriction_mutations()
        print("\nAll possible restriction sites mutations:", *possible_restriction_mutations, sep="\n")

        if possible_restriction_mutations:
            possible_restriction_mutations.sort(key=lambda r_mutation:
            (-r_mutation.number_of_mutations,
             -abs(self.DSB - self.get_mean_of_subtracted_lists(r_mutation)),
             self.check_distance(r_mutation.restriction_site, vicinity=250),
             -self.check_rareness(r_mutation.restriction_site)), reverse=True)
            return possible_restriction_mutations
        else:
            print("No possible restriction mutations found")
            return None

    # receives indexes power set and returns three lists of indexes set: one with sets closest to DSB (up to 10 nt), one
    # up to 15 and one beyond.
    def into_three_groups(self, indexes_sets, mutations):
        up10, up15, up20 = [], [], []
        for index_set in indexes_sets:
            if not index_set:
                continue
            distance = 0
            for index in index_set:
                codon_data, _ = mutations[index]
                tmp_distance = abs(codon_data.codon_sites.get_mean() - self.DSB)
                if tmp_distance > distance:
                    distance = tmp_distance
            if distance <= 10:
                up10.append(index_set)
            elif distance <= 15:
                up15.append(index_set)
            else:
                up20.append(index_set)
        print("divided into three:", len(up10), len(up15), len(up20))
        print("up10:", up10, "up15:", up15, "up20:", up20, sep="\n")
        return [up10, up15, up20]

    def sort_indexes_sets(self, indexes_sets, possible_mutations):
        indexes_sets.sort(key=lambda index_set: (CrisprPlanner.get_group_distance_from_dsb(index_set,
                                                                                           possible_mutations,
                                                                                           self.pam_sites,
                                                                                           self.mutation_direction),
                                                 sum([possible_mutations[i][1].number_of_mutations for i in
                                                      index_set])))

    # sorts the subsets of mutations and divided into fixed sites of groups
    # outdated
    def chunk(self, indexes_sets, possible_mutations, size):
        self.sort_indexes_sets(indexes_sets, possible_mutations)
        print("sorted indexes sets:", indexes_sets)
        it = iter(indexes_sets)
        return iter(lambda: tuple(islice(it, size)), ())

    @staticmethod
    def get_num_of_independent_codons(possible_mutations):
        codons_info = set()
        for possible_mutation in possible_mutations:
            codon_info = possible_mutation[0]
            codons_info.add(codon_info)
        return len(codons_info)

    # try to replace get_possible_restriction_sites
    # aims to find an inserted restriction site and return a list of restriction mutation objects, each contains data to
    # achieve an inserted restriction site. sorting out those who were already found without special mutations added
    def get_possible_restriction_mutations(self):
        valid_restriction_mutations = []
        section_to_mutate = MutateSection(1, DNASection.MUTATION_ZONE, self.mutation_zone)
        possible_mutations, _ = self.get_possible_codon_mutations(section_to_mutate)
        index_power_set = CrisprPlanner.get_power_set(lst=list(range(len(possible_mutations))),
                                                      max_num=self.get_num_of_independent_codons(possible_mutations))
        if [] in index_power_set:
            index_power_set.remove([])
        print("Length of power set:", len(index_power_set))
        for subset in index_power_set[:]:
            # check that mutation in subset don't override each other
            if CrisprPlanner.are_codons_overriden(subset, possible_mutations):
                index_power_set.remove(subset)
        print("Length of power set after filtration:", len(index_power_set))

        # groups = self.into_three_groups(index_power_set, possible_mutations)
        # groups = list(self.chunk(index_power_set, possible_mutations, 3))
        print("Sorting power set...")
        self.sort_indexes_sets(index_power_set, possible_mutations)
        print("Power set sorted!")

        original_strand = self.mutated_strand
        original_mutated_sites = self.mutated_sites[:]
        # print("Original mutated sites:", original_mutated_sites)
        for index_subset in index_power_set:

            self.mutated_strand = original_strand
            self.mutated_sites = original_mutated_sites[:]
            print("original mutated sites:", original_mutated_sites)
            point_mutations = self.get_point_mutations(possible_mutations, index_subset)
            print("For restriction site, now working on index subset:", index_subset, "and point mutations:",
                  point_mutations)
            self.mutated_strand = self.change_chars_in_string(self.mutated_strand, point_mutations)
            self.mutated_sites.extend(point_mutations)
            self.restriction_site_mutations = point_mutations
            # print("all mutated sites:", self.mutated_sites)
            restriction_sites = self.get_distinctive_restriction_sites(RestrictionSiteType.INSERTED)
            if not restriction_sites:
                continue
            for restriction_site in restriction_sites:
                if self.check_distance(restriction_site) > 150 and self. \
                        is_restriction_site_new(restriction_site, valid_restriction_mutations):
                    print("PASSED: rest site:", restriction_site)
                    mutated_strand = self.change_chars_in_string(self.mutated_strand,
                                                                 self.restriction_site_mutations)
                    valid_restriction_mutations.append(
                        RestrictionMutation(restriction_site,
                                            len(index_subset),
                                            self.restriction_site_mutations,
                                            mutated_strand,
                                            self.codon_mutations,
                                            self.reattachment_mutations,
                                            self.get_repair_template(self.mutated_sites,
                                                                     mutated_strand)))
                else:
                    print("FAILED: rest site:", restriction_site)
            print("valid restriction mutations:", len(valid_restriction_mutations))
            if len(valid_restriction_mutations) + len(self.result.inserted_mutations) + len(
                    self.result.no_extra_inserted_mutations) >= self.max_results:
                return valid_restriction_mutations
        return valid_restriction_mutations

    def is_restriction_site_new(self, restriction_site, valid_restriction_mutations):
        for r_mutation in valid_restriction_mutations:
            if restriction_site == r_mutation.restriction_site:
                return False
        for r_mutation in self.result.no_extra_inserted_mutations:
            if restriction_site == r_mutation.restriction_site:
                return False
        for r_mutation in self.result.inserted_mutations:
            if restriction_site == r_mutation.restriction_site:
                return False
        return True

    # this function extracts all possible codons in the mutation zone section and tries to insert silent mutation in
    # them, recursively, and collects all options that add a new restriction site to the strand. number of max_mutations
    # is chosen
    # out of use function
    # def get_possible_restriction_sites(self, mutated_strand, mutated_sites, ssODN_mutation_codon,
    #                                    mutations_so_far: int = 0):
    #     possible_results = []
    #     section_to_mutate = MutateSection(-1,
    #                                       DNASection.MUTATION_ZONE,
    #                                       self.mutation_zone)
    #     codons_data = CrisprPlanner.get_relevant_codons(section_to_mutate, ssODN_mutation_codon, mutated_strand)
    #     print("mutated sites so far:", mutated_sites)
    #     print("all possible codons for restriction mutations:", [data.codon for data in codons_data])
    #     for codon_data in codons_data:
    #         print("for codon:", codon_data.codon, "at", codon_data.codon_sites)
    #         # all codons_data that translate to the same amino acid as the relevant codon does
    #         same_aa_codons = list(CrisprPlanner.get_similar_codons(codon_data.codon))
    #         print("same codons:", same_aa_codons)
    #         CrisprPlanner.check_already_mutated_sites(codon_data, same_aa_codons, mutated_sites)
    #         # removes codons in which the difference is in part that in outside of the mutation zone
    #         CrisprPlanner.check_outside_codon_mutations(codon_data, same_aa_codons, section_to_mutate)
    #
    #         for same_aa_codon in same_aa_codons:
    #             num_of_mutations = mutations_so_far
    #             mutated_codon_details = CrisprPlanner.how_to_get_b_codons_from_a(codon_data.codon, [same_aa_codon])
    #             print("adding mutation:", mutated_codon_details, codon_data.codon_sites)
    #             dict_of_changes = mutated_codon_details.dict_of_mutations
    #             current_mutated_strand = mutated_strand
    #             temporary_point_mutations = []
    #             for key in dict_of_changes:
    #                 current_mutated_strand = CrisprPlanner.change_char_in_string(current_mutated_strand,
    #                                                                              codon_data.codon_sites.start + key,
    #                                                                              dict_of_changes[key])
    #                 temporary_point_mutations.append(
    #                     PointMutation(codon_data.codon_sites.start + key, dict_of_changes[key]))
    #                 num_of_mutations += 1
    #             restriction_sites = self.get_distinctive_restriction_sites(mutated_sites + temporary_point_mutations,
    #                                                                        current_mutated_strand)
    #             if restriction_sites:
    #                 print("restriction sites found:", restriction_sites, codon_data, mutated_codon_details)
    #                 for restriction_site in restriction_sites:
    #                     if self.check_distance(restriction_site) > 150 and \
    #                                     restriction_site.rest_site_type == RestrictionSiteType.INSERTED:
    #                         possible_results.append(RestrictionMutation(restriction_site,
    #                                                                     num_of_mutations,
    #                                                                     mutated_sites + temporary_point_mutations))
    #             print("going recursive!", num_of_mutations)
    #             possible_results.extend(self.get_possible_restriction_sites(current_mutated_strand,
    #                                                                         mutated_sites + temporary_point_mutations,
    #                                                                         ssODN_mutation_codon,
    #                                                                         num_of_mutations))
    #     return possible_results

    # gets two dictionaries of the format of RestrictionSite:RestrictionSiteType and checks if both dictionaries have the same restriction sites
    @staticmethod
    def do_restriction_sites_lists_match(sites_in_original_list, sites_in_mutated_list):
        if len(sites_in_original_list) != len(sites_in_mutated_list):
            return False
        else:
            for restriction_site in sites_in_original_list:
                filtered_list = list(filter(lambda item: item.index == restriction_site.index and
                                                         item.enzyme == restriction_site.enzyme, sites_in_mutated_list))
                # if list is empty, no restriction sites in that index
                if not filtered_list:
                    return False
        return True

    # receives the original sequence (with only codon mutations) and the optional sequence with more mutations that has
    # more\less restrictions sites than the original, and returns the added restriction site and indexes, and whether
    # the change is added restriction site or removed. format: (2, {'TCGA': (19, 23), 'CCGG': (0, 4)})
    @staticmethod
    def find_extra_restriction_sites(original_seq_restriction_sites, mutant_seq_restriction_sites,
                                     rest_site_type: RestrictionSiteType):
        extra_restrictions_sites = []

        # restriction site was added
        if rest_site_type == RestrictionSiteType.INSERTED:
            for restriction_site in mutant_seq_restriction_sites:
                filtered_list = list(filter(lambda item: item.index == restriction_site.index and
                                                         item.enzyme == restriction_site.enzyme,
                                            original_seq_restriction_sites))
                if not filtered_list:
                    extra_restrictions_sites.append(restriction_site)

        # restriction site was removed
        if rest_site_type == RestrictionSiteType.REMOVED:
            for restriction_site in original_seq_restriction_sites:
                filtered_list = list(filter(lambda item: item.index == restriction_site.index and
                                                         item.enzyme == restriction_site.enzyme,
                                            mutant_seq_restriction_sites))
                if not filtered_list:
                    extra_restrictions_sites.append(restriction_site)

        print("Found restriction sites:", *extra_restrictions_sites, sep="\n")
        return extra_restrictions_sites

    def remove_utr(self):
        seq_without_utr = ''
        for letter in self.sense_strand:
            if 'A' <= letter <= 'Z':
                seq_without_utr += letter
            else:
                continue
        return seq_without_utr

    def find_sites(self, seq):
        sites = []
        print("Finding sites in seq:", seq)
        for i in range(len(seq)):
            for j in range(i + 1, len(seq) + 1):
                for enzyme in self.restriction_enzymes:
                    for derivative in enzyme.derivatives:
                        if seq[i:j] == derivative:
                            print(enzyme.name, enzyme.site, derivative, i, j)
                            sites.append((enzyme.name, enzyme.site, derivative, i, j))
        return sites

    def get_favourite_restriction_enzymes(self, favourite_enzymes_names: list):
        all_enzymes = FileReader(FileReader.research_path,
                                 r"/Code/CRISPR/parsed_restriction_enzymes.txt").get_parsed_restriction_enzymes_list()
        if not favourite_enzymes_names:
            return all_enzymes
        favourite_enzymes = list(
            filter(lambda r_enzyme: r_enzyme.name in favourite_enzymes_names, all_enzymes))
        print("favourite enzymes:", *favourite_enzymes, sep="\n")
        if not favourite_enzymes:
            return all_enzymes
        self.result.enzymes = favourite_enzymes
        return favourite_enzymes

    def get_repair_template(self, point_mutations, mutated_strand, arm_length=35):
        min_index, max_index = self.get_min_max_index_of_point_mutations(point_mutations)
        if min_index > max_index:
            print("Could not find minimum and maximum point mutations")
            return [None, None, None]
        mutated_section = mutated_strand[min_index:max_index + 1]
        if min_index >= arm_length and max_index <= len(self.sense_strand) - 1 - arm_length:
            left_arm = self.mutated_strand[min_index - arm_length:min_index]
            right_arm = self.mutated_strand[max_index + 1:max_index + 1 + arm_length]
            return [left_arm, mutated_section, right_arm]
        else:  # need to extract fuller sequence
            print("extracting sequence with padding...")
            seq_with_padding = HttpRequester().get_transcript(gene_id=self.gene_id, with_padding=True)
            if not seq_with_padding:
                return [None, None, None]
            index = seq_with_padding.find(self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand)
            if index < 0:
                seq_with_padding = self.get_complementary_sequence(seq_with_padding)
                index = seq_with_padding.find(self.sense_strand if self.ssODN_direction > 0 else self.anti_sense_strand)
                if index < 0:
                    return [None, None, None]
            print("index:", index)
            min_index += index
            max_index += index
            print("min:", min_index, "max:", max_index)
            if min_index >= arm_length and max_index <= len(seq_with_padding) - 1 - arm_length:
                left_arm = seq_with_padding[min_index - arm_length:min_index]
                right_arm = seq_with_padding[max_index + 1:max_index + 1 + arm_length]
                return [left_arm, mutated_section, right_arm]
            else:
                return [None, None, None]

    def get_min_max_index_of_point_mutations(self, point_mutations):
        max_index = 0
        min_index = len(self.sense_strand) - 1
        for point_mutation in point_mutations:
            if point_mutation.index > max_index:
                max_index = point_mutation.index
            if point_mutation.index < min_index:
                min_index = point_mutation.index
        print("min", min_index, "max:", max_index)
        return min_index, max_index

    @staticmethod
    def get_formatted_crrna(crrna, strand, crrna_length=20):
        new_sequence = crrna[0]  # potential crrna
        sequence_sites = crrna[1]
        new_start = len(strand) - sequence_sites.start - crrna_length
        new_sequence_sites = SequenceSites(new_start, new_start + crrna_length - 1)
        new_crrna = (new_sequence, new_sequence_sites)
        return new_crrna

    @staticmethod
    def get_formatted_point_mutations(point_mutations, strand):
        new_point_mutations = []
        for point_mutation in point_mutations:
            new_index = len(strand) - point_mutation.index - 1
            new_point_mutation = PointMutation(new_index,
                                               point_mutation.old_nucleotide,
                                               point_mutation.new_nucleotide)
            new_point_mutations.append(new_point_mutation)
        return new_point_mutations

    @staticmethod
    def get_formatted_restriction_site(restriction_site: RestrictionSite, strand):
        sequence_sites = restriction_site.index
        new_sequence_sites = SequenceSites(len(strand)-sequence_sites.end-1,
                                           len(strand)-sequence_sites.start-1)
        new_restriction_site = RestrictionSite(index=new_sequence_sites,
                                               enzyme=restriction_site.enzyme,
                                               rest_site_type=restriction_site.rest_site_type)
        return new_restriction_site

    @staticmethod
    def get_formatted_restriction_mutations(restriction_mutations, cp, strand):
        new_restriction_mutations = []
        for rest_mutation in restriction_mutations:
            new_restriction_site = cp.get_formatted_restriction_site(rest_mutation.restriction_site, strand)
            new_mutated_sites = cp.get_formatted_point_mutations(rest_mutation.mutated_sites, strand)
            new_codon_mutations = cp.get_formatted_point_mutations(rest_mutation.codon_mutations, strand)
            new_reattachment_mutations = cp.get_formatted_point_mutations(rest_mutation.reattachment_mutations, strand)
            new_restriction_mutation = RestrictionMutation(restriction_site=new_restriction_site,
                                                           number_of_mutations=rest_mutation.number_of_mutations,
                                                           mutated_sites=new_mutated_sites,
                                                           mutated_strand=rest_mutation.mutated_strand,
                                                           codon_mutations=new_codon_mutations,
                                                           reattachment_mutations=new_reattachment_mutations,
                                                           repair_template=rest_mutation.repair_template)
            new_restriction_mutations.append(new_restriction_mutation)
        return new_restriction_mutations


    @staticmethod
    def change_index_to_match_sense(result: Result):
        cp = CrisprPlanner
        new_anti_sense_crrnas = []
        for crrna in result.anti_sense_crrnas:
            new_crrna = cp.get_formatted_crrna(crrna, result.sense_strand)
            new_anti_sense_crrnas.append(new_crrna)
        result.anti_sense_crrnas = new_anti_sense_crrnas
        if result.crRNA_strand < 0:  # chosen crRNA is from anti-sense
            new_crrna = cp.get_formatted_crrna(result.crRNA, result.sense_strand)
            result.crRNA = new_crrna
            result.pam_sites = SequenceSites(new_crrna[1].start-3, new_crrna[1].start-1)
        if result.ssODN_strand < 0:
            new_no_extra_inserted = cp.get_formatted_restriction_mutations(result.no_extra_inserted_mutations, cp,
                                                                           result.sense_strand)
            result.no_extra_inserted_mutations = new_no_extra_inserted
            new_inserted = cp.get_formatted_restriction_mutations(result.inserted_mutations, cp,
                                                                           result.sense_strand)
            result.inserted_mutations = new_inserted
            new_no_extra_removed = cp.get_formatted_restriction_mutations(result.no_extra_removed_mutations, cp,
                                                                           result.sense_strand)
            result.no_extra_removed_mutations = new_no_extra_removed
            new_removed = cp.get_formatted_restriction_mutations(result.removed_mutations, cp,
                                                                           result.sense_strand)
            result.removed_mutations = new_removed