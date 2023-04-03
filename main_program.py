from Code.Enum.FileType import FileType
from Executors.executor import executor
import sys


def main(program):

    if program == 'Homology':
        ###### HOMOLOGOUS PIPELINE SECTION ########
        past_input = ['MDK', 'PTN', 'IFN', 'JAK1', 'JAK2', 'JAK3', 'JAK4', 'STAT3', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4',
                 'HRAS', 'KRAS', 'NRAS', 'RAF1', 'MEK1', 'MEK2', 'ERK1', 'ERK2', 'PI3K', 'AKT1', 'AKT2', 'AKT3', 'FOXO1',
                 'FOXO3', 'FOXO4', 'FOXO6', 'mTOR', 'RSK', 'INSR', 'IRS', 'PLCG1', 'PLCG2', 'DAG1', 'PKC', 'EGFR', 'GRB2',
                 'SOS1', 'SOS2', 'RAS', 'RAF', 'MEK', 'MAPK', 'CREB', 'MYC', 'EGF', 'BMP', 'SMAD', 'P38', 'JNK']
        former_input = ['IFNG', 'JAK1', 'JAK2', 'JAK3', 'MAP2K1', 'MAP2K2', 'MAPK3', 'MAPK1', 'PIK3CA', 'PIK3CG', 'PIK3CD',
                 'PIK3CB', 'RPS6KB1', 'RPS6KA3', 'RPS6KA1', 'RPS6KB2', 'RPS6KA2', 'RPS6KA6', 'RPS6KA5', 'IRS1', 'PRKCA',
                 'KRAS', 'NRAS', 'HRAS', 'MAPK1', 'MAPK3', 'CREB1', 'CREBBP', 'SMAD2', 'SMAD4', 'SMAD1']
        input = ['SMAD3', 'MAPK14', 'MAPK8']
        true_matches = executor.find_me_orthologs_for_human(input,
                                                            genes_in_names=True,
                                                            sources_bar=3,
                                                            length_bar=10,
                                                            domains_range=(0.3, 2))
        print("We are left with " + str(len(true_matches)) + " matches")
        print(str(true_matches))

    elif program == 'Variants':
        ### VARIANTS ###
        executor.get_variants_data(FileType.CONSOLE,
                                   {'TCP1': ['Asn284Ser', 'Ala453Glu']})
    else:
        print("No program called", sys.argv[1])
        exit()

# main(sys.argv[1])
main("Homology")
