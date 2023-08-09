from Code.Enum.FileMode import FileMode
from Code.Enum.FileType import FileType


class FileMaker:

    # receives a set and a desired file name and writes all items in the set to a new file with said name
    @staticmethod
    def fromSetToFile(s: set, fileName: str):
        f = open(fileName, FileMode.WRITE.value, encoding="utf-8")
        for value in s:
            f.write(value)
            f.write("\n")

        f.close()

    @staticmethod
    def fromListToFile(l: list, file_name: str):
        f = open(file_name, FileMode.WRITE.value, encoding="utf-8")
        for value in l:
            f.write(value + "\n")
        f.close()

    # receives a dict and a file name and write it down to a new file with said name
    @staticmethod
    def fromDictToUniqueFile(d: dict, fileName: str):
        f = open(fileName, FileMode.WRITE.value, encoding="utf-8")
        d_keys = list(d.keys())
        for key in d_keys:
            f.write(key + "\t")
            values = d[key]
            for i in range(len(values)):
                if isinstance(values[i], tuple) and i == 0:
                    f.write(values[i][0] + ":" + str(values[i][1]))
                elif isinstance(values[i], tuple) and i > 0:
                    f.write(", " + values[i][0] + ":" + str(values[i][1]))
                else:
                    f.write('\t')
                    f.write(values[i])
            f.write("\n")

        f.close()

    # receives a simple key-value dictionary and writes it to a file
    @staticmethod
    def from_dict_to_file(d: dict, file_name: str):
        f = open(file_name, FileMode.WRITE.value, encoding="utf-8")
        for key in d.keys():
            f.write(str(key) + "\t" + str(d[key]) + "\n")
        f.close()

    # receives two dict and combines them to a new file
    @staticmethod
    def fromTwoDictToFile(d1: dict, d2: dict, file_name: str):
        f = open(file_name, FileMode.WRITE.value, encoding="utf-8")
        keys = list(d1.keys())
        for key in keys:
            f.write(key + "\t" + str(d1[key]).strip('[]') + "\t" + d2[key] + "\n")
        f.close()

    @staticmethod
    def from_plural_valued_dict_to_file(d: dict, file_name: str):
        f = open(file_name, FileMode.WRITE.value)
        keys = d.keys()
        for key in keys:
            values = ", ".join(d[key])
            f.write(key + "\t" + values + "\n")
        f.close()

    @staticmethod
    def fromOneFileToMany(file_name, new_file_name, lines_in_file: int):
        lines_index = 0
        files_index = 0
        out_file = open(file_name)
        lines = out_file.readlines()
        out_file.close()
        length = len(lines)
        while lines_index < length:
            in_file = open(new_file_name + str(files_index), FileMode.WRITE.value)
            print("file number " + str(files_index) + " has opened!")
            files_index += 1
            for i in range(lines_in_file):
                try:
                    in_file.write(lines[lines_index])
                    lines_index += 1
                except IndexError:
                    in_file.close()
                    return
            print("lines index: " + str(lines_index))
            in_file.close()

    @staticmethod
    def from_two_dict_of_different_keys_to_one_file(d1: dict, d2: dict, converter: dict, file_name: str):
        f = open(file_name, FileMode.WRITE.value)
        keys = d1.keys()
        for d1_key in keys:
            try:
                d2_key = converter[d1_key]
            except Exception as e:
                print("Exception in from_two_dict_of_different_keys_to_one_file:", e)
                print(d1_key + " doesn't have an equivalent in the second dictionary")
                pass
            d1_val = d1[d1_key]
            try:
                d2_val = d2[d2_key]
            except Exception as e:
                print("Exception in from_two_dict_of_different_keys_to_one_file:", e)
                print(d2_key + " doesn't exist as a key in d2")
                pass
            f.write(d1_key + "\t" + d2_key + "\t" + str(d1_val).strip("[]") + "\t" + d2_val + "\n")
        f.close()

    @staticmethod
    def fixTabbedFile(file_path, new_file_path):
        f = open(file_path)
        new_f = open(new_file_path, FileMode.WRITE.value)
        for line in f:
            lst = line.strip("\n").split("\t")
            gene = lst[0]
            gene_id = "".join(lst[2:])
            new_f.write(gene + "\t" + gene_id + "\n")
        f.close()
        new_f.close()

    @staticmethod
    def separateOrthologsInfoInData(data_file_path, new_data_file_path):
        in_file = open(data_file_path)
        out_file = open(new_data_file_path, FileMode.WRITE.value)
        for line in in_file:
            lst = line.strip("\n").split("\t")
            gene_id = lst[0]
            gene_name = lst[1]
            orthologs_info = lst[2].split(", ")
            phenotype = lst[3]
            description = lst[4]
            out_file.write(gene_id + "\t" + gene_name + "\t" + orthologs_info[0] + "\t" + orthologs_info[1]
                           + "\t" + orthologs_info[2] + "\t" + phenotype + "\t" + description + "\n")
        in_file.close()
        out_file.close()

    # receives (1) a file type which differentiate writing to console from writing to a file, (2) data to be writen and
    # (3) the file if there is one, and writes to the file or to the console the data necessary
    @staticmethod
    def write_to(file_type: FileType, data, f = ""):
        if file_type == FileType.CONSOLE:
            print(data)
        else:
            f.write(data)



