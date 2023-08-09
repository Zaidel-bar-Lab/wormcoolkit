import requests
from requests.exceptions import HTTPError
from base64 import b64encode
import json
from urllib import request, parse

from Code.Files.FileReader import FileReader
from Code.Utils.Strings import Strings


class HttpRequester:

    def __init__(self, url=""):
        self.url = url
        self.id_upi_converter = FileReader(FileReader.research_path + r"/Data",
                                           r"/human_gene_id_upi.txt").from_file_to_dict_with_plural_values(0, 1, True)
        # self.id_upa_converter = FileReader(FileReader.research_path + r"\Data",
        #                                    r"\human_gene_id_upa.txt").from_file_to_dict_with_plural_values(0, 1, True)

    def make_request(self):
        # sending get request and saving the response as response object
        try:
            response = requests.get(url=self.url)
            # If the response was successful, no Exception will be raised
            response.raise_for_status()
        except HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')
            return None
        except Exception as err:
            print(f'Other error occurred: {err}')
            return None
        else:
            # success
            return str(response.content)

    @staticmethod
    def get_human_uniprot_html(gene_id):
        try:
            r = requests.get("https://www.uniprot.org/uniprot/?query=" + gene_id +
                             "&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+reviewed%3Ayes&sort=score")
        except Exception as e:
            print("Exception in get_human_uniprot_html:", e)
            print("Communication failure, please check your internet connection")
            return None
        if not r.ok:
            # r.raise_for_status()
            print("Something went wrong with get_human_uniprot_html while trying to extract reviewed ids, please check")
            return None
        return r.text

    @staticmethod
    def get_worm_uniprot_html(gene_id):
        try:
            r = requests.get("https://www.uniprot.org/uniprot/?query=" + gene_id +
                             "&fil=organism%3A%22Caenorhabditis+elegans+%5B6239%5D%22+AND+reviewed%3Ayes&sort=score")
        except Exception as e:
            print("Exception in get_worm_uniprot_html:", e)
            print("Communication failure, please check your internet connection")
            return None
        if not r.ok:
            # r.raise_for_status()
            print("Something went wrong with get_worm_uniprot_html while trying to extract reviewed ids, please check")
            return None
        return r.text

    def get_longest_human_protein_sequence_from_uniprot(self, gene_id):
        chosen_seq = ''
        if gene_id not in self.id_upi_converter:
            return None
        for upi in self.id_upi_converter[gene_id]:
            if upi == "":
                continue
            request_url = "https://www.ebi.ac.uk/proteins/api/uniparc/upi/" + upi + "?rfTaxId=9606"
            try:
                r = requests.get(request_url, headers={"Accept": "text/x-fasta"})
            except Exception as e:
                print("Exception in get_longest_human_protein_sequence_from_uniprot:", e)
                print("Communication failure, please check your internet connection")
                return None
            if not r.ok:
                # r.raise_for_status()
                print(
                    "Something went wrong with get_longest_human_protein_sequence_from_uniprot() while trying to extract sequence"
                    ", please check")
                return None

            response_body = r.text
            optional_seq = Strings.from_fasta_seq_to_seq(response_body)
            if len(optional_seq) > len(chosen_seq):
                chosen_seq = optional_seq
        return chosen_seq

    def get_protein_sequence_from_ensembl(self, gene_id):
        if not gene_id:
            return None
        request_url = self.url + gene_id + "?type=protein;multiple_sequences=1"
        try:
            r = requests.get(request_url, headers={"Accept": "text/x-fasta"})
        except Exception as e:
            print("Communication failure, please check your internet connection:", e)
            return None
        if not r.ok:
            # r.raise_for_status()
            print(
                "Something went wrong with get_longest_human_protein_sequence_from_uniprot() while trying to extract sequence"
                ", please check")
            return None
        return r.text

    # returns a transcript of the sense strand
    @staticmethod
    def get_transcript(gene_id, isoform ,type , result=None, with_padding=False):
        print("def get_transcript")
        chosen_transcript = ''
        chosen_transcript_id = ''
        print("isoform: ", isoform)
        if (isoform != None):
            transcript = HttpRequester(url="http://rest.wormbase.org/rest/field/transcript/"). \
                get_transcript_by_transcript_id(isoform,type, with_padding)
            chosen_transcript = transcript
            chosen_transcript_id = isoform
            print("chosen_transcript: ", chosen_transcript)
            print("chosen_transcript_id: ", chosen_transcript_id)
            if chosen_transcript == None:
                raise Exception(
                    "Check the isoform field again, it seems to be incorrect. If it is indeed written correctly, NCBI cannot provide the gene sequence now,\n insert it yourself in the appropriate field or try again in a few minutes.")
        else:
            transcript_ids = HttpRequester(url="http://rest.wormbase.org/rest/widget/gene/"). \
                get_list_of_transcript_ids(gene_id)
            print("transcript_ids1 ", transcript_ids)
            if transcript_ids == None:
                raise Exception(
                    "Check the gene field again, it seems to be incorrect. If it is indeed written correctly, NCBI cannot provide the gene sequence now,\n insert it yourself in the appropriate field or try again in a few minutes.")
            if not transcript_ids or not len(transcript_ids):
                return None
            print("transcript_ids2 ", transcript_ids)
            for transcript_id in transcript_ids:
                transcript = HttpRequester(url="http://rest.wormbase.org/rest/field/transcript/"). \
                    get_transcript_by_transcript_id(transcript_id,type, with_padding)
                if transcript and len(transcript) > len(chosen_transcript):
                    chosen_transcript = transcript
                    chosen_transcript_id = transcript_id
                    print(transcript_id, "is better:", len(transcript))
            print("chosen_transcript: ", chosen_transcript)
            print("chosen_transcript_id: ", chosen_transcript_id)
            if not chosen_transcript:
                print("Couldn't find transcript")
                return None
            if not with_padding:
                result.request_url = "http://rest.wormbase.org/rest/field/transcript/" + chosen_transcript_id + "/unspliced_sequence_context"
            print("final", chosen_transcript)
        return (chosen_transcript, chosen_transcript_id)

    def get_transcript_by_transcript_id(self, transcript_id ,type, with_padding=False):
        if not transcript_id:
            return None
        with_padding = True
        print("trying!!!!!!!!!!!!!!!!")
        if with_padding or type == "insertion":
            request_url = self.url + transcript_id + "/unspliced_sequence_context_with_padding"
        else:
            request_url = self.url + transcript_id + "/unspliced_sequence_context"
        try:
            r = requests.get(request_url)
        except Exception as e:
            print("Communication failure in get_transcript_by_gene_id, please check your internet connection:", e)
            return None
        if not r.ok:
            print("Something went wrong with get_transcript_by_gene_id while trying to extract transcript"
                  ", please check")
            # r.raise_for_status()
            return None
        try:
            if with_padding or type == "insertion":
                sequence = r.json()['unspliced_sequence_context_with_padding']['data']['positive_strand']['sequence']
            else:
                sequence = r.json()['unspliced_sequence_context']['data']['positive_strand']['sequence']
        except Exception as e:
            print("Couldn't extract sequence from worm base:", e)
            return None
        return sequence

    def get_list_of_transcript_ids(self, gene_id):
        transcript_ids = set()
        # url = "http://rest.wormbase.org/rest/widget/gene/"
        if not gene_id:
            return None
        request_url = self.url + gene_id + "/sequences"
        try:
            r = requests.get(request_url)
        except Exception as e:
            print("Communication failure in get_list_of_transcript_ids, please check your internet connection:", e)
            return None
        if not r.ok:
            # r.raise_for_status()
            print("Something went wrong with get_list_of_transcript_ids while trying to extract transcript ids"
                  ", please check")
            return None
        try:
            table = r.json()['fields']['gene_models']['data']['table']
            for table_item in table:
                records = table_item['model']
                if type(records) == list:
                    for record in records:
                        transcript_ids.add(record['label'])
                if type(records) == dict:  # one dict, no list of dicts
                    transcript_ids.add(records['label'])
        except Exception as e:
            print("Error in get_list_of_transcript_ids: couldn't extract transcript ids:", e)
            return None
        print("list of transcript ids:", *transcript_ids, sep="\n")
        return transcript_ids

    @staticmethod
    def get_protein_seq_by_uniprot_swissprot_id(uniprot_swissprot_id):
        try:
            r = requests.get("https://www.uniprot.org/uniprot/" + uniprot_swissprot_id + ".fasta")
            response_body = r.text
        except Exception as e:
            print("Problem in function get_protein_seq_by_uniprot_swissprot_id while trying to extract sequence for",
                  uniprot_swissprot_id, ":", e)
            return None
        seq = Strings.from_fasta_seq_to_seq(response_body)
        #print(seq)
        return seq if seq else None

    @staticmethod
    def check_crRNA(sequences):
        default_result = {seq: (None, None) for seq in sequences}
        try:
            token = HttpRequester.get_authentication(client_id="wormcoolkit",
                                                     client_secret='ef91a593-8f12-4085-a458-2f6f491bb149',
                                                     idt_username="wormcoolkit",
                                                     idt_password="rzblab2020")
        except Exception as e:
            # authentication failed
            return default_result, e
        if not token:
            # authentication failed
            return default_result, "Token could not be extracted: authentication process failed"
        results_values = {}
        sequences_list = []
        format_sequence_dic = dict(zip([seq.upper() for seq in sequences], sequences))
        print("format dic:", format_sequence_dic)
        for sequence in sequences:
            sequence_dic = {"Name": "seq-" + sequence, "Sequence": sequence}
            sequences_list.append(sequence_dic)
        url = 'https://eu.idtdna.com/restapi/v1/Crispr/Check'
        headers = {
            'Content-Type': 'application/json',
            'Accept': 'application/json',
            'Authorization': 'Bearer ' + token,
        }
        data = {"species": "ROUNDWORM", "data": sequences_list, "libraryType": "CAS9"}
        r = requests.post(url, headers=headers, data=str(data))
        if not r.ok:
            return default_result, "Could not carry the request" + r.reason
        r_json = r.json()
        if r_json['hasModelErrors']:
            return default_result, r_json['modelErrors']
        results = r_json['output']['crisprResults']
        try:
            for result in results:
                input_seq = result['input']['data'][0]['Sequence']
                designs = result['output']['designs']
                for design in designs:
                    try:
                        on_target_score = design['onTargetPotential']
                    except KeyError:
                        on_target_score = None
                    try:
                        off_target_score = design['offTargetRiskSpecificity']
                    except KeyError:
                        off_target_score = None
                    results_values[format_sequence_dic[input_seq]] = (round(on_target_score), round(off_target_score))
            return results_values, None
        except:
            return default_result, "Could not extract on-target and\\or off-target scores for sequences"

    @staticmethod
    def get_authentication(client_id, client_secret, idt_username, idt_password):
        """
        Create the HTTP request and send it, parsing the response for the access token.

        The body_dict will also contain the fields "expires_in" (value of 3600) and "token_type" (value of "Bearer").

        If the request fails for any reason, an exception will be thrown that contains debugging information
        """
        authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
        request_headers = {"Content-Type": "application/x-www-form-urlencoded",
                           "Authorization": "Basic " + authorization_string}

        data_dict = {"grant_type": "password",
                     "scope": "test",
                     "username": idt_username,
                     "password": idt_password}
        request_data = parse.urlencode(data_dict).encode()

        post_request = request.Request("https://eu.idtdna.com/Identityserver/connect/token",
                                       data=request_data,
                                       headers=request_headers,
                                       method="POST")

        response = request.urlopen(post_request)
        body = response.read().decode()

        if response.status != 200:
            raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)

        body_dict = json.loads(body)
        if "access_token" in body_dict:
            return body_dict["access_token"]
        return None

    @staticmethod
    def get_authentication_test(client_id, client_secret, idt_username, idt_password):
        authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
        request_headers = {"Content-Type": "application/x-www-form-urlencoded",
                           "Authorization": "Basic " + authorization_string}

        data_dict = {"grant_type": "password",
                     "scope": "test",
                     "username": idt_username,
                     "password": idt_password}
        request_data = parse.urlencode(data_dict).encode()

        post_request = request.Request("https://eu.idtdna.com/Identityserver/connect/token",
                                       data=request_data,
                                       headers=request_headers,
                                       method="POST")

        print("Headers:{}\nData:{}\nMethod:{}\nURL:{}\n".format(post_request.headers, post_request.data,
                                                                post_request.method, post_request.full_url))
        response = request.urlopen(post_request)
        return response
