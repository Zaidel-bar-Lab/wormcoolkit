from flask import Flask, render_template, request, send_from_directory

from Code.CRISPR.Enum.AminoAcid import AminoAcid
from Code.Utils.BioPython import BioPython

from Code.Http import HttpRequester
from Executors.executor import executor
from Code.CRISPR.CrisprPlanner import CrisprPlanner
import re
import os

app = Flask(__name__)

buffer = {
    "worm_gene_name": "",
    "nt_seq": "",
    "site": "",
    "given_isoform": "",
    "from_aa": "",
    "to_aa": "",
    "favourite_enzymes": "",
    "max_results": "",
    "reverse_linker" : "",
    "forward_linker" : "" ,
    "primer_reverse" : "" ,
    "primer_forward" : "",
    "insertion" : "",
    "end" : ""
    }


aa_names_dict = {'A':"ALANINE", 'R':"ARGININE",'N':"ASPARAGINE" ,'D':"ASPARTIC_ACID",'C':"CYSTEINE",
    'E':"GLUTAMIC_ACID",'Q':"GLUTAMINE",'G':"GLYCINE",'H':"HISTIDINE",'I':"ISOLEUCINE",'L':"LEUCINE",
    'K':"LYSINE",'M':"METHIONINE", 'F':"PHENYLALANINE",'P':"PROLINE",'S':"SERINE",'T':"THREONINE",
    'W':"TRYPTOPHAN",'Y':"TYROSINE",'V':"VALINE"}


@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')


@app.route("/")
def home():
    return render_template("home.html")


@app.route("/about")
def about():
    return render_template("about.html")


@app.route("/leo")
def salvador():
    return "Hello, Leo. Who's a good boy?"


@app.route('/c_elegans_orthologs')
def get_c_elegans_orthologs_input():
    return render_template('query_form.html')


@app.route('/c_elegans_orthologs', methods=['POST'])
def return_c_elegans_orthologs():
    human_genes = None
    try:
        text = request.form['text']
        text = text.replace(" ", "")
        human_genes = text.split(",")
        print("human genes:", human_genes)
        genes_in_names = request.form.get('type_select')
        sources_bar = request.form.get('sources_bar')
        length_ratio_down = float(request.form.get('length_ratio_from'))
        length_ratio_top = float(request.form.get('length_ratio_to'))
        print("length ratio:", (length_ratio_down, length_ratio_top))
        print("genes in names:", genes_in_names)
        results, error = executor.find_me_orthologs_for_human(human_genes=human_genes,
                                                              genes_in_names=genes_in_names,
                                                              sources_bar=sources_bar,
                                                              length_range=(length_ratio_down, length_ratio_top))
    except Exception as e:
        query = ", ".join(human_genes)
        error = "Something went wrong: " + str(e)
        return render_template('failure_response.html', query=query, error=error)
    true_results, false_results = results
    return render_template('orthologs_table_response.html', true_results=true_results, false_results=false_results)


@app.route('/human_orthologs')
def get_human_orthologs_input():
    return render_template('query_form.html')


@app.route('/human_orthologs', methods=['POST'])
def return_human_orthologs():
    try:
        text = request.form['text']
        text = text.replace(" ", "")
        c_elegans_genes = text.split(",")
        print("C.elegans genes:", c_elegans_genes)
        genes_in_names = request.form.get('type_select')
        sources_bar = request.form.get('sources_bar')
        length_ratio_down = float(request.form.get('length_ratio_from'))
        length_ratio_top = float(request.form.get('length_ratio_to'))
        print("genes in names:", genes_in_names)
        results, error = executor.find_me_orthologs_for_worm(worm_genes=c_elegans_genes,
                                                             genes_in_names=genes_in_names,
                                                             sources_bar=sources_bar,
                                                             length_range=(length_ratio_down, length_ratio_top))
    except Exception as e:
        query = ", ".join(c_elegans_genes)
        error = "Something went wrong: " + str(e)
        return render_template('failure_response.html', query=query, error=error)
    true_results, false_result = results
    return render_template('orthologs_table_response.html', true_results=true_results, false_results=false_result)


@app.route('/variants')
def get_variants_input():
    return render_template('query_form.html')


@app.route('/variants', methods=['POST'])
def return_variants_data():
    try:
        text = request.form['text']
        text = text.replace(" ", "")
        print(text)
        genes_and_variants = executor.parse_input(text)
        sources_bar = request.form.get('sources_bar')
        length_ratio_down = float(request.form.get('length_ratio_from'))
        length_ratio_top = float(request.form.get('length_ratio_to'))
        print(genes_and_variants)
        true_results, false_results = executor().get_variants_data_for_server(genes_and_variants,
                                                                              sources_bar,
                                                                              (length_ratio_down, length_ratio_top))
    except Exception as e:
        query = executor.dictionary_output_parser(genes_and_variants)
        return render_template('failure_response.html', query=query, error=e)
    return render_template('variants_table_response.html', true_results=true_results, false_results=false_results)


@app.route("/crispr")
def point_mutation_planner():
    return render_template("crispr_form.html")

@app.route("/crispr_insertion")
def insertion_planner():
    return render_template("crispr_form_2.html")


@app.route('/crispr', methods=['POST'])
def return_crispr_plan():
    global buffer

    try:
        crrna_done = request.form['crrna_done']
        if crrna_done == 'true':
            worm_gene_name = buffer['worm_gene_name']
            nt_seq = buffer['nt_seq']
            site = buffer['site']
            given_isoform = buffer['given_isoform']
            end = buffer['end']
            from_aa = buffer['from_aa']
            to_aa = buffer['to_aa']
            favourite_enzymes = buffer['favourite_enzymes']
            max_results = buffer['max_results']
            crrna = request.form['crrna']
            crrna_strand = int(request.form.get('crrna_strand'))

        else:
            print("crispr")
            worm_gene_name = request.form['name']
            nt_seq = request.form['seq']
            site = int(request.form['site'])
            given_isoform = request.form['given_isoform']
            from_aa = AminoAcid[request.form.get('from_aa')]

            to_aa = AminoAcid[request.form.get('to_aa')]
            favourite_enzymes = re.split('; |, | |,|;|\t|\n', request.form['enzymes']) if request.form['enzymes'] else None
            max_results = int(request.form.get('max_results'))
            crrna = request.form['crrna']
            crrna_strand = int(request.form.get('crrna_strand'))
            end = False


        print("CRISPR Request:", worm_gene_name,given_isoform, site, nt_seq, from_aa, to_aa, max_results, crrna, crrna_strand)

        result, error = CrisprPlanner(gene_name=worm_gene_name,
                                      aa_mutation_site=site,
                                      end = end,
                                      given_isoform=given_isoform,
                                      sense_strand=nt_seq,
                                      favourite_enzymes_names=favourite_enzymes,
                                      type= "point",
                                      max_results=max_results).plan_my_crispr(from_aa=from_aa,
                                                                              to_aa=to_aa,
                                                                              crrna=crrna,
                                                                              crrna_strand=crrna_strand)
        if not result or error:
            print(error)
            return render_template('failure_response.html', query=worm_gene_name, error=error)
    except Exception as e:
        error = "Something went wrong: " + str(e)
        print(error)
        return render_template('failure_response.html', query=worm_gene_name, error=error)
    # executor.increment_point_mutation_index(result)
    if crrna:
        print("Result:", result)
        CrisprPlanner.change_index_to_match_sense(result)
        print("Formatted result:", result)
        ssodn_strand = "sense" if result.ssODN_strand > 0 else "anti-sense"
        return render_template('crispr_response.html', result=result, ssodn_strand=ssodn_strand)
    else:
        CrisprPlanner.change_index_to_match_sense(result)


        buffer = {
            "worm_gene_name": worm_gene_name,
            "nt_seq": nt_seq,
            "site": site,
            "given_isoform": given_isoform,
            "from_aa": from_aa,
            "to_aa": to_aa,
            "end" : end,
            "favourite_enzymes": favourite_enzymes,
            "max_results": max_results,
        }

        return render_template('choose_crrna.html', result=result)


@app.route('/crispr_insertion', methods=['POST'])
def return_crispr_plan_insertion():
    global buffer

    try:
        crrna_done = request.form['crrna_done']
        if crrna_done == 'true':
            print("kkkk")
            worm_gene_name = buffer['worm_gene_name']
            nt_seq = buffer['nt_seq']
            site = buffer['site']
            given_isoform = buffer['given_isoform']
            from_aa = buffer['from_aa']
            to_aa = buffer['to_aa']
            end = buffer['end']
            favourite_enzymes = buffer['favourite_enzymes']
            max_results = buffer['max_results']
            reverse_linker = buffer['reverse_linker']
            forward_linker = buffer['forward_linker']
            primer_reverse = buffer['primer_reverse']
            primer_forward = buffer['primer_forward']
            insertion = buffer['insertion']
            crrna = request.form['crrna']
            crrna_strand = int(request.form.get('crrna_strand'))
            end = buffer['end']

        else:
            print("ddddddddd")
            worm_gene_name = request.form['name']
            nt_seq = None
            end = False
            if request.form['site'] == "start" or request.form['site'] == "START":
                site = 1
            elif request.form['site'] == "":
                return render_template('failure_response.html', query=worm_gene_name, error="The site field is empty.")
            elif request.form['site'] == "end" or request.form['site'] == "END":
                end = True
                site = len(BioPython().get_aa_seq_by_c_elegans_gene_name(worm_gene_name))
                print("new line is in heree" , site)
                print(BioPython().get_aa_seq_by_c_elegans_gene_name(worm_gene_name))
            else:
                site = int(request.form['site'])
            given_isoform = request.form['given_isoform']
            #from_aa = AminoAcid[(aa_names_dict[BioPython().get_aa_seq_by_c_elegans_gene_name(worm_gene_name)[site - 1]])]                  ##junk value - just for the code to work
            from_aa = AminoAcid.ALANINE
            to_aa = AminoAcid.ALANINE                            ##junk value - just for the code to work
            favourite_enzymes = None
            max_results = int(request.form.get('max_results'))
            crrna = request.form['crrna']
            crrna_strand = int(request.form.get('crrna_strand'))
            reverse_linker = request.form['reverse_linker']
            forward_linker = request.form['forward_linker']
            primer_reverse = request.form['reverse_primer']
            primer_forward = request.form['forward_primer']
            insertion = request.form['insertion']


        print("CRISPR Request:",worm_gene_name,given_isoform, site, from_aa, to_aa, max_results, crrna, crrna_strand)

        result, error = CrisprPlanner(gene_name=worm_gene_name,
                                      aa_mutation_site=site,
                                      end = end,
                                      given_isoform=given_isoform,
                                      sense_strand=nt_seq,
                                      favourite_enzymes_names=favourite_enzymes,
                                      type="insertion",
                                      max_results=max_results).plan_my_crispr_insertion(from_aa=from_aa,
                                                                              to_aa=to_aa,
                                                                              crrna=crrna,
                                                                              crrna_strand=crrna_strand,
                                                                            reverse_linker = reverse_linker,
                                                                            forward_linker = forward_linker,
                                                                            primer_reverse = primer_reverse,
                                                                          primer_forward = primer_forward,
                                                                            insertion = insertion)
        if not result or error:
            print(error)
            return render_template('failure_response.html', query=worm_gene_name, error=error)
    except Exception as e:
        error = "Something went wrong: " + str(e)
        print(error)
        return render_template('failure_response.html', query=worm_gene_name, error=error)
    # executor.increment_point_mutation_index(result)

    if crrna:
        print("Result:", result)
        CrisprPlanner.change_index_to_match_sense(result)
        print("Formatted result:", result)
        ssodn_strand = "sense" if result.ssODN_strand > 0 else "anti-sense"
        return render_template('crispr_insertion_response.html', result=result, ssodn_strand=ssodn_strand)

    else:
        CrisprPlanner.change_index_to_match_sense(result)


        buffer = {
            "worm_gene_name": worm_gene_name,
            "nt_seq": nt_seq,
            "site": site,
            "given_isoform": given_isoform,
            "from_aa": from_aa,
            "to_aa": to_aa,
            "favourite_enzymes": favourite_enzymes,
            "max_results": max_results,
            "reverse_linker": reverse_linker,
            "forward_linker": forward_linker,
            "primer_reverse": primer_reverse,
            "primer_forward": primer_forward,
            "insertion": insertion,
            "end": end
        }

        return render_template('choose_crrna_2.html', result=result)


@app.route("/faq")
def faq():
    return render_template("faq.html")


@app.route("/protocol")
def protocol():
    return render_template("Our-injection-mix-and-other-tips.html")


@app.route("/references")
def references():
    return render_template("references.html")


if __name__ == "__main__":
    app.run(host="127.0.0.1", debug=True, port=5005)
