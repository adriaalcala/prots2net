import csv
import json
import os
import random
import requests
import subprocess
import platform
import tarfile
from time import sleep

from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO
import Server.VisualizeNN as VisNN
import numpy as np


def create_secuences_file(input_file, species1):
    print(input_file, species1)
    if input_file is None:
        get_secuences(species1)
    else:
        bashcommand = f'cp {input_file} data/secuences/{species1}.faa'
        subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)

def read_info(input_file):
    pass


def get_secuences(species: str, use_cache=True):
    if os.path.isfile(f'data/secuences/{species}.protein.sequences.v11.5.fa') and use_cache:
        return None
    url = f"https://stringdb-static.org/download/protein.sequences.v11.5/{species}.protein.sequences.v11.5.fa.gz"
    bashcommand = f"curl {url} --output {species}.fa.gz"
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    print(output)
    bashcommand = f"gunzip {species}.fa.gz"
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    print(output)
    bashcommand = f"mv {species}.fa data/secuences/{species}.protein.sequences.v11.5.fa"
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    print(output)


def download_edges(species: str):
    url = f"https://stringdb-static.org/download/protein.links.v11.5/{species}.protein.links.v11.5.txt.gz"
    bashcommand = f"curl {url} --output {species}_links.txt.gz"
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    print(output)
    bashcommand = f"gunzip {species}_links.txt.gz"
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    print(output)
    bashcommand = f"mv {species}_links.txt data/edges/{species}.protein.links.v11.5.txt"
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    print(output)
    

def get_edges(species: str, use_cache=True):
    if not os.path.isfile(f'data/edges/{species}.protein.links.v11.5.txt') and use_cache:
        download_edges(species)
    with open(f'data/edges/{species}.protein.links.v11.5.txt') as f:
        reader = csv.reader(f, delimiter=' ')
        rows = list(reader)
    return {f"{e[0]}_{e[1]}": float(e[2]) for e in  rows[1:]}


def build_dataset(species1, species2):
    blast_matrix = compute_blast(species1, species2)
    prot_bitscore = {r['protein_1']: float(r['bit_score']) for r in blast_matrix}
    prot_prot = {r['protein_1']: r['protein_2'] for r in blast_matrix}
    edges = get_edges(species1)
    edges_2 = get_edges(species2)
    prots = list(prot_bitscore.keys())
    print('start build dataset')
    training_data = []
    result = []
    train_prots = []
    i = 0
    total = len(prots)
    for p11 in prots:
        i += 1
        if random.random() > 0.8:
            print(i, total)
        p21 = prot_prot[p11]
        for p12 in prots:
            p22 = prot_prot[p12]
            if random.random() > 0.9:
                edge_2_value = edges_2.get(f"{p21}_{p22}") or edges_2.get(f"{p22}_{p21}", 0)
                training_data.append([prot_bitscore[p11], prot_bitscore[p12], edge_2_value])
                edge_1_value = edges.get(f"{p11}_{p12}") or edges.get(f"{p12}_{p11}", 0)
                result.append(1 if edge_1_value > 700 else 0)
                train_prots.append([p11, p12])

    print('end build dataset')

    return training_data, result


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def compute_blast_new(species1, species2, write_result=True):
    print('blast', species1, species2)
    if os.path.isfile(f'data/blast/result_{species1}_{species2}.csv'):
        with open(f'data/blast/result_{species1}_{species2}.csv') as f:
            reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            return list(reader)
    print('start compute blast')
    fasta_sequences = SeqIO.parse(open(f"data/secuences/{species1}.protein.sequences.v11.5.fa"),'fasta')
    names = [f.id for f in fasta_sequences]
    triads = []
    i = 0
    total = len(names) / 100
    species1 = species1.split('_')[0]
    for chunk in chunks(names, 100):
        identifiers = '%0d'.join(chunk)
        url = f'https://string-db.org/api/json/homology_best?identifiers={identifiers}&species={species1}&species_b={species2}'
        #print('url', url)
        response = requests.get(url)
        #print(response)
        result = response.json()
        #print(result)
        triads.extend([{'protein_1': r['stringId_A'], 'protein_2': r['stringId_B'], 'bit_score': r['bitscore']} for r in result])
        i += 1
        print(i, total)
    print(triads[-1])
    if write_result:
        with open(f'data/blast/result_{species1}_{species2}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            writer.writerows(triads)
    return triads


def compute_blast(species1, species2, write_result=True):
    if os.path.isfile(f'data/blast/result_{species1}_{species2}.csv'):
        with open(f'data/blast/result_{species1}_{species2}.csv') as f:
            reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            return list(reader)
    print('start compute blast')
    bashcommand = f'psiblast -query data/secuences/{species1}.protein.sequences.v11.5.fa -subject data/secuences/{species2}.protein.sequences.v11.5.fa -outfmt 15'
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    print('end compute blast')
    output_dict = json.loads(output)
    triads = [get_most_similar_sequence(b['report']['results']['bl2seq'][0]) for b in output_dict['BlastOutput2']]
    if write_result:
        with open(f'data/blast/result_{species1}_{species2}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            writer.writerows(triads)
    return triads

def compute_blast_splitted(species1, species2, write_result=True):
    if os.path.isfile(f'data/blast/result_{species1}_{species2}.csv'):
        with open(f'data/blast/result_{species1}_{species2}.csv') as f:
            reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            return list(reader)
    print('start compute blast')
    all_triads = []
    for secuence_file in os.listdir(f'data/secuences/{species1}'):
        print(secuence_file)
        if os.path.isfile(f'data/blast/{species1}/result_{secuence_file.split(".")[0]}_{species2}.csv'):
            with open(f'data/blast/{species1}/result_{secuence_file.split(".")[0]}_{species2}.csv') as f:
                reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
                all_triads.extend(list(reader))
            continue
        bashcommand = f'psiblast -query data/secuences/{species1}/{secuence_file} -subject data/secuences/{species2}.faa -outfmt 15'
        process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
        output, _ = process.communicate()
        print('end compute blast')
        output_dict = json.loads(output)
        triads = [get_most_similar_sequence(b['report']['results']['bl2seq'][0]) for b in output_dict['BlastOutput2']]
        all_triads.extend(triads)
        if write_result:
            with open(f'data/blast/{species1}/result_{secuence_file.split(".")[0]}_{species2}.csv', 'w') as f:
                writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
                writer.writerows(triads)
    if write_result:
        with open(f'data/blast/result_{species1}_{species2}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            writer.writerows(all_triads)
    return all_triads

def get_most_similar_sequence(result):
    protein1 = result['query_title']
    if result['hits']:
        if result['hits'][0]['description']:
            protein2 = result['hits'][0]['description'][0]['title']
        else:
            protein2 = None
        if result['hits'][0]['hsps']:
            bitscore = result['hits'][0]['hsps'][0]['bit_score']
        else:
            bitscore = 0
    else:
        protein2 = None
        bitscore = 0
    return {'protein_1': protein1, 'protein_2': protein2, 'bit_score': bitscore}



def train_model(species1, species2):
    training_data, result = build_dataset(species1, species2)
    scaler = StandardScaler()
    scaler.fit(training_data)
    X_train = scaler.transform(training_data)
    clf = MLPClassifier(solver='adam', alpha=1e-6, max_iter=30000,hidden_layer_sizes=(5, 2), random_state=1, tol=1e-8)
    print('start fit')
    clf.fit(X_train, result).score(X_train, result)
    print('model trained')
    return clf, scaler


def set_bar_text(bar_text, text):
    bar_text.configure(text=text)
    bar_text.text = text

def get_size(info_file, attribute_filter):
    with open(info_file) as f:
        reader = csv.DictReader(f)
        info = list(reader)
    sizes = {r['ID'] : r[attribute_filter] for r in info}
    attributes = list(set(sizes.values()))
    return attributes, sizes

def compute_network(main_species_id, aux_species_id_1, aux_species_id_2, secuences_file, info_file, attribute_filter, bar=None, bar_text=None):
    if os.path.isfile(f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv'):
        return None, f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv'
    bar.step(0.1)
    set_bar_text(bar_text, f"Getting secuences for {aux_species_id_1}")
    get_secuences(aux_species_id_1)
    bar.step(5)
    set_bar_text(bar_text, f"Getting secuences for {aux_species_id_2}")
    get_secuences(aux_species_id_2)
    bar.step(5)


    set_bar_text(bar_text, f"Training model")
    clf, scaler = train_model(aux_species_id_2, aux_species_id_1)
    bar.step(25)

    create_secuences_file(secuences_file, main_species_id)
    set_bar_text(bar_text, f"Computing blast between {main_species_id} and {aux_species_id_1}")
    blast_matrix = compute_blast(main_species_id, aux_species_id_1)
    bar.step(25)
    set_bar_text(bar_text, f"Getting {aux_species_id_1} edges")
    edges_output = get_edges(aux_species_id_1)
    prot_bitscore = {r['protein_1']: r['bit_score'] for r in blast_matrix}
    prot_prot = {r['protein_1']: r['protein_2'] for r in blast_matrix}
    prots = list(prot_prot.keys())
    prots_idx = {p:i for i,p in enumerate(prots)}
    new_edges = []
    set_bar_text(bar_text, f"Predicting edges")
    for i in range(len(prots)):
        p11 = prots[i]
        print(f"{i} / {len(prots)}")
        bar.step(40/len(prots))
        for j in range(i, len(prots)):
            p12 = prots[j]
            p21 = prot_prot[p11]
            p22 = prot_prot[p12]
            edge2_value = edges_output.get(f'{p21}_{p22}') or edges_output.get(f'{p22}_{p21}', 0)
            v = [prot_bitscore[p11], prot_bitscore[p12], edge2_value]
            predict = clf.predict(scaler.transform([v]))[0]
            if predict > 0:
                new_edges.append(
                    {
                        'protein_1': p11, 'protein_1_idx': prots_idx[p11],
                        'protein_2': p12, 'protein_2_idx': prots_idx[p12],
                        'confidence': clf.predict_proba(scaler.transform([v]))[0][1],
                    }
                )
    with open(f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_1_idx', 'protein_2', 'protein_2_idx', 'confidence'])
            writer.writeheader()
            writer.writerows(new_edges)


    return new_edges, f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv'
