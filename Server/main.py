import csv
import json
import os
import random
import requests
import subprocess

from joblib import dump, load
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO

from logger import logger
from Server.constants import DEFAULT_BLAST_METHOD, INTERACTION_THRESHOLD, STRING_API, STRING_URL, STRING_VERSION


def prepare_paths():
    path_used = ['blast', 'edges', 'graphs', 'nets', 'sequences', 'models']
    for path in path_used:
        if not os.path.exists(f"data/{path}"):
            logger.debug(f"Creating path data/{path}")
            os.makedirs(f"data/{path}")


def create_sequences_file(input_file, species1):
    if input_file is None:
        get_sequences(species1)
    else:
        bash_command = f'cp {input_file} data/sequences/{species1}.faa'
        subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)


def get_file(url: str, output: str, destiny: str, compressed=True):
    output_curl = output if not compressed else f"{output}.gz"
    bash_command = f"curl {url} --output {output_curl}"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    if compressed:
        bash_command = f"gunzip {output_curl}"
        process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
        output, _ = process.communicate()
    bash_command = f"mv {output} {destiny}"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()


def get_sequences(species: str, use_cache=True):
    if os.path.isfile(f'data/sequences/{species}.protein.sequences.{STRING_VERSION}.fa') and use_cache:
        return None
    url = f"{STRING_URL}/protein.sequences.{STRING_VERSION}/{species}.protein.sequences.{STRING_VERSION}.fa.gz"
    get_file(url, "{species}.fa", f"data/sequences/{species}.protein.sequences.{STRING_VERSION}.fa")


def get_edges(species: str, use_cache=True):
    if not os.path.isfile(f'data/edges/{species}.protein.links.{STRING_VERSION}.txt') and use_cache:
        url = f"{STRING_URL}/protein.links.{STRING_VERSION}/{species}.protein.links.{STRING_VERSION}.txt.gz"
        get_file(url, f"{species}_links.txt", f"data/edges/{species}.protein.links.{STRING_VERSION}.txt")
    with open(f'data/edges/{species}.protein.links.{STRING_VERSION}.txt') as f:
        reader = csv.reader(f, delimiter=' ')
        rows = list(reader)
    return {f"{e[0]}_{e[1]}": float(e[2]) for e in rows[1:]}


def build_dataset(species1, species2):
    blast_matrix = compute_blast(species1, species2)
    prot_bit_score = {r['protein_1']: float(r['bit_score']) for r in blast_matrix}
    prot_prot = {r['protein_1']: r['protein_2'] for r in blast_matrix}
    edges = get_edges(species1)
    edges_2 = get_edges(species2)
    prots = list(prot_bit_score.keys())
    training_data = []
    result = []
    train_prots = []
    for p11 in prots:
        p21 = prot_prot[p11]
        for p12 in prots:
            p22 = prot_prot[p12]
            if random.random() > 0.9:
                edge_2_value = edges_2.get(f"{p21}_{p22}") or edges_2.get(f"{p22}_{p21}", 0)
                training_data.append([prot_bit_score[p11], prot_bit_score[p12], edge_2_value])
                edge_1_value = edges.get(f"{p11}_{p12}") or edges.get(f"{p12}_{p11}", 0)
                result.append(1 if edge_1_value > INTERACTION_THRESHOLD else 0)
                train_prots.append([p11, p12])

    return training_data, result


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def compute_blast(species1, species2, write_result=True, method=DEFAULT_BLAST_METHOD):
    method_dispatcher = {
        'STRING': compute_blast_string,
        'LOCAL': compute_blast_local,
        'SPLIT': compute_blast_split
    }
    if method not in method_dispatcher:
        logger.error(f"Method {method} not available, using default method {DEFAULT_BLAST_METHOD}")
        method = DEFAULT_BLAST_METHOD

    return method_dispatcher[method](species1, species2, write_result=write_result)


def compute_blast_string(species1, species2, write_result=True):
    if os.path.isfile(f'data/blast/result_{species1}_{species2}.csv'):
        with open(f'data/blast/result_{species1}_{species2}.csv') as f:
            reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            return list(reader)
    fasta_sequences = SeqIO.parse(open(f"data/sequences/{species1}.protein.sequences.{STRING_VERSION}.fa"), 'fasta')
    names = [f.id for f in fasta_sequences]
    triads = []
    species1 = species1.split('_')[0]
    for chunk in chunks(names, 100):
        identifiers = '%0d'.join(chunk)
        url = f'{STRING_API}/homology_best?identifiers={identifiers}&species={species1}&species_b={species2}'
        response = requests.get(url)
        result = response.json()
        triads.extend([
            {'protein_1': r['stringId_A'], 'protein_2': r['stringId_B'], 'bit_score': r['bitscore']}
            for r in result
        ])

    if write_result:
        with open(f'data/blast/result_{species1}_{species2}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            writer.writerows(triads)
    return triads


def compute_blast_aux(query: str, subject: str):
    bash_command = f'psiblast -query {query} -subject {subject} -outfmt 15'
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    output_dict = json.loads(output)
    return [get_most_similar_sequence(b['report']['results']['bl2seq'][0]) for b in output_dict['BlastOutput2']]


def compute_blast_local(species1, species2, write_result=True):
    if os.path.isfile(f'data/blast/result_{species1}_{species2}.csv'):
        with open(f'data/blast/result_{species1}_{species2}.csv') as f:
            reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            return list(reader)
    sequences1 = f"data/sequences/{species1}.protein.sequences.{STRING_VERSION}.fa"
    sequences2 = f"data/sequences/{species2}.protein.sequences.{STRING_VERSION}.fa"
    triads = compute_blast_aux(sequences1, sequences2)

    if write_result:
        with open(f'data/blast/result_{species1}_{species2}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            writer.writerows(triads)
    return triads


def compute_blast_split(species1, species2, write_result=True):
    if os.path.isfile(f'data/blast/result_{species1}_{species2}.csv'):
        with open(f'data/blast/result_{species1}_{species2}.csv') as f:
            reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            return list(reader)
    all_triads = []
    sequences2 = f"data/sequences/{species2}.faa"
    for sequence_file in os.listdir(f'data/sequences/{species1}'):
        if os.path.isfile(f'data/blast/{species1}/result_{sequence_file.split(".")[0]}_{species2}.csv'):
            with open(f'data/blast/{species1}/result_{sequence_file.split(".")[0]}_{species2}.csv') as f:
                reader = csv.DictReader(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
                all_triads.extend(list(reader))
            continue
        sequences1 = f"data/sequences/{species1}/{sequence_file}"
        triads = compute_blast_aux(sequences1, sequences2)
        all_triads.extend(triads)
        if write_result:
            with open(f'data/blast/{species1}/result_{sequence_file.split(".")[0]}_{species2}.csv', 'w') as f:
                writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
                writer.writerows(triads)
    if write_result:
        with open(f'data/blast/result_{species1}_{species2}.csv', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_2', 'bit_score'])
            writer.writerows(all_triads)
    return all_triads


def get_most_similar_sequence(result):
    protein1 = result['query_title']
    hits = result.get('hits', [{}])[0]
    protein2 = hits.get('description', [{}])[0].get('title')
    bit_score = hits.get('hsps', [{}])[0].get('bit_score', 0)
    return {'protein_1': protein1, 'protein_2': protein2, 'bit_score': bit_score}


def train_model(species1, species2, write_result=True):
    model_path = f'data/models/model_{species1}_{species2}.joblib'
    scaler_path = f'data/models/scaler_{species1}_{species2}.joblib'
    if os.path.isfile(model_path) and os.path.isfile(scaler_path):
        return load(model_path), load(scaler_path)

    training_data, result = build_dataset(species1, species2)
    scaler = StandardScaler()
    scaler.fit(training_data)
    train_set = scaler.transform(training_data)
    clf = MLPClassifier(solver='adam', alpha=1e-6, max_iter=30000, hidden_layer_sizes=(5, 2), random_state=1, tol=1e-8)
    clf.fit(train_set, result).score(train_set, result)
    if write_result:
        dump(clf, model_path)
        dump(scaler, scaler_path)
    return clf, scaler


def set_bar_text(bar_text, text):
    if bar_text:
        bar_text.configure(text=text)
        bar_text.text = text
    else:
        logger.info(text)


def set_bar_step(bar, step, old_step):
    new_step = step + old_step
    if bar:
        bar.step(step)
    else:
        logger.info(new_step)
    return new_step


def compute_network(main_species_id, aux_species_id_1, aux_species_id_2, sequences_file, bar=None, bar_text=None):
    prepare_paths()
    if os.path.isfile(f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv'):
        return None, f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv'
    bar_step = 0
    bar_step = set_bar_step(bar, 0.1, bar_step)
    set_bar_text(bar_text, f"Getting sequences for {aux_species_id_1}")
    get_sequences(aux_species_id_1)
    bar_step = set_bar_step(bar, 5, bar_step)
    set_bar_text(bar_text, f"Getting sequences for {aux_species_id_2}")
    get_sequences(aux_species_id_2)
    bar_step = set_bar_step(bar, 5, bar_step)

    set_bar_text(bar_text, f"Training model")
    clf, scaler = train_model(aux_species_id_2, aux_species_id_1)
    bar_step = set_bar_step(bar, 25, bar_step)

    create_sequences_file(sequences_file, main_species_id)
    set_bar_text(bar_text, f"Computing blast between {main_species_id} and {aux_species_id_1}")
    blast_matrix = compute_blast(main_species_id, aux_species_id_1)
    bar_step = set_bar_step(bar, 25, bar_step)
    set_bar_text(bar_text, f"Getting {aux_species_id_1} edges")
    edges_output = get_edges(aux_species_id_1)
    prot_bitscore = {r['protein_1']: r['bit_score'] for r in blast_matrix}
    prot_prot = {r['protein_1']: r['protein_2'] for r in blast_matrix}
    prots = list(prot_prot.keys())
    prots_idx = {prot: idx for idx, prot in enumerate(prots)}
    new_edges = []
    set_bar_text(bar_text, f"Predicting edges")
    total_prots = len(prots)
    for i in range(total_prots):
        p11 = prots[i]
        bar_step = set_bar_step(bar, 40/total_prots, bar_step)
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
        fieldnames = ['protein_1', 'protein_1_idx', 'protein_2', 'protein_2_idx', 'confidence']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(new_edges)

    return new_edges, f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv'
