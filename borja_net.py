from Server.main import *
main_species_id, aux_species_id_1, aux_species_id_2 = 'aaborjaM7', '518766', '309807'
print('start model train')
clf, scaler = train_model(aux_species_id_2, aux_species_id_1)
print('model train finished')
blast_matrix = compute_blast_splitted(main_species_id, aux_species_id_1)
edges_output = get_edges(aux_species_id_1)
prot_bitscore = {r['protein_1']: float(r['bit_score']) * 20 for r in blast_matrix}
prot_prot = {r['protein_1']: r['protein_2'] for r in blast_matrix}
prots = list(prot_prot.keys())
prots_idx = {p:i for i,p in enumerate(prots)}
new_edges = []
for i in range(len(prots)):
    if i%1000 == 0:
        print(f"{i} / {len(prots)}")
    p11 = prots[i]
    p21 = prot_prot[p11]
    if p21 == '':
        continue
    print('Protein found', i)

    p12_list = [prots[j] for j in range(i, len(prots))]
    p22_list = [prot_prot[p12] for p12 in p12_list]
    edge2_value_list = [edges_output.get(f'{p21}_{p22}') or edges_output.get(f'{p22}_{p21}', 0) for p22 in p22_list]
    if set(edge2_value_list) == {0}:
        continue
    print('edge_2 found', set(edge2_value_list))
    v_list = [(prot_bitscore[p11], prot_bitscore[p12], edge2_value) for p12, edge2_value in zip(p12_list, edge2_value_list)]
    predict_list = clf.predict(scaler.transform(v_list))
    if set(predict_list) != {0}:
        new_edges_aux = [{'protein_1': p11, 'protein_1_idx': prots_idx[p11], 'protein_2': p12, 'protein_2_idx': prots_idx[p12]} for p12, p in zip(p12_list, predict_list) if p > 0]
        new_edges.extend(new_edges_aux)
        print('edges found!!!!')
        print(len(new_edges))
with open(f'data/nets/Net_{main_species_id}_{aux_species_id_1}.csv', 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['protein_1', 'protein_1_idx', 'protein_2', 'protein_2_idx', 'confidence', 'attributes', 'attribute_p1', 'attribute_p2'])
        writer.writeheader()
        writer.writerows(new_edges)


# i = 0
# p11 = prots[i]
# p21 = prot_prot[p11]


# p12_list = [prots[j] for j in range(i, len(prots))]
# p22_list = [prot_prot[p12] for p12 in p12_list]
# edge2_value_list = [edges_output.get(f'{p21}_{p22}') or edges_output.get(f'{p22}_{p21}', 0) for p22 in p22_list]

# print('edge_2 found', set(edge2_value_list))
# v_list = [(prot_bitscore[p11], prot_bitscore[p12], edge2_value) for p12, edge2_value in zip(p12_list, edge2_value_list)]
# predict_list = clf.predict(scaler.transform(v_list))
