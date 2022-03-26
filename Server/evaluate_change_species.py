from unittest.mock import Mock
from Server.main import compute_network
from Server.evaluation_2 import roc_curve, mcc

main_species_id = "9606"
#secuences_file = 'data/secuences/633813_example.protein.sequences.v11.5.fa'
other_species = [('3702', '59689'), ('224325', '190192')]# (9597, 9593), (4920, 4932), (309807, 518766)]
#aux_species_id_1 = 518766
#aux_species_id_2 = 309807
#for aux_species_id_1, aux_species_id_2 in other_species:
#    print('starting', aux_species_id_1, aux_species_id_2)
#    compute_network(main_species_id, aux_species_id_1, aux_species_id_2, secuences_file, None, None, bar=Mock(), bar_text=Mock())

for aux_species_id_1, aux_species_id_2 in other_species:
    r = roc_curve('data/edges/9606.protein.links.v11.5.txt', f'data/nets/Net_9606_{aux_species_id_1}.csv')
    mcc_values = mcc(r)
    print(aux_species_id_2, max(mcc_values))

# compute_network(main_species_id, aux_species_id_1, aux_species_id_2, secuences_file, None, None, bar=Mock(), bar_text=Mock())