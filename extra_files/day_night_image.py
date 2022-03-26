import networkx as nx
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
import csv

file_name = '../data/nets/Net_sal_518766.csv'
image_file = 'data/graphs/Net_sal_518766'
with open(file_name) as f:
    reader = csv.DictReader(f)
    edges = list(reader)
edge_list = [(r['protein_1'], r['protein_2'], 1.0) for r in edges]

G = nx.Graph()
G.add_weighted_edges_from(edge_list)
pos = nx.spring_layout(G, iterations=100)
edge_list_dict = {}
for attribute in eval(edges[0]['attributes']):
    edge_list_dict[attribute] = [(r['protein_1'], r['protein_2'], 1.0) for r in edges if r['attribute_p1'] == r['attribute_p1'] == attribute]
    with open(f'data/nets/Net_sal_518766_{attribute21}.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(edge_list_dict[attribute])

figure = plt.figure(figsize=(100,100), dpi=100)
G_aux = nx.Graph()
G_aux.add_weighted_edges_from(edge_list)
nx.draw(G_aux, pos, node_size=5, with_labels=True)
plt.savefig(f"{image_file}_ALL.png", dpi=100)

for attribute in eval(edges[0]['attributes']):
    figure = plt.figure(figsize=(100,100), dpi=100)
    G_aux = nx.Graph()
    G_aux.add_weighted_edges_from(edge_list_dict[attribute])
    nx.draw(G_aux, pos, node_size=5, with_labels=True)
    plt.savefig(f"{image_file}_{attribute}.png", dpi=100)
