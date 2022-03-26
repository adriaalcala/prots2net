from functools import partial

from ttkthemes import ThemedTk
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from threading import Thread
from Server.main import compute_network
import logging as _logging
import platform
import sys
import networkx as nx
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
import csv
from GUI.visualization import ScrollableImage
from GUI.tkentrycomplete import AutocompleteCombobox

# Globals
logger = _logging.getLogger("tkinter_.py")

# Constants

selected_input = {
    'secuences_file': None,
    'info_file': None
}
finished = True
with open('data/species.v11.5.txt') as f:
    reader = csv.DictReader(f, delimiter='\t')
    species = {r['official_name_NCBI']: r['#taxon_id'] for r in reader}

all_species_names = list(species.keys())

def init() -> ThemedTk:
    window = ThemedTk()
    themes = window.get_themes()
    print(themes)
    window.set_theme('itft1')

    window.title("Prots2Net app")
    window.geometry('1080x1080')
    window.resizable(width=True, height=True)

    tab_control = ttk.Notebook(window)
    tab1 = ttk.Frame(tab_control)
    tab2 = ttk.Frame(tab_control)
    tab3 = ttk.Frame(tab_control)
    tab4 = ttk.Frame(tab_control)
    tab_control.add(tab1, text='Select data')
    tab_control.add(tab2, text='Visualize Network')
    tab_control.add(tab3, text='Node info')
    tab_control.add(tab4, text='Edge info')

    tree_node = ttk.Treeview(tab3, columns=["accession"])
    tree_node.heading("#0", text="Protein idx",anchor=tk.W)
    tree_node.heading(0, text="Protein accession", anchor=tk.W)
    tree_node.pack(side=tk.TOP, fill=tk.X)

    tree_edge = ttk.Treeview(tab4, columns=["idx_1", "idx_2"])
    tree_edge.heading("#0", text="Interaction idx", anchor=tk.W)
    tree_edge.heading(0, text="Protein idx 1", anchor=tk.W)
    tree_edge.heading(1, text="Protein idx 2", anchor=tk.W)
    tree_edge.pack(side=tk.TOP, fill=tk.X)

    sec_btn = ttk.Button(tab1, text='Upload sequences', command=partial(upload_file, 'secuences_file'))
    sec_btn.grid(column=0, row=0, sticky='W')
    # info_btn = Button(tab1,text='Upload info csv', anchor='w', command=partial(upload_file, 'info_file'))
    # info_btn.grid(column=0,row=1, sticky='W')
    # lbl = Label(tab1, text="Attribute to split proteins",  anchor='w')
    # lbl.grid(column=1, row=1, sticky='W')
    # attribute_txt = Entry(tab1, width=10)
    # attribute_txt.grid(column=2, row=1, sticky='W')

    lbl = ttk.Label(tab1, anchor='w', text="Main species id")
    lbl.grid(column=0, row=2, sticky='W')
    sec_txt = ttk.Entry(tab1, width=10)
    sec_txt.grid(column=1, row=2, sticky='W')
    sec_txt = ttk.Entry(tab1, width=10)
    sec_txt.grid(column=1, row=2, sticky='W')


    lbl = ttk.Label(tab1, anchor='w', text="Similar species")
    lbl.grid(column=0, row=3, sticky='W')
    combo_1 = AutocompleteCombobox(tab1)
    combo_1["values"] = all_species_names
    combo_1.set_completion_list(all_species_names)
    combo_1.grid(column=1, row=3, sticky='W')
    #var_1 = tk.StringVar()
    #entry = ttk.Entry(tab1, textvariable=var_1)
    #entry.pack()
    # search_button = Button(tab)
    # sec_aux_1_txt = ttk.Entry(tab1, width=10)
    # sec_aux_1_txt.grid(column=1, row=3, sticky='W')
    lbl = ttk.Label(tab1, anchor='w', text="   Similar species   ")
    lbl.grid(column=2, row=3, sticky='W')
    combo_2 = AutocompleteCombobox(tab1)  #ttk.Combobox(tab1) #, state="readonly")
    combo_2["values"] = all_species_names
    combo_2.set_completion_list(all_species_names)
    combo_2.grid(column=4, row=3, sticky='W')
#    sec_aux_2_txt = ttk.Entry(tab1, width=10)
#    sec_aux_2_txt.grid(column=4, row=3, sticky='W')
    bar = ttk.Progressbar(tab1, length=100)
    bar['value'] = 0
    bar.grid(column=1, row=6)
    bar_text = ttk.Label(tab1, anchor="s", text="")
    bar_text.grid(column=1, row=7)

    graph_btn = ttk.Button(tab2,text='Plot graph', command=partial(draw_graph, tab2, tree_node, tree_edge))
    graph_btn.pack()


    tab_control.pack(expand=2, fill='both')
    btn2 = ttk.Button(tab1 ,text='Compute Network', command=partial(clicked, sec_txt, combo_1, combo_2, bar, bar_text))
    btn2.grid(column=0,row=4, sticky='W')
    return window

def run(window: ThemedTk) -> None:
    window.mainloop()

def upload_file(type_file: str):
    selected_input[type_file] = filedialog.askopenfilename()

def clicked(sec_txt, combo_1, combo_2, bar, bar_text):
    selected_input['main_species_id'] = sec_txt.get()
    # selected_input['aux_species_id_1'] = sec_aux_1_txt.get()
    selected_input['aux_species_id_1'] = species.get(combo_1.get())
    selected_input['aux_species_id_2'] = species.get(combo_2.get())
    res = messagebox.askyesno('Message title', selected_input)
    print(res)
    # if res:

      #  Thread(target=partial(launch_compute_network, bar=bar, bar_text=bar_text, **selected_input)).start()


def launch_compute_network(bar=None, bar_text=None, main_species_id=None, aux_species_id_1=None, aux_species_id_2=None, secuences_file=None, info_file=None, attribute_filter=None):
    _, net_file_name = compute_network(main_species_id, aux_species_id_1, aux_species_id_2, secuences_file, info_file,  attribute_filter, bar=bar, bar_text=bar_text)
    bar_text.configure(text=f"The network is computed and you can find in {net_file_name}")
    bar_text.text=net_file_name

def create_image(file_name, image_file, tree_node, tree_edge):
    print(file_name)
    with open(file_name) as f:
        reader = csv.DictReader(f)
        edges = list(reader)
    edge_list = [(int(r['protein_1_idx']), int(r['protein_2_idx']), float(r['confidence'])) for r in edges]
    visited_nodes = set()
    nodes = list()
    # tree.delete(tree.get_children()[1:])
    edges_tree = list()
    for idx, edge in enumerate(edges):
        if int(edge['protein_1_idx']) not in visited_nodes:
            nodes.append([int(edge['protein_1_idx']), (edge['protein_1'], )])
            visited_nodes.add(int(edge['protein_1_idx']))
        if int(edge['protein_2_idx']) not in visited_nodes:
            nodes.append([int(edge['protein_2_idx']), (edge['protein_2'], )])
            visited_nodes.add(int(edge['protein_2_idx']))
        edges_tree.append([idx, (edge['protein_1_idx'], edge['protein_2_idx'])])
    
    for node, values in sorted(nodes, key=lambda x: x[0]):
        tree_node.insert("", "end", text=node, values=values)
    tree_node.pack(side=tk.TOP, fill=tk.X)
    for edge, values in sorted(edges_tree, key=lambda x: x[0]):
        tree_edge.insert("", "end", text=edge, values=values)
    tree_edge.pack(side=tk.TOP, fill=tk.X)
    G = nx.Graph()
    G.add_weighted_edges_from(edge_list)
    pos = nx.spring_layout(G, iterations=100)
    figure = plt.figure(figsize=(20,20), dpi=100)

    nx.draw(G, pos, node_size=5, with_labels=True)
    plt.savefig(image_file, dpi=100)


def draw_graph(tab2, tree_node, tree_edge):
    global finished
    file_name = filedialog.askopenfilename()
    main_species, species_1 = file_name.split('/')[-1].split('.')[0].split('_')[1:]
    image_file = f"data/graphs/Net_{main_species}_{species_1}.png"
    create_image(file_name, image_file, tree_node, tree_edge)
    if '!scrollableimage' in tab2.children:
        image_window = tab2.children['!scrollableimage']
        image_window.grid_forget()
        image_window.destroy()
    img = ImageTk.PhotoImage(Image.open(image_file))
    image_window = ScrollableImage(tab2, image=img, scrollbarwidth=6, width=1000, height=1000)
    image_window.pack()

if __name__ == '__main__':
    logger.setLevel(_logging.INFO)
    stream_handler = _logging.StreamHandler()
    formatter = _logging.Formatter("[%(filename)s] %(message)s")
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    app = init()
    run(app)
    # # Tk must be initialized before CEF otherwise fatal error (Issue #306)
    # cef.Initialize()
    
    # app.mainloop()
    # cef.Shutdown()
