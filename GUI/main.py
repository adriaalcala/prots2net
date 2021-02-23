from functools import partial
import tkinter as tk
from tkinter import Button, Entry, filedialog, Label, messagebox, Tk, ttk
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

# Globals
logger = _logging.getLogger("tkinter_.py")

# Constants

selected_input = {
    'secuences_file': '/Users/adria/projects/UIB/proteomica/S.rub_M8.faa',
    'info_file': None
}
finished = True

def init() -> Tk:
    window = Tk()
    
    window.title("Prots2Net app")
    window.geometry('1080x1080')
    window.resizable(width=True, height=True)

    tab_control = ttk.Notebook(window)
    tab1 = ttk.Frame(tab_control)
    tab2 = ttk.Frame(tab_control)
    tab3 = ttk.Frame(tab_control)
    tab_control.add(tab1, text='Select data')
    tab_control.add(tab2, text='Visualize Network')
    tab_control.add(tab3, text='Node info')
    
    tree = ttk.Treeview(tab3, columns=["accession", "attribute"])
    tree.heading("#0",text="Protein idx",anchor=tk.W)
    tree.heading(0, text="Protein accession",anchor=tk.W)
    tree.heading(1, text="Protein attribute",anchor=tk.W)
    tree.pack(side=tk.TOP,fill=tk.X)
    
    sec_btn = Button(tab1,text='Upload sequences csv', anchor='w', command=partial(upload_file, 'secuences_file'))
    sec_btn.grid(column=0,row=0, sticky='W')
    info_btn = Button(tab1,text='Upload info csv', anchor='w', command=partial(upload_file, 'info_file'))
    info_btn.grid(column=0,row=1, sticky='W')
    lbl = Label(tab1, text="Attribute to split proteins",  anchor='w')
    lbl.grid(column=1, row=1, sticky='W')
    attribute_txt = Entry(tab1, width=10)
    attribute_txt.grid(column=2, row=1, sticky='W')
    
    lbl = Label(tab1, anchor='w', text="Main specie id")
    lbl.grid(column=0, row=2, sticky='W')
    sec_txt = Entry(tab1, width=10)
    sec_txt.grid(column=1, row=2, sticky='W')
    sec_txt = Entry(tab1, width=10)
    sec_txt.grid(column=1, row=2, sticky='W')
    
    
    lbl = Label(tab1, anchor='w', text="Similar specie id")
    lbl.grid(column=0, row=3, sticky='W')
    sec_aux_1_txt = Entry(tab1, width=10)
    sec_aux_1_txt.grid(column=1, row=3, sticky='W')
    lbl = Label(tab1, anchor='w', text="Similar specie id")
    lbl.grid(column=2, row=3, sticky='W')
    sec_aux_2_txt = Entry(tab1, width=10)
    sec_aux_2_txt.grid(column=3, row=3, sticky='W')
    bar = ttk.Progressbar(tab1, length=100)
    bar['value'] = 0
    bar.grid(column=1, row=6)
    bar_text = Label(tab1, anchor="s", text="")
    bar_text.grid(column=1, row=7)

    graph_btn = Button(tab2,text='Plot graph', anchor='w', command=partial(draw_graph, tab2, tree))
    graph_btn.pack()
    

    tab_control.pack(expand=2, fill='both')
    btn2 = Button(tab1 ,text='Compute Network', anchor='w', command=partial(clicked, sec_txt, sec_aux_1_txt, sec_aux_2_txt, attribute_txt, bar, bar_text))
    btn2.grid(column=0,row=4, sticky='W')
    return window


def run(window: Tk) -> None:
    window.mainloop()

def upload_file(type_file: str):
    selected_input[type_file] = filedialog.askopenfilename()

def clicked(sec_txt, sec_aux_1_txt, sec_aux_2_txt, attribute_txt, bar, bar_text):
    selected_input['main_species_id'] = 'M7_borja' #sec_txt.get()
    selected_input['aux_species_id_1'] = '518766' # sec_aux_1_txt.get()
    selected_input['aux_species_id_2'] =  '309807' #sec_aux_2_txt.get()
    selected_input['attribute_filter'] = 'day/night' #attribute_txt.get()
    res = messagebox.askyesno('Message title', selected_input)
    print(res)
    if res:

        Thread(target=partial(launch_compute_network, bar=bar, bar_text=bar_text, **selected_input)).start()


def launch_compute_network(bar=None, bar_text=None, main_species_id=None, aux_species_id_1=None, aux_species_id_2=None, secuences_file=None, info_file=None, attribute_filter=None):
    _, net_file_name = compute_network(main_species_id, aux_species_id_1, aux_species_id_2, secuences_file, info_file,  attribute_filter, bar=bar, bar_text=bar_text)
    bar_text.configure(text=f"The network is computed and you can find in {net_file_name}")
    bar_text.text=net_file_name

def create_image(file_name, image_file, tree):
    with open(file_name) as f:
        reader = csv.DictReader(f)
        edges = list(reader)
    edge_list = [(int(r['protein_1_idx']), int(r['protein_2_idx']), float(r['confidence'])) for r in edges]
    visited_nodes = set()
    nodes = list()
    # tree.delete(tree.get_children()[1:])
    for edge in edges:
        if int(edge['protein_1_idx']) not in visited_nodes:
            nodes.append([int(edge['protein_1_idx']), (edge['protein_1'], edge['attribute_p1'])])
            visited_nodes.add(int(edge['protein_1_idx']))
        if int(edge['protein_2_idx']) not in visited_nodes:
            nodes.append([int(edge['protein_2_idx']), (edge['protein_2'], edge['attribute_p2'])])
            visited_nodes.add(int(edge['protein_2_idx']))
    
    for node, values in sorted(nodes, key=lambda x: x[0]):
        tree.insert("", "end",text=node, values=values)
    tree.pack(side=tk.TOP,fill=tk.X)
    G = nx.Graph()
    G.add_weighted_edges_from(edge_list)
    pos = nx.spring_layout(G, iterations=100)
    edge_list_dict = {}
    for attribute in eval(edges[0]['attributes']):
        edge_list_dict[attribute] = [(int(r['protein_1_idx']), int(r['protein_2_idx']), float(r['confidence'])) for r in edges if r['attribute_p1'] == r['attribute_p1'] == attribute]
    figure = plt.figure(figsize=(20,20), dpi=100)

    i = 1
    aux = figure.add_subplot(220 + i)
    i += 1
    aux.set_title("ALL", fontsize=32)
    nx.draw(G, pos, node_size=5, with_labels=True)
    aux.axis('on')
    subplots = [aux]
    for attribute in eval(edges[0]['attributes']):
        subplots.append(figure.add_subplot(220 + i))
        G_aux = nx.Graph()
        G_aux.add_weighted_edges_from(edge_list_dict[attribute])
        
        i += 1
        nx.draw(G_aux, pos, node_size=5, with_labels=True)
        subplots[-1].set_title(attribute, fontsize=32)
        subplots[-1].axis('on')
    plt.savefig(image_file, dpi=100)


def draw_graph(tab2, tree):
    global finished
    file_name = filedialog.askopenfilename()
    image_file = "data/graphs/Net_1234_518766.png"
    create_image(file_name, image_file, tree)
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
