import os
import networkx as nx
from sklearn.metrics import adjusted_mutual_info_score

from scripts_aux import MOD_Read

def write_rows_to_file(rows, output_path):
    try:
        with open(output_path, 'w') as o:
            for row in rows:
                row = [str(elem) for elem in row]
                o.write('\t'.join(row))
                o.write('\n')
            o.close()
    except IOError:
        print('Can\'t write data on file.')

# Read all DCLs files
dcl_dir = "amp-ad/results/"
dcl_files = os.listdir(dcl_dir)
dcl_files = [dcl_dir + dcl_file for dcl_file in dcl_files if dcl_file.count("DCL") > 0]

target_genes = MOD_Read.read_rows_from_tab_file("scripts_aux/genes_diffcoexp_input.txt", header = False)
target_genes = [gene[0] for gene in target_genes]

# Createing the AD-DiffCoexpNet
GFG = nx.Graph()
dcls_fusiform = MOD_Read.read_rows_from_tab_file("amp-ad/results/GEO_ GSE125583_fusiform gyrus_DCLs.txt", header = True)
for row in dcls_fusiform:
    GFG.add_edge(row[1], row[2])

glist = []
tissue_list = []

for dcl_file in dcl_files:
    tissue = dcl_file.split("_")
    tissue = tissue[1] + " :: " + tissue[2]
    G = nx.Graph()
    dcls = MOD_Read.read_rows_from_tab_file(dcl_file, header = True)
    if len(dcls) > 0:
        for row in dcls:
            G.add_edge(row[1], row[2])
        glist.append(G)
        tissue_list.append(tissue)


def jaccard_coeficient_nodes(G1, G2):
    """
    Calculate the Jaccard Index between two graphs regarding their nodes.
    Parameters
    ----------
        G1 : a graph from networkx class nx.Graph().
        G2 : a graph from networkx class nx.Graph().
    """
    try:
        print("Calculating JC index ...")
        if isinstance(G1, nx.Graph) and isinstance(G2, nx.Graph):
            if len(G1.nodes()) > 0 and len(G2.nodes()) > 0:
                G1_nodes = set(G1.nodes())
                G2_nodes = set(G2.nodes())
                N = len(G1_nodes.union(G2_nodes))
                I = len(G1_nodes.intersection(G2_nodes))
                JC = I / N
                return(JC)
            return(0)
        else:
            raise("Verify the class of your graph.")
    except:
        print("Verify the class of your graph.")
    finally:
        print("Done!")

def jaccard_coeficient_edges(G1, G2):
    """
    Calculate the Jaccard Index between two graphs regarding their nodes.
    Parameters
    ----------
        G1 : a graph from networkx class nx.Graph().
        G2 : a graph from networkx class nx.Graph().
    """
    try:
        print("Calculating JC index ...")
        if isinstance(G1, nx.Graph) and isinstance(G2, nx.Graph):
            if len(G1.edges()) > 0 and len(G2.edges()) > 0:
                G1_nodes = set(G1.edges())
                G2_nodes = set(G2.edges())
                N = len(G1_nodes.union(G2_nodes))
                I = len(G1_nodes.intersection(G2_nodes))
                JC = I / N
                return(JC)
            return(0)
        else:
            raise("Verify the class of your graph.")
    except:
        print("Verify the class of your graph.")
    finally:
        print("Done!")


results = []
for i in range(0, len(tissue_list)):
    g1 = glist[i]
    t1 = tissue_list[i]
    for j in range(0, len(tissue_list)):
        g2 = glist[j]
        t2 = tissue_list[j]
        jc1 = jaccard_coeficient_nodes(g1, g2)
        jc2 = jaccard_coeficient_edges(g1, g2)
        results.append([t1, t2, jc1, jc2 ])
write_rows_to_file(results, "results/Jaccard.txt")

all_nodes = []
for i in range(0, len(tissue_list)):
    g = glist[i]
    all_nodes.extend(list(g.nodes()))
all_nodes = sorted(all_nodes)

binary_matrix = []
for i in range(0, len(tissue_list)):
    tissue = [tissue_list[i]]
    for gene in all_nodes:
        g = glist[i]
        if gene in list(GFG.nodes()):
            tissue.append(1)
        else:
            tissue.append(0)
    binary_matrix.append(tissue)

results = []

# sklearn.metrics.adjusted_mutual_info_score

print("Calculating the Adjusted Mutual Info Score")
for i in range(0, len(tissue_list)):
    for j in range(0, len(tissue_list)):
        i_aux = binary_matrix[i]
        j_aux = binary_matrix[j]
        mri = adjusted_mutual_info_score(i_aux[1:], j_aux[1:])
        results.append([i_aux[0], j_aux[0], mri])
write_rows_to_file(results, "results/MRI.txt")
print("Done!")
