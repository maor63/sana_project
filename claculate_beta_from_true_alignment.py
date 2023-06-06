import networkx as nx
import pandas as pd

g1_name = 'syeast0_delete_0.15nodes'
g2_name = 'syeast0_rewire_37edges'
G1 = nx.read_leda('networks/{0}/{0}.gw'.format(g1_name))
G2 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(g2_name))
true_alignment = pd.DataFrame.from_csv('y0_y2.align', sep='\t').reset_index()
induced_nodes = {v2 for v1, v2 in true_alignment.to_records(index=False)}

assert isinstance(G1, nx.Graph)
assert isinstance(G2, nx.Graph)

Ea_hat = len(G2.subgraph(induced_nodes).edges)
E1 = len(G1.edges) * 1.0
beta = E1 / Ea_hat
print "beta: {}, beta^2: {}, beta sqrt: {}".format(beta, beta**2, beta**0.5)

