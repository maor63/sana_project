import networkx as nx
import random
import os


def write_leda_to_file(G, graph_name='syeast-05'):
    header = """LEDA.GRAPH\nstring\nlong\n-2\n"""
    nodes = G.nodes()
    node_index = {n: i + 1 for i, n in enumerate(nodes)}
    edges = G.edges()
    path = 'output_graphs/{0}'.format(graph_name)
    if not os.path.exists(path):
        os.makedirs(path)

    output_f = open(os.path.join(path, '/{0}.gw'.format(graph_name)), 'wb')
    output_f.write(header)
    rows = []
    rows.append('{0}'.format(len(nodes)))
    for node in nodes:
        rows.append('|{' + str(node) + '}|')
    rows.append('{0}'.format(len(edges)))
    for v, u in edges:
        m = node_index[v]
        j = node_index[u]
        rows.append('{0} {1} 0 '.format(m, j) + '|{0}|')

    output_f.write('\n'.join(rows))
    output_f.write('\n')
    output_f.close()


# write_leda_to_file(nx.complement(nx.read_leda('networks/syeast0/syeast0.gw')), 'syeast0_complement')
# write_leda_to_file(nx.complement(nx.read_leda('networks/syeast05/syeast05.gw')), 'syeast05_complement')
# write_leda_to_file(nx.complement(nx.read_leda('networks/syeast25/syeast25.gw')), 'syeast25_complement')
# write_leda_to_file(nx.complement(nx.read_leda('networks/yeast/yeast.gw')), 'yeast_complement')
# G_origin = nx.read_leda('networks/syeast0/syeast0.gw')
# edges = G_origin.edges()
# edges_count = len(edges)
# assert isinstance(G, nx.Graph)
# for i in xrange(1, 14):
#     edges = list(nx.non_edges(G))
#     # G.remove_edges_from(random.sample(edges, int(edges_count * 0.05)))
#     G.add_edges_from(random.sample(edges, int(edges_count * 0.05)))
#     dif = str(5 * i + 25)
#     if len(dif) == 1:
#         dif = '0' + dif
#     write_leda_to_file(G, 'syeast{0}'.format(dif))

def add_random_edges(G, edge_count):
    nodes = G.nodes
    added = 0
    while added < edge_count:
        v, u = random.sample(nodes, 2)
        if not G.has_edge(v, u):
            G.add_edge(v, u)
            added += 1
    return G


def add_random_edges_2_graph(G1, G2, edge_count1, edge_count2):
    nodes = G1.nodes
    added = 0
    while added < edge_count1:
        v, u = random.sample(nodes, 2)
        if not G1.has_edge(v, u) and not G2.has_edge(v, u):
            G1.add_edge(v, u)
            G2.add_edge(v, u)
            added += 1
    add_random_edges(G2, edge_count2)
    return G1


def delete_random_edges(G, edge_count):
    edge_list = list(G.edges)
    edge_to_remove = random.sample(edge_list, edge_count)
    G.remove_edges_from(edge_to_remove)
    return G


def delete_random_edges_2_graphs(G1, G2, edge_count1, delta):
    edge_list1 = list(G1.edges)
    edge_to_remove = random.sample(edge_list1, edge_count1)
    G1.remove_edges_from(edge_to_remove)
    G2.remove_edges_from(edge_to_remove)

    edge_list2 = list(G2.edges)
    edge_to_remove = random.sample(edge_list2, delta)
    G2.remove_edges_from(edge_to_remove)


# # generate small graphs
# n = 20
# m = 95
# G1 = nx.gnm_random_graph(n, m)
#
# G1_1 = nx.double_edge_swap(G1.copy(), 5)
# G1_2 = add_random_edges(G1.copy(), 15)
# G1_3 = delete_random_edges(G1.copy(), 15)
#
# G1 = nx.relabel_nodes(G1, {v: 'O{}'.format(v) for v in G1.nodes})
# G1_1 = nx.relabel_nodes(G1_1, {v: 'Rw{}'.format(v) for v in G1_1.nodes})
# G1_2 = nx.relabel_nodes(G1_2, {v: 'AddRw{}'.format(v) for v in G1_2.nodes})
# G1_3 = nx.relabel_nodes(G1_3, {v: 'DelRw{}'.format(v) for v in G1_3.nodes})
# G2 = nx.compose_all([G1_1, G1_2, G1_3])
#
# write_leda_to_file(G1, 'ER20')
# write_leda_to_file(G2, 'ER60_compose')


# input_graph_name = 'syeast0_del_0.15n_fixed_0.75_density'
# G1 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(input_graph_name))
#
# input_graph_name = 'syeast0_fixed_0.75_density'
# G2 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(input_graph_name))
# assert isinstance(G1, nx.Graph)
# base_g1 = 3642
# base_g2 = 5035
#
# precent = 15
# add1 = 3642 * precent
# add2 = 5035 * precent - add1
#
# add_random_edges_2_graph(G1, G2, add1, add2)
# write_leda_to_file(G1, 'syeast0_del_0.15n_fixed_0.90_density')
# write_leda_to_file(G2, 'syeast0_fixed_0.90_density')
# print 'G1 ({},{}) density: {}'.format(len(G1.nodes), len(G1.edges), nx.density(G1))
# print 'G2 ({},{}) density: {}'.format(len(G2.nodes), len(G2.edges), nx.density(G2))


# G1 = nx.complete_graph(854)
# G2 = nx.complete_graph(1004)
# input_graph_name = 'syeast0_del_0.15n_fixed_0.98_density'

# G1 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(input_graph_name))

# input_graph_name = 'syeast0_fixed_0.98_density'
# G2 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(input_graph_name))
# remove1 = 3642 * precent
# remove2 = 5035 * precent - remove1
# delete_random_edges_2_graphs(G1, G2, remove1, remove2)

seed = 100
n = 1000
densety_precent = 1
m = 5 * densety_precent
del_precent = 0.1
# G2 = nx.barabasi_albert_graph(n, m, seed)


input_graph_name = 'syeast0'
G2 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(input_graph_name))
G1_dense = nx.Graph(G2)

nodes = list(G1_dense.nodes)
# nodes_to_delete = int(len(nodes) * del_precent)
nodes_to_delete = int(50)
G1_dense.remove_nodes_from(random.sample(nodes, nodes_to_delete))
# G1_dense.remove_nodes_from(random.sample(nodes, 10))


# total_posiple_edges = (len(G1_dense) * (len(G1_dense) - 1)) / 2.0
#
# G1_equal = nx.Graph(G1_dense)
# delta = nx.density(G1_dense) - nx.density(G2)
# delete_random_edges(G1_equal, int(total_posiple_edges * delta))
#
# G1_sparse = nx.Graph(G1_equal)
# delete_random_edges(G1_sparse, int(total_posiple_edges * delta))

# write_leda_to_file(G1_dense, 'ba_n{}_m{}_del_{}nodes_dense'.format(n, m, del_precent))
# write_leda_to_file(G1_equal, 'ba_n{}_m{}_del_{}nodes_equal'.format(n, m, del_precent))
# write_leda_to_file(G1_sparse, 'ba_n{}_m{}_del_{}nodes_sparse'.format(n, m, del_precent))
# write_leda_to_file(G2, 'ba_n{}_m{}'.format(n, m))

output_graph_name = 'ba_n{}_m{}'.format(n, m)
write_leda_to_file(G1_dense, '{}_del{}'.format(output_graph_name, nodes_to_delete))
write_leda_to_file(G2, '{}'.format(output_graph_name))

print 'G1_dense ({},{}) density: {}'.format(len(G1_dense.nodes), len(G1_dense.edges), nx.density(G1_dense))
# print 'G1_equal ({},{}) density: {}'.format(len(G1_equal.nodes), len(G1_equal.edges), nx.density(G1_equal))
# print 'G1_sparse ({},{}) density: {}'.format(len(G1_sparse.nodes), len(G1_sparse.edges), nx.density(G1_sparse))
print 'G2 ({},{}) density: {}'.format(len(G2.nodes), len(G2.edges), nx.density(G2))

# write_leda_to_file(nx.complement(G2), '%s_complement' % output_graph_name)
# write_leda_to_file(nx.complement(G1_dense), '{}_del{}_complement'.format(output_graph_name, nodes_to_delete))
