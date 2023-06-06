from __future__ import print_function

import collections
from collections import Counter

import matplotlib.pyplot as plt
import timeit
from multiprocessing import Pool, Process
import networkx as nx
import random
import pandas as pd
import os
import numpy as np

from eval_true_alignment import get_node_data


def write_leda_to_file(G, graph_name, path):
    header = """LEDA.GRAPH\nstring\nlong\n-2\n"""
    nodes = G.nodes()
    node_index = {n: i + 1 for i, n in enumerate(nodes)}
    edges = G.edges()
    if not os.path.exists(path):
        os.makedirs(path)

    output_f = open(os.path.join(path, '{0}.gw'.format(graph_name)), 'wb')
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


def get_nc(alignment, true_alignment):
    correct = sum(1 for v in alignment.index if alignment[v] == true_alignment[v])
    return correct * 1.0 / len(alignment)


def update_Ea_change(v, u_old, u_new, G1, G2, A):
    res = 0
    for n in G1.neighbors(v):
        res -= 1 if G2.has_edge(u_old, A[n]) else 0
        res += 1 if G2.has_edge(u_new, A[n]) else 0
    return res


def update_Ea_hat(u_old, u_new, G2, set_free_nodes):
    res = 0
    for n in G2.neighbors(u_old):
        res -= 1 if n not in set_free_nodes else 0
    for n in G2.neighbors(u_new):
        res += 1 if n not in set_free_nodes else 0

    res -= 1 if G2.has_edge(u_old, u_new) else 0
    return res


def add_random_edges(G, k):
    assert isinstance(G, nx.Graph)
    edges_to_add = random.sample(list(nx.non_edges(G)), k)
    G.add_edges_from(edges_to_add)
    return G


def remove_random_edges(G, k):
    assert isinstance(G, nx.Graph)
    edges_to_remove = random.sample(list(G.edges), k)
    G.remove_edges_from(edges_to_remove)
    return G


def get_wattz_graphs():
    n = 1000
    WS = nx.watts_strogatz_graph(n, 18, 0.5)
    WS.graph['name'] = 'ws_1000_18_0.5'
    WS_bal = nx.watts_strogatz_graph(n, 500, 0.5)
    WS_bal.graph['name'] = 'ws_1000_500_0.5'
    WS_dense = nx.watts_strogatz_graph(n, 980, 0.5)
    WS_dense.graph['name'] = 'ws_1000_982_0.5'
    G2s = [WS, WS_bal, WS_dense]
    return G2s


def get_erdos_graphs():
    n = 1000
    ER = nx.gnm_random_graph(n, 8500)
    ER.graph['name'] = 'er_1000_8500'
    ER_bal = nx.gnm_random_graph(n, 250000)
    ER_bal.graph['name'] = 'er_1000_250000'
    ER_dense = nx.gnm_random_graph(n, 491000)
    ER_dense.graph['name'] = 'er_1000_491000'
    G2s = [ER, ER_bal, ER_dense]
    return G2s


def get_barabasi_graphs():
    n = 1000
    BA = nx.barabasi_albert_graph(n, 9)
    BA.graph['name'] = 'ba_1000_9'
    BA_bal = nx.barabasi_albert_graph(n, 500)
    BA_bal.graph['name'] = 'ba_1000_500'
    BA_compliment = nx.complement(nx.barabasi_albert_graph(n, 9))
    BA_compliment.graph['name'] = 'ba_1000_9_compliment'
    G2s = [BA, BA_bal, BA_compliment]
    return G2s


def syeast0_graphs():
    G2_name = 'syeast0'
    Y0 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(G2_name))
    Y0.graph['name'] = G2_name
    Y0_complement = nx.complement(Y0)
    Y0_complement.graph['name'] = G2_name + "_complement"
    write_leda_to_file(Y0_complement, Y0_complement.graph['name'], 'output_alignment/graphs')
    Y0_V = len(Y0.nodes)
    edges_count_to_add = int(Y0_V * (Y0_V - 1) / 4.0 - len(Y0.edges))
    Y0_balanced = add_random_edges(nx.Graph(Y0), edges_count_to_add)
    Y0_balanced.graph['name'] = 'syeast0_balanced'
    write_leda_to_file(Y0_balanced, Y0_balanced.graph['name'], 'output_alignment/graphs')
    G2s = [Y0, Y0_complement, Y0_balanced]
    return G2s


def generate_alignments(G2s, edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, shuffles,
                        output_file_name='alignmnets_y0_com_bal'):
    iter_counter = 1
    for i in range(full_iterations):
        for j, G2 in enumerate(G2s):
            for nodes_to_delete in nodes_to_delete_list:
                for edge_frac_To_modify in edge_frac_To_modify_list:
                    start = timeit.default_timer()
                    G2_name = G2.graph['name']
                    print('full iteration {}/{}, using G2 {}, nodes_to_delete {}, edge_to_delete {} {}/{}'.format(
                        str(i + 1), full_iterations, G2_name, nodes_to_delete, edge_frac_To_modify, iter_counter,
                        str(len(G2s) * len(nodes_to_delete_list) * len(edge_frac_To_modify_list) * full_iterations)))
                    iter_counter += 1
                    G1 = generate_G1(G2, nodes_to_delete, edge_frac_To_modify)
                    G1_name = '{}_del_nodes_{}_modify_{}_iter_{}'.format(G2.graph['name'], nodes_to_delete,
                                                                         edge_frac_To_modify, i)
                    write_leda_to_file(G1, G1_name, 'output_alignment/graphs')

                    generate(G1, G1_name, G2, G2_name, edge_frac_To_modify, i, nodes_to_delete, output_file_name,
                             shuffles, start)


def generate(G1, G1_name, G2, G2_name, edge_frac_To_modify, i, nodes_to_delete, output_file_name, shuffles, start):
    V1 = list(G1.nodes)
    V2 = list(G2.nodes)
    E2 = len(list(G2.edges))
    true_A = pd.Series(data=V1, index=V1)
    E1, Ea, Ea_hat = get_align_edges(G1, G2, true_A)
    Ea_hat_star = Ea_hat
    Ea_star = Ea
    omega = (len(V1) * (len(V1) - 1)) / 2.0
    try:
        nc = get_nc(A, true_A)
        measures = get_measures(E1, Ea, Ea_hat, omega)
        acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures

        row = (
            G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, '', 0, '', 0,
            '', 0,
            nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk, mcc, f01, f033, f3, f10)
        rows = [row]
    except Exception as e:
        rows = []
    original_Ea = Ea
    original_Ea_hat = Ea_hat
    set_free_nodes = set(G2.nodes) - set(G1.nodes)
    for shuffle in range(shuffles):
        print('\r shuffle iteration {}/{}'.format(shuffle, shuffles), end='')
        A = pd.Series(true_A.copy())
        random.shuffle(V1)
        G1_name = '{}_del_nodes_{}_modify_{}_iter_{}_shuffle_{}'.format(G2.graph['name'],
                                                                        nodes_to_delete,
                                                                        edge_frac_To_modify, i,
                                                                        shuffle)
        Ea = original_Ea
        Ea_hat = original_Ea_hat

        for v in V1:
            u = A[v]
            u2 = random.choice(list(set_free_nodes))

            ea_add = update_Ea_change(v, u, u2, G1, G2, A)
            ea_hat_add = update_Ea_hat(u, u2, G2, set_free_nodes)

            Ea += ea_add
            Ea_hat += ea_hat_add

            change(A, v, u2, set_free_nodes)
            try:
                nc = get_nc(A, true_A)
                measuers = get_measures(E1, Ea, Ea_hat, omega)
                acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measuers

                row = (
                    G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, v,
                    G1.degree[v], u,
                    G2.degree[u], u2, G2.degree[u2], nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk,
                    mcc,
                    f01,
                    f033, f3, f10)
                rows.append(row)
            except Exception as e:
                print(e)
    stop = timeit.default_timer()
    print('\r shuffle iteration {}/{} time: {}sec'.format(shuffles, shuffles, stop - start), end='')
    print()
    df = pd.DataFrame(rows,
                      columns=['G1', 'V1', 'E1', 'G2', 'V2', 'E2', 'omega', 'Ea_hat*', 'Ea', 'Ea_hat',
                               'G1 changed node', 'Degree of changed node', 'G2 old node',
                               'Degree of old node', 'G2 new node', 'Degree of new node', 'NC', 'TP',
                               'FP', 'TN', 'FN', 'EC', 'ICS', 'S3', 'F1', 'ACC', 'BM', 'MK', 'MCC',
                               'F0.1', 'F0.33', 'F3', 'F10'])
    output_path = 'output_alignment/%s.csv' % output_file_name
    if os.path.isfile(output_path):
        df.to_csv(output_path, mode='a', header=None)
    else:
        df.to_csv(output_path)


def change(A, v, u2, set_free_nodes):
    u = A[v]
    A[v] = u2
    set_free_nodes.remove(u2)
    set_free_nodes.add(u)


def get_measures(E1, Ea, Ea_hat, omega):
    ec, f1, ics, s3, bm, mk, mcc = compute_measures(E1, Ea, Ea_hat, omega)
    fn, fp, tn, tp = get_confution_matrix(E1, Ea, Ea_hat, omega)
    acc = (tp + tn) / omega * 1.0
    f01 = get_f_beta(E1, Ea, Ea_hat, 0.1)
    f033 = get_f_beta(E1, Ea, Ea_hat, 1.0 / 3)
    f3 = get_f_beta(E1, Ea, Ea_hat, 3)
    f10 = get_f_beta(E1, Ea, Ea_hat, 10)
    measures = [acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp]
    measure_names = ['acc', 'bm', 'ec', 'f01', 'f033', 'f1', 'f10', 'f3', 'fn', 'fp', 'ics', 'mcc', 'mk', 's3',
                     'tn', 'tp']
    return pd.Series(index=measure_names, data=measures)


def get_f_beta(E1, Ea, Ea_hat, beta):
    f_beta = ((1 + beta ** 2) * Ea) / float(E1 + beta ** 2 * Ea_hat)
    return f_beta


def get_confution_matrix(E1, Ea, Ea_hat, omega):
    tp = Ea
    fp = E1 - Ea
    fn = Ea_hat - Ea
    tn = omega - (tp + fp + fn)
    return fn, fp, tn, tp


def compute_measures(E1, Ea, Ea_hat, omega):
    ec = Ea / E1
    ics = Ea / Ea_hat
    s3 = Ea / (E1 + Ea_hat - Ea)
    f1 = (2 * Ea) / (E1 + Ea_hat)
    bm = (omega * Ea - E1 * Ea_hat + 1.0) / ((omega - Ea_hat) * Ea_hat + 1.0)
    mk = (omega * Ea - E1 * Ea_hat + 1.0) / ((omega - E1) * E1 + 1.0)
    mcc = (omega * Ea - E1 * Ea_hat + 1.0) / ((E1 * Ea_hat * (omega - E1) * (omega - Ea_hat)) ** 0.5 + 1.0)
    # mcc_geometric = (bm + mk) / 2.0
    return ec, f1, ics, s3, bm, mk, mcc


def generate_G1(Y0, nodes_to_delete, edges_to_modify):
    G1 = nx.Graph(Y0)
    v1 = list(G1.nodes)
    deg_view = G1.degree()
    deg_median = np.median([d for n, d in deg_view])
    node_with_low_deg = [n for n, d in deg_view if d < deg_median]
    G1.remove_nodes_from(random.sample(node_with_low_deg, min(nodes_to_delete, len(node_with_low_deg))))

    if edges_to_modify > 0:
        total_edges_to_change = int(len(v1) * edges_to_modify)
        edges_to_add = random.randint(0, total_edges_to_change)
        edges_to_remove = total_edges_to_change - edges_to_add
        add_random_edges(G1, edges_to_add)
        remove_random_edges(G1, edges_to_remove)
    return G1


def get_align_edges(G1, G2, A):
    Ea_hat = len(G2.subgraph(A.values).edges) * 1.0
    Ea = sum([1 if G2.has_edge(A[v1], A[v2]) else 0 for v1, v2 in G1.edges]) * 1.0
    E1 = len(G1.edges) * 1.0
    return E1, Ea, Ea_hat


def update_Ea_swap(v1, v2, u1, u2, G1, G2, A):
    res = 0
    for n in G1.neighbors(v1):
        res -= 1 if G2.has_edge(u1, A[n]) else 0
        res += 1 if G2.has_edge(u2, A[n]) else 0

    for n in G1.neighbors(v2):
        res -= 1 if G2.has_edge(u2, A[n]) else 0
        res += 1 if G2.has_edge(u1, A[n]) else 0

    res += 2 if G2.has_edge(u1, u2) and G1.has_edge(v1, v2) else 0
    return res


def generate_alignments_for_objective(G2s, edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, objective,
                                      output_file_name='alignmnets_y0_com_bal'):
    iter_counter = 1
    for i in range(full_iterations):
        for j, G2 in enumerate(G2s):
            for nodes_to_delete in nodes_to_delete_list:
                for edge_frac_To_modify in edge_frac_To_modify_list:

                    G2_name = G2.graph['name']
                    print('full iteration {}/{}, using G2 {}, nodes_to_delete {}, edge_to_delete {} {}/{}'.format(
                        str(i + 1), full_iterations, G2_name, nodes_to_delete, edge_frac_To_modify, iter_counter,
                        str(len(G2s) * len(nodes_to_delete_list) * len(edge_frac_To_modify_list) * full_iterations)))
                    iter_counter += 1
                    G1 = generate_G1(G2, nodes_to_delete, edge_frac_To_modify)
                    G1_name = '{}_del_nodes_{}_modify_{}_iter_{}'.format(G2.graph['name'], nodes_to_delete,
                                                                         edge_frac_To_modify, i)
                    V1 = list(G1.nodes)
                    V2 = list(G2.nodes)
                    E2 = len(list(G2.edges))

                    write_leda_to_file(G1, G1_name, 'output_alignment/graphs')

                    true_A = pd.Series(data=V1, index=V1)

                    A = generate_random_alignment(V1, V2)

                    E1, Ea, Ea_hat = get_align_edges(G1, G2, A)
                    Ea_hat_star = Ea_hat
                    Ea_star = Ea
                    omega = (len(V1) * (len(V1) - 1)) / 2.0

                    nc = get_nc(A, true_A)
                    measures = get_measures(E1, Ea, Ea_hat, omega)
                    acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures
                    row = (
                        G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, '', 0, '', 0,
                        '', 0, nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk, mcc, f01, f033, f3, f10)
                    rows = [row]

                    original_Ea = Ea
                    original_Ea_hat = Ea_hat
                    objective_function = measures[objective]
                    set_free_nodes = set(G2.nodes) - set(G1.nodes)
                    correct_nodes = nc * len(V1)
                    tries = 10000000
                    start = timeit.default_timer()
                    for x in range(tries):
                        if random.random() < 0:
                            v1 = random.choice(V1)
                            u1 = A[v1]
                            u2 = random.choice(list(set_free_nodes))

                            ea_add = update_Ea_change(v1, u1, u2, G1, G2, A)
                            ea_hat_add = update_Ea_hat(u1, u2, G2, set_free_nodes)

                            Ea += ea_add
                            Ea_hat += ea_hat_add

                            measures = get_measures(E1, Ea, Ea_hat, omega)
                            acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures

                            if measures[objective] > objective_function:
                                change(A, v1, u2, set_free_nodes)
                                objective_function = measures[objective]

                                correct_nodes += update_nc(v1, u1, u2)
                            else:
                                Ea -= ea_add
                                Ea_hat -= ea_hat_add
                        else:
                            v1, v2 = random.sample(V1, 2)
                            u1, u2 = A[v1], A[v2]

                            ea_add = update_Ea_swap(v1, v2, u1, u2, G1, G2, A)

                            Ea += ea_add

                            measures = get_measures(E1, Ea, Ea_hat, omega)
                            acc, bm, ec, f01, f033, f1, f10, f3, fn, fp, ics, mcc, mk, s3, tn, tp = measures

                            if measures[objective] > objective_function:
                                # swap(A, u1, u2, v1, v2)
                                A[v1], A[v2] = u2, u1
                                objective_function = measures[objective]

                                correct_nodes += update_nc(v1, u1, u2)
                                correct_nodes += update_nc(v2, u2, u1)
                            else:
                                Ea -= ea_add

                        nc = correct_nodes / len(V1)

                        # assert nc == get_nc(A, true_A)

                        row = (
                            G1_name, len(V1), E1, G2_name, len(V2), E2, omega, Ea_hat_star, Ea, Ea_hat, v1,
                            G1.degree[v1], u1,
                            G2.degree[u1], u2, G2.degree[u2], nc, tp, fp, tn, fn, ec, ics, s3, f1, acc, bm, mk,
                            mcc,
                            f01,
                            f033, f3, f10)
                        rows.append(row)

                        stop = timeit.default_timer()
                        print('\r sub iteration {}/{}, obj {} time: {}sec'.format(x, tries, objective_function,
                                                                                  stop - start),
                              end='')

                    df = pd.DataFrame(rows,
                                      columns=['G1', 'V1', 'E1', 'G2', 'V2', 'E2', 'omega', 'Ea_hat*', 'Ea',
                                               'Ea_hat', 'G1 changed node', 'Degree of changed node',
                                               'G2 old node', 'Degree of old node', 'G2 new node',
                                               'Degree of new node', 'NC', 'TP', 'FP', 'TN', 'FN', 'EC', 'ICS',
                                               'S3', 'F1', 'ACC', 'BM', 'MK', 'MCC', 'F0.1', 'F0.33', 'F3',
                                               'F10'])

                    output_path = 'output_alignment/%s.csv' % output_file_name
                    if os.path.isfile(output_path):
                        df.to_csv(output_path, mode='a', header=None)
                    else:
                        df.to_csv(output_path)
                    print()


def update_nc(v, u_old, u_new):
    res = 0
    if v == u_old:
        res -= 1
    if v == u_new:
        res += 1
    return res


def swap(A, u1, u2, v1, v2):
    A[v1], A[v2] = u2, u1


def generate_random_alignment(V1, V2):
    return pd.Series(data=random.sample(V2, len(V1)), index=V1)


def star_barabasi_albert_graph(n, m, seed=None):
    if m < 1 or m >= n:
        raise nx.NetworkXError("BA network must have m >= 1"
                               " and m < n, m = %d, n = %d" % (m, n))
    if seed is not None:
        random.seed(seed)

    # Add m initial nodes (m0 in barabasi-speak)
    G = nx.star_graph(m)
    # Target nodes for new edges
    targets = list(range(m))
    # List of existing nodes, with nodes repeated once for each adjacent edge
    repeated_nodes = []
    # Start adding the other n-m nodes. The first node is m.
    source = m
    while source < n:
        # Add edges to m nodes from the source.
        G.add_edges_from(zip([source] * m, targets))
        # Add one node to the list for each new edge just created.
        repeated_nodes.extend(targets)
        # And the new node "source" has m edges to add to the list.
        repeated_nodes.extend([source] * m)
        # Now choose m unique nodes from the existing nodes
        # Pick uniformly from repeated_nodes (preferential attachement)
        targets = _random_subset(repeated_nodes, m)
        source += 1
    return G


def _random_subset(seq, m):
    """ Return m unique elements from seq.

    This differs from random.sample which can return repeated
    elements if seq holds repeated elements.
    """
    targets = set()
    while len(targets) < m:
        x = random.choice(seq)
        targets.add(x)
    return targets


def main():
    '''
    graphs to generate: Barabasi, Wattz, Erdos size 1000 nodes
    barabasi_albert: n: 1000, m: [9, 500] and ba(1000, 9) complement
    erdos: n: 1000, m: [8500, 250000, 491000]
    wattz: n: 1000, k: [18, 500, 982], p: 0.5

    delete [10, 50, 100] nodes
    full_iterations: 2
    shuffles: 2
    edge_frac_To_modify_list = [0.01, 0.05, 0.1]
    nodes_to_delete_list = [10, 50, 100]
    :return:
    '''

    # nx.gnm_random_graph(n, m) erdos
    # nx.watts_strogatz_graph()

    full_iterations = 5
    shuffles = 5
    # edge_frac_To_modify_list = [0.01, 0.05, 0.1]
    edge_frac_To_modify_list = [0.0]
    # nodes_to_delete_list = [10, 50, 100]
    nodes_to_delete_list = [100]

    G2s = []
    # G2s += get_barabasi_graphs()
    # G2s += get_erdos_graphs()
    # G2s += get_wattz_graphs()

    G2s += syeast0_graphs()

    # jobs = []
    # for G2 in G2s:
    #     experiment_name = '{}_iter_{}_shuffle_{}'.format(G2.graph['name'], full_iterations, shuffles)
    #     # args = ([G2], edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, shuffles, experiment_name)
    #     args = ([G2], edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, 'mcc', experiment_name)
    #     # p = Process(target=generate_alignments, args=args)
    #     p = Process(target=generate_alignments_for_objective, args=args)
    #     jobs.append(p)
    #     p.start()
    #
    # for p in jobs:
    #     p.join()
    generate_alignments(G2s, edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, shuffles, 'test')
    # generate_alignments_for_objective(G2s[:1], edge_frac_To_modify_list, full_iterations, nodes_to_delete_list, 'mcc',
    #                                   'test_hill')

    pass


def generate_true_alignment(G, g_name, output_path):
    graph_true_alignment_path = os.path.join(output_path, '%s.align' % g_name)
    # if not os.path.exists(graph_true_alignment_path):
    #     os.makedirs(graph_true_alignment_path)
    pd.DataFrame(zip(list(G.nodes), list(G.nodes))).to_csv(graph_true_alignment_path, header=None, sep='\t',
                                                           index=False)
    # with open(graph_true_alignment_path, 'wb') as ouput_file:
    #     for node in list(G.nodes):
    #         ouput_file.write(node + '\t' + node + '\n')


def write_grapha_and_true_alignment(G1_add, G1_add_name, output_path):
    generate_true_alignment(G1_add, G1_add_name, output_path)
    write_leda_to_file(G1_add, G1_add_name, os.path.join(output_path, '{0}/'.format(G1_add_name)))


def generate_graphs_for_beta_experiments():
    output_path = 'roni_test_graphs_and_data/'
    nodes_to_delete = 100
    # np.arange(0.140, 0.162, 0.002)
    g2_names = set()
    for i in range(10):
        for target_dense in [0.017]:
            G2_name = 'syeast0'
            input_g2_file = 'networks/{0}/{0}.gw'.format(G2_name)
            G2 = nx.read_leda(input_g2_file)
            G2 = nx.complement(G2)
            G2_name = 'syeast0_complement'
            # show_graph_degree_histogram(G2)
            G2.graph['name'] = G2_name
            V2 = len(G2.nodes)
            g2_dens = nx.density(nx.read_leda(input_g2_file))
            print(g2_dens)
            # G2_base_name = G2_name.split('_')[0]
            k = int(V2 * (V2 - 1) / 2 * target_dense) - len(G2.edges)
            # # G2 = add_random_edges(G2, k)
            # G2 = remove_random_edges(G2, -k)
            # G2_base = '{0}_dens_{1:.4f}'.format(G2_name, nx.density(G2))
            G2_base = G2_name
            g2_names.add(G2_base)
            write_leda_to_file(G2, G2_base, os.path.join(output_path, '{0}/'.format(G2_base)))
            G_name_template = '{}_node_del_{}_{}_edges_{}_iter_{}'
            print('G2 density {}'.format(nx.density(G2)))
            for modify in [0.30]:
                G1_name = G_name_template.format(G2_base, nodes_to_delete, 'add', modify, i)
                G1 = generate_modified_graph(G2, add_random_edges, G1_name, modify, nodes_to_delete, output_path)

                # modify = 0.5
                # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'remove', modify, i)
                # generate_modified_graph(G2, remove_random_edges, G1_name, modify, nodes_to_delete, output_path)
                #
                # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'add', modify, i)
                # generate_modified_graph(G2, add_random_edges, G1_name, modify, nodes_to_delete, output_path)
                #
                # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'rewire', modify, i)
                # generate_modified_graph(G2, rewire_edges, G1_name, modify, nodes_to_delete, output_path)
    print(' '.join(g2_names))


def get_erdos(n, density):
    assert 0 <= density <= 1
    k = int(n * (n - 1) / 2 * density)
    ER = nx.gnm_random_graph(n, k)
    ER.graph['name'] = 'er_{}_{}'.format(n, k)
    return ER


def get_barabasi_by_dense(n, density):
    assert 0 <= density <= 0.5
    k = int(n * (n - 1) / 2 * density)
    m = int(round((n - (n ** 2 - 4 * k) ** 0.5) / 2))
    BA = create_barabasi(n, m)
    print(nx.density(BA))
    return BA


def create_barabasi(n, m):
    BA = nx.barabasi_albert_graph(n, m)
    BA.graph['name'] = 'ba_{}_{}'.format(n, m)
    print(nx.density(BA))
    return BA


def get_star_barabasi_by_dense(n, density):
    assert 0 <= density <= 0.5
    k = int(n * (n - 1) / 2 * density)
    m = int(round((n - (n ** 2 - 4 * k) ** 0.5) / 2))
    BA = create_star_barabasi(n, m)
    print(nx.density(BA))
    return BA


def create_star_barabasi(n, m):
    BA = nx.barabasi_albert_graph(n, m)
    BA.graph['name'] = 'ba_star_{}_{}'.format(n, m)
    print(nx.density(BA))
    return BA


def generate_graphs_for_beta_experiments_ba_er():
    output_path = 'roni_test_graphs_and_data/'
    n = 1000

    g2_names = list()
    for j in range(1):
        G2s = []
        for target_dense in [0.02]:
            # G2 = get_star_barabasi_by_dense(n, target_dense)
            m = 10
            p = 0.2
            G2 = nx.powerlaw_cluster_graph(n, m, p)
            show_graph_degree_histogram(G2)
            G2.graph['name'] = 'powerlaw_cluster_n{}_m{}_p{}'.format(n, m, p)
            write_grapha_and_true_alignment(G2, G2.graph['name'], output_path)
            G2s.append(G2)

            # BA_complement = get_barabasi_by_dense_complement(n, target_dense)
            G2_complement = nx.complement(G2)
            G2_complement.graph['name'] = G2.graph['name'] + '_comp'

            G2s.append(G2_complement)
            for G2 in G2s:
                G2_name = G2.graph['name'] + '_ver_{}'.format(j)
                g2_dens = nx.density(G2)
                print('G2 {}, dense {} '.format(G2_name, g2_dens))
                g2_names.append(G2_name)
                write_leda_to_file(G2, G2_name, os.path.join(output_path, '{0}/'.format(G2_name)))
                print('G2 density {}'.format(nx.density(G2)))
                G_name_template = '{}_node_del_{}_{}_{}'
                for nodes_to_delete in [100]:
                    for modify in [0.0]:
                        # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'add', modify)
                        # G1 = generate_modified_graph(G2, add_random_edges, G1_name, modify, nodes_to_delete,
                        #                              output_path)
                        #
                        # G1_name = G_name_template.format(G2_name, nodes_to_delete, 'remove', modify)
                        # G1 = generate_modified_graph(G2, remove_random_edges, G1_name, modify, nodes_to_delete,
                        #                              output_path)

                        G1_name = G_name_template.format(G2_name, nodes_to_delete, 'rewire', modify)
                        G1 = generate_modified_graph(G2, maor_rewire_edges, G1_name, modify, nodes_to_delete,
                                                     output_path)

    print(' '.join(g2_names))


def get_barabasi_by_dense_complement(n, target_dense):
    BA = get_barabasi_by_dense(n, target_dense)
    BA_complement = nx.complement(BA)
    BA_complement.graph['name'] = BA.graph['name'] + '_comp'
    return BA_complement


def show_graph_degree_histogram(G):
    degree_sequence = np.array([d for n, d in G.degree()])
    # bins = np.logspace(degree_sequence[0], degree_sequence[-1], 10, base=2)
    plt.hist(degree_sequence)  # arguments are passed to np.histogram
    plt.title("Histogram log 2")
    plt.show()


def rewire_edges(G1, k):
    return nx.double_edge_swap(nx.Graph(G1), k, max_tries=k * 100)


def maor_rewire_edges(G1, k):
    G1 = remove_random_edges(G1, k)
    return add_random_edges(G1, k)


def generate_modified_graph(G2, edge_modify_fn, G1_name, modify, nodes_to_delete, output_path):
    G1 = generate_G1(G2, nodes_to_delete, 0)
    V1 = len(G1.nodes)
    E1 = len(G1.edges)
    edges_to_modify = int(modify * E1)
    G1_remove = edge_modify_fn(nx.Graph(G1), edges_to_modify)
    density = nx.density(G1)
    print('G1 n: {}, m: {}, density {}, {}'.format(V1, len(G1_remove.edges), density, G1_name))
    write_grapha_and_true_alignment(G1_remove, G1_name, output_path)
    return G1


def convert_contak_grapgs_to_leda():
    output_path = 'roni_test_graphs_and_data/'
    graph_name = 'bn-mouse_retina_1'
    tsv_file = 'rat_network/%s.edges' % graph_name
    df = pd.read_csv(tsv_file, sep=' ', header=None)
    edge_list = []
    for index, row in df.iterrows():
        v1, v2 = row[0], row[1]
        edge_list.append((v1, v2))
        pass
    G2 = nx.from_edgelist(edge_list)
    G2.graph['name'] = graph_name

    print(list(nx.selfloop_edges(G2)))
    G2.remove_edges_from(nx.selfloop_edges(G2))
    print('density', nx.density(G2))
    print('edges', len(G2.edges()))
    print('nodes', len(G2.nodes()))
    G2s = []
    G2s.append(G2)
    # exit()
    # G2_complement = nx.complement(G2)
    # G2_complement.graph['name'] = G2.graph['name'] + '_comp'
    # G2s.append(G2_complement)

    for G2 in G2s:
        G2_name = G2.graph['name']
        write_leda_to_file(G2, G2_name, os.path.join(output_path, '{0}/'.format(G2_name)))
        write_grapha_and_true_alignment(G2, G2_name, output_path)

        G2_complement = nx.complement(G2)
        G2_complement.graph['name'] = G2.graph['name'] + '_comp'
        write_leda_to_file(G2_complement, G2_complement.graph['name'],
                           os.path.join(output_path, '{0}/'.format(G2_name)))
        write_grapha_and_true_alignment(G2_complement, G2_complement.graph['name'], output_path)

        V2 = len(G2.nodes)
        G_name_template = '{}_node_del_{}_{}_{}'
        for nodes_to_delete in [int(V2 * 0.1)]:
            for modify in [0.1]:
                G1_name = G_name_template.format(G2_name, nodes_to_delete, 'rewire', modify)
                G1 = generate_modified_graph(G2, maor_rewire_edges, G1_name, modify, nodes_to_delete,
                                             output_path)

                G1_comp_name = G_name_template.format(G2_name + '_comp', nodes_to_delete, 'rewire', modify)
                write_grapha_and_true_alignment(nx.complement(G1), G1_comp_name, output_path)


def convert_snaps_grapgs_to_leda():
    output_path = 'roni_test_graphs_and_data/'
    input_path = 'networks/'
    G1 = nx.read_leda(os.path.join(input_path, '{0}/{0}.gw'.format('syeast0')))
    G2 = nx.read_leda(os.path.join(input_path, '{0}/{0}.gw'.format('yeast')))
    assert isinstance(G2, nx.Graph)
    node_core_dict = nx.core_number(G2)
    g2_core_df = pd.DataFrame(node_core_dict.items(), columns=['node', 'core'])
    g2_core_count = g2_core_df.groupby('core')['node'].count().reset_index()
    # print(g2_core_count)

    g1_nodes_set = set(G1.nodes)
    node_core_dict = {n: c for n, c in node_core_dict.iteritems() if n in g1_nodes_set}
    g1_subset_of_g2_core_df = pd.DataFrame(node_core_dict.items(), columns=['node', 'core'])
    g1_subset_g2_core_count = g1_subset_of_g2_core_df.groupby('core')['node'].count().reset_index()
    # print(g1_subset_g2_core_count)

    assert isinstance(g2_core_count, pd.DataFrame)
    res = g2_core_count.merge(g1_subset_g2_core_count, on='core')
    print(res)

    # G2_sub_graph = G2.subgraph(G1.nodes)
    # nx.write_gexf(G2_sub_graph, os.path.join(output_path, "yeast_sub_graph_of_syeast0.gexf"))

    # tsv_file = 'comm-f2f-Resistance/email-Eu-core.txt'
    # df = pd.read_csv(tsv_file, sep=' ', header=None)
    # edge_list = list(zip(df[0], df[1]))
    # G2 = nx.from_edgelist(edge_list, nx.DiGraph)
    # G2 = G2.to_undirected()
    # G2.remove_edges_from(G2.selfloop_edges())
    # G2_name = 'email_eu_core'
    # write_leda_to_file(G2, G2_name, os.path.join(output_path, '{0}/'.format(G2_name)))
    # write_grapha_and_true_alignment(G2, G2_name, output_path)

    # degree_sequence = sorted([d for n, d in G2.degree()], reverse=True)  # degree sequence
    # G_name_template = '{}_node_del_{}_{}_{}'
    # for nodes_to_delete in [100]:
    #     for modify in [0.1]:
    #         G1_name = G_name_template.format(G2_name, nodes_to_delete, 'rewire', modify)
    #         G1 = generate_modified_graph(G2, maor_rewire_edges, G1_name, modify, nodes_to_delete,
    #                                      output_path)

    # G1_name = 'email_eu_core_node_del_100_rewire_0.1'
    # G1 = nx.read_leda(os.path.join(output_path, '{0}/{0}.gw'.format(G1_name)))
    # G2_name = 'email_eu_core'
    # G2 = nx.read_leda(os.path.join(output_path, '{0}/{0}.gw'.format(G2_name)))
    #
    # G1 = nx.complement(G1)
    # G1_name += '_complement'
    # G2 = nx.complement(G2)
    # G2_name += '_complement'
    #
    # write_leda_to_file(G1, G1_name, os.path.join(output_path, '{0}/'.format(G1_name)))
    # write_grapha_and_true_alignment(G1, G1_name, output_path)
    # print('g1 comp density: {}'.format(nx.density(G1)))
    #
    # write_leda_to_file(G2, G2_name, os.path.join(output_path, '{0}/'.format(G2_name)))


def create_bfs_graphs():
    output_path = 'roni_test_graphs_and_data/'
    input_path = 'networks/'
    G1 = nx.read_leda(os.path.join(input_path, '{0}/{0}.gw'.format('syeast0')))
    G2 = nx.read_leda(os.path.join(input_path, '{0}/{0}.gw'.format('yeast')))
    G2.graph['name'] = 'yeast'

    G2s = [G2]

    for G2 in G2s:
        G2_name = G2.graph['name']
        write_leda_to_file(G2, G2_name, os.path.join(output_path, '{0}/'.format(G2_name)))
        write_grapha_and_true_alignment(G2, G2_name, output_path)

        G2_complement = nx.complement(G2)
        G2_complement.graph['name'] = G2.graph['name'] + '_comp'
        write_leda_to_file(G2_complement, G2_complement.graph['name'],
                           os.path.join(output_path, '{0}/'.format(G2_name)))
        write_grapha_and_true_alignment(G2_complement, G2_complement.graph['name'], output_path)

        assert isinstance(G2, nx.Graph)
        while True:
            random_node = random.choice(list(G2.nodes))
            T = nx.bfs_tree(G2, random_node)
            bfs_nodes = list(T.nodes)
            if len(bfs_nodes) > len(G2.nodes) * 0.75:
                break

        for node_portion in [0.25, 0.5, 0.75]:
            g1_nodes = bfs_nodes[:int(len(G2.nodes) * node_portion)]
            G_name_template = '{}_from_{}_nodes_{}_{}_{}'
            nodes_count = len(g1_nodes)
            for modify in [0.01, 0.05, 0.1]:
                G1_name = G_name_template.format(G2_name, random_node, nodes_count, 'rewire', modify)
                E1 = len(G1.edges)
                edges_to_modify = int(modify * E1)

                G1 = nx.Graph(G2.subgraph(g1_nodes))
                G1 = maor_rewire_edges(G1, edges_to_modify)
                write_grapha_and_true_alignment(G1, G1_name, output_path)

                G1_comp_name = G_name_template.format(G2_name + '_comp', random_node, nodes_count, 'rewire', modify)
                write_grapha_and_true_alignment(nx.complement(G1), G1_comp_name, output_path)


def print_densities():
    input_path = 'roni_test_graphs_and_data/'
    for graph_name in os.listdir(input_path):
        if os.path.isdir(os.path.join(input_path, graph_name)):
            G2 = nx.read_leda(os.path.join(input_path, '{0}/{0}.gw'.format(graph_name)))
            density = nx.density(G2)
            if density < (3 - 2 ** 0.5) / 4:
                b = 1 / (4 * (1 - density))
            else:
                b = (1 + (-7 + 24 * density - 16 * density ** 2) ** 0.5) / (4 * (1 - density))
            print('{} density: {}, Fb# = F {}'.format(graph_name, density, b))


def print_F_beta_star():
    graph_names = [
        ('econ-psmigr1_node_del_314_rewire_0.1', 'econ-psmigr1'),
        ('econ-psmigr1_comp_node_del_314_rewire_0.1', 'econ-psmigr1_comp'),
        ('email_eu_core_node_del_100_rewire_0.1', 'email_eu_core'),
        ('email_eu_core_node_del_100_rewire_0.1_complement', 'email_eu_core_complement'),
        ('facebook_combined_comp_node_del_403_rewire_0.1', 'facebook_complement'),
        ('facebook_combined_node_del_403_rewire_0.1', 'facebook'),
        ('syeast0', 'yeast'),
        ('syeast0_complement', 'yeast_comp'),
    ]

    input_path = 'roni_test_graphs_and_data/'
    for g1_name, g2_name in graph_names:
        G1 = nx.read_leda(os.path.join(input_path, '{0}/{0}.gw'.format(g1_name)))
        G2 = nx.read_leda(os.path.join(input_path, '{0}/{0}.gw'.format(g2_name)))

        V1 = list(G1.nodes())
        true_A = pd.Series(data=V1, index=V1)
        E1, Ea, Ea_hat = get_align_edges(G1, G2, true_A)
        Ea_hat_star = Ea_hat
        omega = (len(V1) * (len(V1) - 1)) / 2.0

        TP = Ea
        FP = E1 - Ea
        FN = Ea_hat - Ea
        TN = omega - (TP + FP + FN)
        ilia_beta = (E1 * TN) / ((omega - E1) * TP)
        D1 = nx.density(G1)
        D2 = nx.density(G2)
        print('################')
        print('G1: {} \t--\t G2: {}'.format(g1_name, g2_name))
        print('D1: {} \t--\t D2: {}'.format(D1, D2))
        print('Ea_hat_star: {}'.format(Ea_hat_star))
        print('F_beta_star: {}'.format(E1 / float(Ea_hat_star)))
        print('F#: {}'.format(
            (1 / (4 * (1 - D1))) if D1 < 0.42 else (1 + (-7 + 24 * D1 - 16 * D1 ** 2) ** 0.5) / (4 * (1 - D1))))
        print('F_beta Ilia measureT: {}'.format(ilia_beta))


# def create_rat_network():
#     G1_name = G_name_template.format(G2_name, nodes_to_delete, 'rewire', modify)
#     G1 = generate_modified_graph(G2, maor_rewire_edges, G1_name, modify, nodes_to_delete,
#                                  output_path)

import matplotlib.pyplot as plt
from tqdm import tqdm


def compute_shell():
    graph_name = 'bn-mouse_retina_1'
    tsv_file = 'rat_network/%s.edges' % graph_name
    # df = pd.read_csv(tsv_file, sep=' ', header=None)
    # edge_list = []
    # for index, row in tqdm(df.iterrows(), total=len(df)):
    #     v1, v2 = row[0], row[1]
    #     edge_list.append((v1, v2))
    #     pass
    # G2 = nx.from_edgelist(edge_list)
    # G2 = nx.read_leda('networks/syeast0/syeast0.gw')
    G2 = nx.from_edgelist(tqdm(read_edgelist('rat_network/soc-LiveJournal1.txt')))
    print('max deg: {}'.format(max(G2.degree.values())))
    # print('compute k shell')
    # v_cores = pd.Series(nx.core_number(G2))
    # print('Max core: {}'.format(v_cores.max()))
    # v_cores.plot.hist()
    # plt.show()


import csv


def read_edgelist(path):
    with open(path) as f:
        for edge in csv.reader(f, delimiter='\t'):
            yield (edge[0], edge[1])


if __name__ == "__main__":
    # main()
    # generate_graphs_for_beta_experiments()
    # generate_graphs_for_beta_experiments_ba_er()
    # convert_contak_grapgs_to_leda()
    # convert_snaps_grapgs_to_leda()
    # create_bfs_graphs()
    # print_densities()
    print_F_beta_star()
    # compute_shell()
