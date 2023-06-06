import networkx as nx
import itertools
import math
import pandas as pd


def eval_scores_for_graphs(graph1, graph2):
    assert isinstance(graph1, nx.Graph)
    assert isinstance(graph2, nx.Graph)
    TP = 0.0
    FP = 0.0
    TN = 0.0
    FN = 0.0
    V1 = list(graph1.nodes())
    V2 = list(graph2.nodes())
    assert set(V1).issubset(set(V2))
    for i, v in enumerate(V1):
        for u in V1[i + 1:]:
            if graph1.has_edge(v, u):
                if graph2.has_edge(v, u):
                    TP += 1
                else:
                    FP += 1
            else:
                if graph2.has_edge(v, u):
                    FN += 1
                else:
                    TN += 1
    return TP, FP, TN, FN


def get_objective_functions(TP, FP, TN, FN):
    omega = TP + FP + TN + FN
    S3 = TP / float(TP + FP + FN)
    EC = TP / float(TP + FP)
    ICS = TP / float(TP + FN)
    ACC = (TP + TN) / float(omega)
    F1 = (2 * TP) / float(2 * TP + FP + FN)
    GED = FP + FN
    MCC = (TN * TP - FP * FN) / float(math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
    return S3, EC, ICS, ACC, F1, GED, MCC


def get_acc_weighted_by_betas(TP, FP, TN, FN, betas):
    results = []
    for beta in betas:
        results.append((beta * TP + TN) / float(beta * (TP + FP) + (TN + FN)))
    return results


g1_name = 'syeast0_delete_0.15nodes'
g2_name = 'syeast0_rewire_37edges'
G1 = nx.read_leda('networks/{0}/{0}.gw'.format(g1_name))
G2 = nx.read_leda('output_graphs/{0}/{0}.gw'.format(g2_name))

TP, FP, TN, FN = eval_scores_for_graphs(G1, G2)
output_df = pd.DataFrame()
output_df['G1'] = [g1_name]
output_df['G2'] = [g2_name]
output_df['TP'] = [TP]
output_df['FP'] = [FP]
output_df['TN'] = [TN]
output_df['FN'] = [FN]
output_df.to_csv('eval_equality.csv')
print('TP: {}, FP: {}, TN: {}, FN: {}'.format(TP, FP, TN, FN))

# G1_name = 'syeast0'
# G1 = nx.read_leda('networks/syeast0/%s.gw' % G1_name)
# target_graphs_names = ['syeast05', 'syeast10', 'syeast15', 'syeast20', 'syeast25', 'syeast-05', 'syeast-10',
#                        'syeast-15', 'syeast-20', 'syeast-25']
# rows = []
# betas = [0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 100]
# for G2_name in target_graphs_names:
#     G2 = nx.read_leda('networks/{0}/{0}.gw'.format(G2_name))
#     assert isinstance(G1, nx.Graph)
#     assert isinstance(G2, nx.Graph)
#     TP, FP, TN, FN = eval_scores_for_graphs(G1, G2)
#     S3, EC, ICS, ACC, F1, GED, MCC = get_objective_functions(TP, FP, TN, FN)
#     weighted_accs = get_acc_weighted_by_betas(TP, FP, TN, FN, betas)
#     omega = TP + FP + TN + FN
#     min_acc = 1 - float(G1.number_of_edges() + G2.number_of_edges()) / float(omega)
#     nor_acc = (ACC - min_acc) / float(1 - min_acc)
#     rows.append([G1_name, G2_name, TP, FP, TN, FN, S3, EC, ICS, ACC, nor_acc, F1, GED, MCC] + weighted_accs)
#
# weighted_acc_columns = ['weighted_acc_' + str(b) for b in betas]
# columns = ['G1_name', 'G2_name', 'TP', 'FP', 'TN', 'FN', 'S3', 'EC', 'ICS', 'ACC', 'nor_acc', 'F1', 'GED',
#            'MCC'] + weighted_acc_columns
# eval_true_alignment_df = pd.DataFrame(rows, columns=columns)
# eval_true_alignment_df.to_csv('eval_true_alignment_with_weighted.csv')
