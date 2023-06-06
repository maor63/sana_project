import random
import numpy as np
from collections import Counter
from deap import algorithms
from deap import base
from deap import creator
import networkx as nx
from deap import tools
import pandas as pd
from functools import partial

G1 = nx.read_leda('networks/{0}/{0}.gw'.format('syeast0'))
# G2 = nx.read_leda('networks/{0}/{0}.gw'.format('yeast'))
G2 = nx.read_leda('networks/{0}/{0}.gw'.format('syeast0'))
assert isinstance(G1, nx.Graph)
assert isinstance(G2, nx.Graph)
nodes = list(G1.nodes)[:500]
G1 = G1.subgraph(nodes)
G2 = G2.subgraph(nodes)

G1_NODES_MAPPING = pd.Series(index=list(G1.nodes), data=range(len(G1.nodes)))

G2_NODES_MAPPING = pd.Series(index=list(G2.nodes), data=range(len(G2.nodes)))
G1 = nx.relabel_nodes(G1, G2_NODES_MAPPING)
G2 = nx.relabel_nodes(G2, G2_NODES_MAPPING)

NB_G1_NODES = len(G1.nodes)
NB_G2_NODES = len(G2.nodes)
G2_NODES_SET = set(range(len(G2.nodes)))

E1 = len(G1.edges) * 1.0

true_A = pd.Series(data=G1_NODES_MAPPING[list(G1.nodes)].values, index=G1_NODES_MAPPING[list(G1.nodes)].values)


def get_align_edges(G1, G2, individual):
    A = individual
    Ea_hat = len(G2.subgraph(A).edges) * 1.0
    Ea = sum([1 if G2.has_edge(A[v1], A[v2]) else 0 for v1, v2 in G1.edges]) * 1.0
    return E1, Ea, Ea_hat


def ind_to_alignment(G1, individual):
    A = G2_NODES_MAPPING[individual]
    A.index = list(G1.nodes)
    return A


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


creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
# Helper function for defining NetworkAlignment permutation
toolbox.register("permutation", random.sample, range(NB_G2_NODES), k=NB_G1_NODES)
# Definition for how to initialize an individual (Alignment)
toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.permutation)
# Definition for the container of all individuals (Alignment)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


def eval_alignment(individual):
    size = len(individual)
    if hasattr(individual, 'Ea'):
        E1, Ea, Ea_hat = individual.E1, individual.Ea, individual.Ea_hat
    else:
        E1, Ea, Ea_hat = update_ind_parameters(individual)
    s3 = calculate_s3(E1, Ea, Ea_hat)
    individual.S3 = s3
    return s3,


def update_ind_parameters(individual):
    E1, Ea, Ea_hat = get_align_edges(G1, G2, individual)
    individual.E1, individual.Ea, individual.Ea_hat = E1, Ea, Ea_hat
    return E1, Ea, Ea_hat


def calculate_s3(E1, Ea, Ea_hat):
    return Ea / (E1 + Ea_hat - Ea)


def get_nc(alignment, true_alignment):
    correct = sum(1 for v in alignment if alignment[v] == true_alignment[v])
    return correct * 1.0 / len(alignment)


def mut_change_swap(individual, indpb, change_prob):
    for j in range(random.randint(1, 3)):
        if random.random() > change_prob or NB_G2_NODES == NB_G1_NODES:
            v1, v2 = random.sample(range(len(individual)), 2)
            u1, u2 = individual[v1], individual[v2]

            ea_add = update_Ea_swap(v1, v2, u1, u2, G1, G2, individual)
            individual.Ea += ea_add

            individual[v1], individual[v2] = individual[v2], individual[v1]
            # e1, ea, ea_hat = get_align_edges(G1, G2, individual)
            # assert ea == individual.Ea
            # if not(calculate_s3(individual.E1, individual.Ea, individual.Ea_hat) > individual.S3
            #        or random.random() < indpb):
            #     individual[v1], individual[v2] = individual[v2], individual[v1]
            #     individual.Ea -= ea_add

        else:
            v1 = random.choice(range(len(individual)))
            u1 = individual[v1]
            set_free_nodes = list(set(range(NB_G2_NODES)) - set(individual))
            u2 = random.choice(set_free_nodes)

            set_free_nodes_m = set(G2_NODES_MAPPING[list(set_free_nodes)])
            ea_add = update_Ea_change(v1, u1, u2, G1, G2, individual)
            ea_hat_add = update_Ea_hat(u1, u2, G2, set_free_nodes_m)

            individual[v1] = u2
            individual.Ea += ea_add
            individual.Ea_hat += ea_hat_add

            # if not(calculate_s3(individual.E1, individual.Ea, individual.Ea_hat) > individual.S3
            #         or random.random() < indpb):
            #     individual[v1] = u1
            #     individual.Ea -= ea_add
            #     individual.Ea_hat -= ea_hat_add
    return individual,


def cxOrdered(ind1, ind2):
    size = min(len(ind1), len(ind2))
    a, b = random.sample(range(size), 2)
    if a > b:
        a, b = b, a

    ind1_candidates = [j for j in ind2 if j not in set(ind1[a:b])]
    ind2_candidates = [j for j in ind1 if j not in set(ind2[a:b])]
    for i in range(size - (b - a)):
        index = (i + b) % size
        ind1[index] = ind1_candidates[i]
        ind2[index] = ind2_candidates[i]

    return ind1, ind2


def cxUniformPartialyMatched(ind1, ind2):
    size = min(len(ind1), len(ind2))
    max_size = max(len(ind1), len(ind2))
    p1, p2 = [-1] * len(G2.nodes), [-1] * len(G2.nodes)

    # Initialize the position of each indices in the individuals
    for i in range(size):
        p1[ind1[i]] = i
        p2[ind2[i]] = i
    # Choose crossover points
    std = int(size * 0.1)
    # k = random.randint(size / 2 - std, size / 2 + std)
    k = random.randint(0, size)
    cxpoint = random.sample(range(size), k)
    cxpoint2 = sorted(cxpoint)
    # Apply crossover between cx points
    for i in cxpoint2:
        # Keep track of the selected values
        temp1 = ind1[i]
        temp2 = ind2[i]
        # Swap the matched value
        ea_add1 = update_Ea_swap(i, p1[temp2], temp1, ind1[p1[temp2]], G1, G2, ind1)
        ea_add2 = update_Ea_swap(i, p2[temp1], temp2, ind2[p2[temp1]], G1, G2, ind2)
        ind1[i], ind1[p1[temp2]] = temp2, temp1
        ind2[i], ind2[p2[temp1]] = temp1, temp2

        ind1.Ea += ea_add1
        ind2.Ea += ea_add2
        # E1, Ea1, Ea_hat = get_align_edges(G1, G2, ind1)
        # E1, Ea2, Ea_hat = get_align_edges(G1, G2, ind2)
        # if Ea1 != ind1.Ea or ind2.Ea != Ea2:
        #     pass
        # Position bookkeeping
        p1[temp1], p1[temp2] = p1[temp2], p1[temp1]
        p2[temp1], p2[temp2] = p2[temp2], p2[temp1]

    # update_ind_parameters(ind1)
    # update_ind_parameters(ind2)
    return ind1, ind2


toolbox.register("evaluate", eval_alignment)
toolbox.register("mate", cxUniformPartialyMatched)
toolbox.register("mutate", mut_change_swap, indpb=0.3, change_prob=0.5)
toolbox.register("select", tools.selRoulette)


def main():
    pop = toolbox.population(n=30)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    # stats_size = tools.Statistics(lambda ind: get_nc(pd.Series(G2_NODES_MAPPING[ind].values, list(G1.nodes)), true_A))
    stats_size = tools.Statistics(lambda ind: 0)
    mstats = tools.MultiStatistics(s3=stats, nc=stats_size)
    mstats.register("Avg", np.mean)
    mstats.register("Std", np.std)
    mstats.register("Min", np.min)
    mstats.register("Max", np.max)

    algorithms.eaSimple(pop, toolbox, cxpb=0.7, mutpb=0.2, ngen=1000, stats=mstats,
                        halloffame=hof, verbose=True)

    return pop, stats, hof


if __name__ == "__main__":
    pop, stats, hof = main()
    print(hof)
    print(Counter(list(hof[0])).most_common(1))
    # series = G2_NODES_MAPPING[hof[0]]
    # series.index = list(G1.nodes)
    # print(series)
    print(get_nc(hof[0], true_A))
    print(calculate_s3(hof[0].E1, hof[0].Ea, hof[0].Ea_hat))
