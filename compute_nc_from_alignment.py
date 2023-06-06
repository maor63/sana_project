import networkx as nx
import random
import pandas as pd
import os


def get_nc(alignment, true_alignment):
    correct = sum(1 for v in alignment if alignment[v] == true_alignment[v])
    return correct * 1.0 / len(alignment)


true_alignment_rw = {'O{}'.format(v): 'Rw{}'.format(v) for v in xrange(20)}
true_alignment_add_rw = {'O{}'.format(v): 'AddRw{}'.format(v) for v in xrange(20)}
true_alignment_delete_rw = {'O{}'.format(v): 'DelRw{}'.format(v) for v in xrange(20)}

base_dir = 'sana_results_ER20_ER60_compose_S3_ACC__F1_ICS_EC_20_min_iter_10'
run_times_dirs_names = os.listdir(base_dir)
rows = []
for run_time in run_times_dirs_names:
    run_time_path = os.path.join(base_dir, run_time)
    objective_functions = os.listdir(run_time_path)
    for objective_function in objective_functions:
        objective_function_path = os.path.join(run_time_path, objective_function)
        experiments = os.listdir(objective_function_path)
        for experiment in experiments:
            experiment_align_path = os.path.join(objective_function_path, experiment,
                                                 'sana_' + objective_function.lower() + '.align')
            # output_file = open(experiment_align_path)
            alignment_tsv = pd.read_csv(experiment_align_path, delimiter='\t')
            alignment_dict = dict(alignment_tsv.to_records(index=False))
            nc_rw = get_nc(alignment_dict, true_alignment_rw)
            nc_add_rw = get_nc(alignment_dict, true_alignment_add_rw)
            nc_delete_rw = get_nc(alignment_dict, true_alignment_delete_rw)
            rows.append([experiment, objective_function, nc_rw, nc_add_rw, nc_delete_rw])

nc_results_df = pd.DataFrame(rows, columns=['experiment', 'objective', 'rw', 'nc_add_rw', 'nc_delete_rw'])
nc_results_df.to_csv('sana_results_ER20_ER60_compose_S3_ACC__F1_ICS_EC_20_min_iter_10_nc_results.csv')
