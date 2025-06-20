import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def combine_assignments(data_dir):
    '''
    Combines assignment files for methods that have been performed on gRNA subsets
    
    Args:
        data_dir (str): directory containing the assignment output
    Returns:
        None
    '''
    # check which perturbation files are in the directory
    all_files = os.listdir(data_dir)
    csv_files = [file for file in all_files if file.startswith('perturbations_') and file.endswith('.csv')]

    perturbations = pd.DataFrame()
    
    # read in data for all gRNA subsets
    for file in csv_files:
        perturbations_g = pd.read_csv(data_dir + file)
        perturbations = pd.concat([perturbations, perturbations_g], ignore_index = True)
    
    # remove duplicates (in case there were overlapping indices)
    perturbations = perturbations.drop_duplicates()
    
    # save file containing all gRNAs
    perturbations.to_csv(data_dir + "perturbations.csv")
    
    
def load_assignments(gc_dict, data_dir):
    '''
    Loads and combines assignment files of specified methods
    
    Args:
        gc_dict (dict): dictionary with method directory and threshold for each assignment
        data_dir (str): directory containing all assignment output folders
        
    Returns:
        A pd DataFrame containing the gRNA-cell assignments of each method
    '''
    perturbations = pd.DataFrame()
    for name, method_params in gc_dict.items():
        method_dir, t = method_params

        if t == None: 
            assignment = pd.read_csv(os.path.join(data_dir, method_dir, 'assignments.csv'))[['cell', 'gRNA']]
        else:
            assignment = pd.read_csv(data_dir + method_dir + '/assignments_t' + str(t) + '.csv')[['cell', 'gRNA']]
        assignment['method'] = name

        # Combine with all others
        perturbations = pd.concat([perturbations, assignment])
    return perturbations


def calculate_jaccard(pert_dict):
    '''
    Calculates pairwise Jaccard index
    
    Args:
        pert_dict (dictionary): dictionary which contains assigned cell_gRNAs for each method
        
    Returns:
        A pd DataFrame with the pairwise similarity scores
    '''
    matrix = np.zeros((len(pert_dict), len(pert_dict)))
    results = pd.DataFrame(matrix, index=pert_dict.keys(), columns=pert_dict.keys())

    for i in list(pert_dict.keys()):
        for j in list(pert_dict.keys()):
            set1 = set(pert_dict[i])
            set2 = set(pert_dict[j])    
            intersection_size = len(set1.intersection(set2))
            union_size = len(set1.union(set2))
            similarity = intersection_size / union_size
            results.loc[i, j] = similarity
    return results

