# -*- coding: utf-8 -*-
import argparse
import root_numpy
import ROOT
import numpy as np
import pandas as pd
import os.path
import yaml
import xgboost as xgb
from shutil import copyfile
ROOT.gROOT.SetBatch(True)

def main(args):
    ### Prepare data for fitting

    # Ensure output directory exists
    out_path = args.out_path+'/'+args.info+'/'
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # Open config file
    f_config = open(args.config)
    global D_CONFIGS
    D_CONFIGS = yaml.load(f_config)

    ## Create a copy of the datafiles in the output directory
    # Destination is within data subfolder
    print('*** Copying files ***')
    print('   *** mc_x ***')
    dst_mc_x = out_path+'mc_x_proba.root'
    print('   *** mc_p ***')
    dst_mc_p = out_path+'mc_p_proba.root'
    print('   *** side ***')    
    dst_side = out_path+'side_proba.root'
    print('   *** data ***')
    dst_data = out_path+'data_proba.root'
    # Move files here
    copyfile(args.mc_x, dst_mc_x)
    copyfile(args.mc_p, dst_mc_p)
    copyfile(args.side, dst_side)
    copyfile(args.data, dst_data)

    ## Load files
    # List required branches (only once)
    s_load_branches = set()
    for run in list(D_CONFIGS.keys()):
        s_load_branches |= set(D_CONFIGS[run]['fit_vars'])
    s_load_branches.add('mjpipi')
    s_load_branches.add('mpipi')
    l_load_branches = list(s_load_branches)
    # Load files into arrays
    print('*** Loading Data ***')
    a_mc_x = root_numpy.root2array(dst_mc_x, treename = args.tree_name, branches = l_load_branches)
    a_mc_p = root_numpy.root2array(dst_mc_p, treename = args.tree_name, branches = l_load_branches)
    a_side = root_numpy.root2array(dst_side, treename = args.tree_name, branches = l_load_branches)
    a_data = root_numpy.root2array(dst_data, treename = args.tree_name, branches = l_load_branches)

    ## Process for training
    print('*** Processing Data ***')
    # Convert to DataFrames
    df_mc_x = pd.DataFrame(a_mc_x)
    df_mc_p = pd.DataFrame(a_mc_p)
    df_side = pd.DataFrame(a_side)
    df_data = pd.DataFrame(a_data)
    # Add categoriastion
    df_mc_x['cat'] = 'mc_x'
    df_mc_p['cat'] = 'mc_p'
    df_side['cat'] = 'side'
    # Add target
    df_mc_x['class'] = 1
    df_mc_p['class'] = 1
    df_side['class'] = 0
    # Combine training sets
    df_train = pd.concat([df_mc_x, df_mc_p, df_side])

    ### Predict probabilities
    ## Train model for each set of parameters
    for run in list(D_CONFIGS.keys()):
        # Initialise model
        max_depth     = D_CONFIGS[run]['max_depth']
        learning_rate = D_CONFIGS[run]['learning_rate']
        n_estimators  = D_CONFIGS[run]['n_estimators']
        xgbc = xgb.XGBClassifier(max_depth=max_depth, learning_rate=learning_rate, n_estimators=n_estimators)
        # Train model
        print('*** Training model for run: {} ***'.format(run))
        l_training_vars = list(D_CONFIGS[run]['fit_vars'])
        xgbc.fit(df_train[l_training_vars], df_train['class'])
        # Predict probabilities
        print('*** Predicting training probabilities ***')
        df_mc_x['prob_'+run] = xgbc.predict_proba(df_train[df_train['cat'] == 'mc_x'][l_training_vars])[:,1]  # Returns (1-p, p)
        df_mc_p['prob_'+run] = xgbc.predict_proba(df_train[df_train['cat'] == 'mc_p'][l_training_vars])[:,1]
        df_side['prob_'+run] = xgbc.predict_proba(df_train[df_train['cat'] == 'side'][l_training_vars])[:,1]
        print('*** Predicting data probabilities ***')
        df_data['prob_'+run] = xgbc.predict_proba(df_data[l_training_vars])[:,1]
        # Write probabilities to root files
        print('*** Writing to files ***')
        a_mc_x_prob = np.array(df_mc_x['prob_'+run], dtype=[('prob_'+run, 'f8')])
        root_numpy.array2root(a_mc_x_prob, dst_mc_x, args.tree_name)
        a_mc_p_prob = np.array(df_mc_p['prob_'+run], dtype=[('prob_'+run, 'f8')])
        root_numpy.array2root(a_mc_p_prob, dst_mc_p, args.tree_name)
        a_side_prob = np.array(df_side['prob_'+run], dtype=[('prob_'+run, 'f8')])
        root_numpy.array2root(a_side_prob, dst_side, args.tree_name)
        a_data_prob = np.array(df_data['prob_'+run], dtype=[('prob_'+run, 'f8')])
        root_numpy.array2root(a_data_prob, dst_data, args.tree_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "--add description--")
    parser.add_argument("mc_x"       ,       default = None                , help = "input MC X(3823) file")
    parser.add_argument("mc_p"       ,       default = None                , help = "input MC psi(2S) file")
    parser.add_argument("side"       ,       default = None                , help = "input sideband file")
    parser.add_argument("data"       ,       default = None                , help = "input data file")
    parser.add_argument("--out_path" , "-o", default = 'data'              , help = "directory for storing results")
    parser.add_argument("--config"   , "-c", default = 'config/xgboost.yml', help =  "config file")
    parser.add_argument("--tree_name", "-t", default = "DecayTree"         , help = "input tree name")
    parser.add_argument("--info"     , "-i", default = "probs"             , help = "string specifying information about the run")
    args = parser.parse_args()
    main(args)