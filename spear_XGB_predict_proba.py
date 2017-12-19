# -*- coding: utf-8 -*-
import argparse
import root_numpy
import ROOT
import numpy as np
import pandas as pd
import os.path
import yaml
from xgboost import XGBClassifier
from sklearn.cross_validation import StratifiedShuffleSplit
from sklearn.metrics import log_loss
from shutil import copyfile
ROOT.gROOT.SetBatch(True)

### Command line arguments - look at how to combine with Spearmint
run_info = "spearmint_prep"
out_path = "data/" + run_info + "/"
loc_mc_x = "data/mcx_cut.root"
loc_mc_p = "data/mcpsi_cut.root"
loc_side = "data/data_sideband_cut.root"
loc_data = "data/data_cut.root"

def predict_probability(max_depth = 6, learning_rate = 0.1, n_estimators = 600):
    ### Prepare data for fitting

    # Ensure output directory exists
    if not os.path.exists(out_path):
        os.makedirs(out_path)

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
    copyfile(loc_mc_x, dst_mc_x)
    copyfile(loc_mc_p, dst_mc_p)
    copyfile(loc_side, dst_side)
    copyfile(loc_data, dst_data)

    ## Load files
    l_fit_vars = ['logDIRA', 'log_bplus_IPCHI2_OWNPV', 'bplus_LOKI_DTF_CHI2NDOF', 'log_bplus_FDCHI2_OWNPV', 'bplus_ETA', 
                  'log_1_IPCHI2_OWNPV', 'log_2_IPCHI2_OWNPV', 'log_3_IPCHI2_OWNPV', 'log_4_IPCHI2_OWNPV', 'log_5_IPCHI2_OWNPV', 
                  'mu_PT_max', 'mu_PT_min']
    # Load files into arrays
    print('*** Loading Data ***')
    a_mc_x = root_numpy.root2array(dst_mc_x, treename = 'DecayTree', branches = l_fit_vars)
    a_mc_p = root_numpy.root2array(dst_mc_p, treename = 'DecayTree', branches = l_fit_vars)
    a_side = root_numpy.root2array(dst_side, treename = 'DecayTree', branches = l_fit_vars)
    a_data = root_numpy.root2array(dst_data, treename = 'DecayTree', branches = l_fit_vars)

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
    # Split into training and validation sets
    sss = StratifiedShuffleSplit(df_train['class'].as_matrix(), n_iter=1, test_size=0.2)
    for trn_index, val_index in sss:
        df_trn = df_train.iloc[trn_index]
        df_val = df_train.iloc[val_index]

    ### Predict probabilities
    xgbc = XGBClassifier(max_depth=max_depth, learning_rate=learning_rate, n_estimators=n_estimators)
    # Train model
    print('*** Training model ***')    
    xgbc.fit(df_trn[l_fit_vars], df_trn['class'])
    # Predict probabilities
    print('*** Predicting training probabilities ***')
    df_mc_x['prob'] = xgbc.predict_proba(df_train[df_train['cat'] == 'mc_x'][l_fit_vars])[:,1]  # Returns (1-p, p)
    df_mc_p['prob'] = xgbc.predict_proba(df_train[df_train['cat'] == 'mc_p'][l_fit_vars])[:,1]
    df_side['prob'] = xgbc.predict_proba(df_train[df_train['cat'] == 'side'][l_fit_vars])[:,1]
    print('*** Predicting data probabilities ***')
    df_data['prob'] = xgbc.predict_proba(df_data[l_fit_vars])[:,1]
    # Write probabilities to root files
    print('*** Writing to files ***')
    a_mc_x_prob = np.array(df_mc_x['prob'], dtype=[('prob', 'f8')])
    root_numpy.array2root(a_mc_x_prob, dst_mc_x, 'DecayTree')
    a_mc_p_prob = np.array(df_mc_p['prob'], dtype=[('prob', 'f8')])
    root_numpy.array2root(a_mc_p_prob, dst_mc_p, 'DecayTree')
    a_side_prob = np.array(df_side['prob'], dtype=[('prob', 'f8')])
    root_numpy.array2root(a_side_prob, dst_side, 'DecayTree')
    a_data_prob = np.array(df_data['prob'], dtype=[('prob', 'f8')])
    root_numpy.array2root(a_data_prob, dst_data, 'DecayTree')

    # Determine log-loss for training and validation sets
    ll_trn = log_loss(df_trn['class'].as_matrix(), xgbc.predict_proba(df_trn[l_fit_vars])[:,1])
    ll_val = log_loss(df_val['class'].as_matrix(), xgbc.predict_proba(df_val[l_fit_vars])[:,1])
    # Print and return
    print("*** Training set log-loss: {:.3f} ***".format(ll_trn))
    print("*** validation set log-loss: {:.3f} ***".format(ll_val))
    return ll_val


# Remove this line when using spearmint
predict_probability()

# Main definition in-line with Spearmint requirements
def main(job_id, params):
    max_depth     = params['max_depth'][0]
    learning_rate = params['learning_rate'][0]
    n_estimators  = params['n_estimators'][0]
    res = predict_probability(max_depth, learning_rate, n_estimators)
    print("Optimising XGBoost Classifier Hyperparameters")
    print('\t(max_depth: {}, learning_rate: {.3f}, n_estimators: {}) -> ll_val: {:.3f}'.format(max_depth, learning_rate, n_estimators, res))
    return predict_probability(max_depth, learning_rate, n_estimators)