# -*- coding: utf-8 -*-
import argparse
import root_numpy
import ROOT
import numpy as np
import pandas as pd
import os.path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml
import xgboost as xgb
from math import sqrt
from sklearn.cross_validation import StratifiedKFold
ROOT.gROOT.SetBatch(True)


def main(args):
    # Open config file
    f_config = open(args.config)
    global D_CONFIGS
    D_CONFIGS = yaml.load(f_config)
    ## Load files
    # List required branches (only once)
    s_load_branches = set()
    for run in list(D_CONFIGS.keys()):
        s_load_branches |= set(D_CONFIGS[run]['fit_vars'])
    s_load_branches.add('mjpipi')
    s_load_branches.add('mpipi')
    s_load_branches.add('scaledmass')
    l_load_branches = list(s_load_branches)
    # Load files into arrays
    print('*** Loading Data ***')
    a_mc_x = root_numpy.root2array(args.mc_x, treename = args.tree_name, branches = l_load_branches)
    a_mc_p = root_numpy.root2array(args.mc_p, treename = args.tree_name, branches = l_load_branches)
    a_side = root_numpy.root2array(args.side, treename = args.tree_name, branches = l_load_branches)

    ## Process for training
    print('*** Processing Data ***')
    # Convert to DataFrames
    df_mc_x = pd.DataFrame(a_mc_x)
    df_mc_p = pd.DataFrame(a_mc_p)
    df_side = pd.DataFrame(a_side)
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
        max_depth     = D_CONFIGS[run]['max_depth']
        learning_rate = D_CONFIGS[run]['learning_rate']
        n_estimators  = D_CONFIGS[run]['n_estimators']
        l_training_vars = list(D_CONFIGS[run]['fit_vars'])
        cv_train_accs = []
        cv_test_accs  = []
        ## Cross validate 10 times
        skf = StratifiedKFold(df_train['class'].as_matrix(), n_folds=10)
        for i, (train_index, test_index) in enumerate(skf):
            # Split into training and test sets
            df_train_cv = df_train.iloc[train_index].copy()
            df_test_cv  = df_train.iloc[test_index].copy()
            # Initialise model
            xgbc = xgb.XGBClassifier(max_depth=max_depth, learning_rate=learning_rate, n_estimators=n_estimators)
            # Train model
            print('*** Training model for run: {}, cv {} of {} ***'.format(run, i+1, len(skf)))
            xgbc.fit(df_train_cv[l_training_vars], df_train_cv['class'])
            # Predict probabilities
            df_train_cv['prob'] = xgbc.predict_proba(df_train_cv[l_training_vars])[:,1]  # Returns (1-p, p)
            df_test_cv['prob'] = xgbc.predict_proba(df_test_cv[l_training_vars])[:,1]
            ### Find optimal cut
            print('   *** Determining optimal cut ***')
            # Optimise the probability cut
            sig_effs_mcj = []  # record signal efficiencies - on MC X(3823) only
            sig_effs_mcx = []  # record signal efficiencies - on MC J(2S) only
            sig_effs_all = []  # record signal efficiencies
            bgr_rejs     = []  # record background rejections
            cut_scores   = []  # record cut optimisation metric
            # Determine cut metric for a range of cuts
            cuts = np.linspace(.0, 1., 200, endpoint=False)
            for prob_threshold in cuts:
                # Determine how many predictions are correct
                signal_efficiency_mcj = float(df_train_cv[(df_train_cv['prob'] > prob_threshold) & (df_train_cv['cat']   == 'mc_x')].shape[0]) / float(df_train_cv[df_train_cv['cat']   == 'mc_x'].shape[0])
                signal_efficiency_mcx = float(df_train_cv[(df_train_cv['prob'] > prob_threshold) & (df_train_cv['cat']   == 'mc_p')].shape[0]) / float(df_train_cv[df_train_cv['cat']   == 'mc_p'].shape[0])
                signal_efficiency_all = float(df_train_cv[(df_train_cv['prob'] > prob_threshold) & (df_train_cv['class'] == 1     )].shape[0]) / float(df_train_cv[df_train_cv['class'] == 1     ].shape[0])
                background_rejection  = float(df_train_cv[(df_train_cv['prob'] < prob_threshold) & (df_train_cv['class'] == 0     )].shape[0]) / float(df_train_cv[df_train_cv['class'] == 0     ].shape[0])
                # Store scores
                sig_effs_all.append(signal_efficiency_all)
                sig_effs_mcj.append(signal_efficiency_mcj)
                sig_effs_mcx.append(signal_efficiency_mcx)
                bgr_rejs    .append(background_rejection)
                # Optimize cut
                eff = signal_efficiency_all
                a   = 5. # expected significance
                # Background events, scaled to 40MeV window about B peak, considering only those in X(3823) region
                B  = df_train_cv[((df_train_cv['prob'] > prob_threshold)) & ((df_train_cv['scaledmass'] > 5400.) & (df_train_cv['scaledmass'] < 5450.)) & ((df_train_cv['mjpipi'] > 3773) & (df_train_cv['mjpipi'] < 3873))].shape[0] * .8            
                s0 = float(args.sig_yield)
                cut_scores.append((s0 * eff) / sqrt((s0 * eff) + B))
            # Find optimal cut
            cut_index = np.argmax(cut_scores)
            prob_threshold  = cuts[cut_index]
            print('   *** Optimised cut: {:.3f} ***'.format(prob_threshold))
            train_acc = float(df_train_cv[((df_train_cv['prob'] > prob_threshold) & (df_train_cv['class'] == 1)) | ((df_train_cv['prob'] < prob_threshold) & (df_train_cv['class'] == 0))].shape[0]) / float(df_train_cv.shape[0])
            test_acc  = float(df_test_cv [(( df_test_cv['prob'] > prob_threshold) & ( df_test_cv['class'] == 1)) | (( df_test_cv['prob'] < prob_threshold) & ( df_test_cv['class'] == 0))].shape[0]) / float(df_test_cv.shape[0])
            print('   *** Training accuracy: {:.3f} ***'.format(train_acc))
            print('   ***  Testing accuracy: {:.3f} ***'.format(test_acc))
            cv_train_accs.append(train_acc)
            cv_test_accs .append(test_acc)

        # Determine mean accuracies
        cv_train_acc = np.mean(cv_train_accs)
        cv_test_acc  = np.mean(cv_test_accs)
        # Print accuracies
        print('Training CV accuracy: {:.3f}'.format(cv_train_acc))
        print(' Testing CV accuracy: {:.3f}'.format(cv_test_acc))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "--add description--")
    parser.add_argument("mc_x"       ,       default = None                                                                     , help = "input MC X(3823) files")
    parser.add_argument("mc_p"       ,       default = None                                                                     , help = "input MC psi(2S) files")
    parser.add_argument("side"       ,       default = None                                                                     , help = "input root background files")
    parser.add_argument("sig_yield"  ,       default = None                                                                     , help = "specify the expected X(3823) signal yield for cut optimisation")
    parser.add_argument("--tree_name", "-t", default = "DecayTree"                                                              , help = "input tree name")
    parser.add_argument("--config"   , "-c", default = '/Disk/ds-sopa-group/PPE/users/s1305440/project_stuff/config/xgboost.yml', help =  "config file")
    args = parser.parse_args()
    main(args)
