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

def fit_doubleCB(a_train, a_data, out_path, s_info=''):
    # Initialise dictionary for storing fit info
    d_fit_info = {}
    # Format arrays
    a_fit_data = a_data.astype(dtype=[('b_mass', np.float)])
    a_fit_mc   = a_train.astype(dtype=[('b_mass_mc', np.float)])
    # Estimate sig/bkg yields
    max_yield_est = a_fit_data.shape[0]                         # max possible signal yield
    bkg_yield_est = len(a_fit_data[a_fit_data > 5350]) * 4.     # estimate background from sideband
    # Create tree for fitting
    t_fit_data = root_numpy.array2tree(a_fit_data)
    t_fit_mc   = root_numpy.array2tree(a_fit_mc)
    
    ## Monte Carlo
    # Inialise parameters for fit
    b_mass_mc     = ROOT.RooRealVar("b_mass_mc"   , "B mass MC [MeV]", 5220. , 5350.)
    mean_mc       = ROOT.RooRealVar("mean_mc"     , "mean_mc"        , 5279. , 5195.  , 5400.)
    sig_1_mc      = ROOT.RooRealVar("sig_1_mc"    , "sig_1_mc"       ,    5. ,     .1 ,   15.)
    alpha_1_mc    = ROOT.RooRealVar("alpha_1_mc"  , "alpha_1_mc"     ,    5. ,     .1 ,   10.)
    n_1_mc        = ROOT.RooRealVar("n_1_mc"      , "n_1_mc"         ,    2. ,    0.  ,   15.)
    r_s1_s2_mc    = ROOT.RooRealVar("r_s1_s2_mc"  , "r_s1_s2_mc"     ,     .5,     .1 ,   10.)
    sig_2_mc      = ROOT.RooFormulaVar("sig_2_mc" , "sig_1_mc*r_s1_s2_mc", ROOT.RooArgList(sig_1_mc, r_s1_s2_mc))
    alpha_2_mc    = ROOT.RooRealVar("alpha_2_mc"  , "alpha_2_mc"     ,   -5. ,  -10.  ,    -.1)
    n_2_mc        = ROOT.RooRealVar("n_2_mc"      , "n_2_mc"         ,    2. ,    0.  ,   15.)
    r_cb1_cb2_mc  = ROOT.RooRealVar("r_cb1_cb2_mc",  "r_cb1_cb2_mc"  ,     .5,     .01,    1.)
    # Initialise fit model
    cb_1_mc      = ROOT.RooCBShape("cb_1_mc"     , "cb_1_mc"     , b_mass_mc, mean_mc, sig_1_mc, alpha_1_mc, n_1_mc)
    cb_2_mc      = ROOT.RooCBShape("cb_2_mc"     , "cb_2_mc"     , b_mass_mc, mean_mc, sig_2_mc, alpha_2_mc, n_2_mc)
    model_sig_mc = ROOT. RooAddPdf("model_sig_mc", "model_sig_mc", ROOT.RooArgList(cb_1_mc, cb_2_mc), ROOT.RooArgList(r_cb1_cb2_mc))
    # Initialise dataset
    dataset_mc = ROOT.RooDataSet("dataset_mc","dataset from tree", t_fit_mc, ROOT.RooArgSet(b_mass_mc))
    # Perform fit
    mean_mc     .setConstant(ROOT.kFALSE)
    sig_1_mc    .setConstant(ROOT.kFALSE)
    alpha_1_mc  .setConstant(ROOT.kFALSE)
    n_1_mc      .setConstant(ROOT.kFALSE)
    r_s1_s2_mc  .setConstant(ROOT.kFALSE)
    alpha_2_mc  .setConstant(ROOT.kFALSE)
    n_2_mc      .setConstant(ROOT.kFALSE)
    r_cb1_cb2_mc.setConstant(ROOT.kFALSE)
    model_sig_mc.fitTo(dataset_mc, ROOT.RooFit.Range(5220., 5350.))
    # Store fitted values
    f_mean_mc       = ROOT.RooRealVar("f_mean_mc"   , "f_mean_mc"   , mean_mc   .getValV())
    f_sig_1_mc      = ROOT.RooRealVar("f_sig_1_mc"  , "f_sig_1_mc"  , sig_1_mc  .getValV())
    f_alpha_1_mc    = ROOT.RooRealVar("f_alpha_1_mc", "f_alpha_1_mc", alpha_1_mc.getValV())
    f_n_1_mc        = ROOT.RooRealVar("f_n_1_mc"    , "f_n_1_mc"    , n_1_mc    .getValV())
    f_sig_2_mc      = ROOT.RooRealVar("f_sig_2_mc"  , "f_sig_2_mc"  , sig_2_mc  .getValV())
    f_alpha_2_mc    = ROOT.RooRealVar("f_alpha_2_mc", "f_alpha_2_mc", alpha_2_mc.getValV())
    f_n_2_mc        = ROOT.RooRealVar("f_n_2_mc"    , "f_n_2_mc"    , n_1_mc    .getValV())
    # Store fitted models
    f_cb_1_mc      = ROOT.RooCBShape("cb_1_mc"     , "cb_1_mc"     , b_mass_mc, f_mean_mc, f_sig_1_mc, f_alpha_1_mc, f_n_1_mc)
    f_cb_2_mc      = ROOT.RooCBShape("cb_2_mc"     , "cb_2_mc"     , b_mass_mc, f_mean_mc, f_sig_2_mc, f_alpha_2_mc, f_n_2_mc)
    f_model_sig_mc = ROOT. RooAddPdf("model_sig_mc", "model_sig_mc", ROOT.RooArgList(f_cb_1_mc, f_cb_2_mc), ROOT.RooArgList(r_cb1_cb2_mc))
    ## Plot
    # Frame for fit
    frame_mc = b_mass_mc.frame()
    # Frame for pulls
    frame_mc_pull = b_mass_mc.frame()
    # Add data and fit to frame
    dataset_mc.plotOn(frame_mc)
    f_model_sig_mc.plotOn(frame_mc)
    # Plot on split canvas
    c = ROOT.TCanvas("X3872MC", "X3872MC", 400, 500)
    c.Divide(1, 2, 0, 0)
    # Plot data and fit
    c.cd(2)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetPad(.01,.01,.95,.77)
    frame_mc.SetTitle("Fitted Monte-Carlo Bmass")
    frame_mc.SetMaximum(frame_mc.GetMaximum()*1.1)
    frame_mc.GetYaxis().SetTitleOffset(1.6)
    frame_mc.Draw()
    # Plot pulls
    c.cd(1)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetPad(.01,.76,.95,.97)
    # Determine pulls and format
    h_pull_mc = frame_mc.pullHist()
    h_pull_mc.SetFillColor(15)
    h_pull_mc.SetFillStyle(3144)
    # Add pulls to frame
    frame_mc_pull.addPlotable(h_pull_mc,'L3')
    frame_mc_pull.GetYaxis().SetNdivisions(505)
    frame_mc_pull.GetYaxis().SetLabelSize(0.20)
    frame_mc_pull.SetTitle("")
    frame_mc_pull.Draw()
    # Save canvas
    c.SaveAs(out_path+'/'+s_info+'_FittedMassDistribution_MonteCarlo.pdf')
    # Store fit variables
    d_fit_info['mc_mean']     = mean_mc.getValV()
    d_fit_info['mc_mean_err'] = mean_mc.getError()
    d_fit_info['mc_sig1']     = sig_1_mc.getValV() 
    d_fit_info['mc_sig1_err'] = sig_1_mc.getError()
    d_fit_info['mc_r_s1_s2']     = r_s1_s2_mc.getValV() 
    d_fit_info['mc_r_s1_s2_err'] = r_s1_s2_mc.getError()
    d_fit_info['mc_alpha1']     = alpha_1_mc.getValV()
    d_fit_info['mc_alpha1_err'] = alpha_1_mc.getError()
    d_fit_info['mc_alpha2']     = alpha_2_mc.getValV()
    d_fit_info['mc_alpha2_err'] = alpha_2_mc.getError()
    d_fit_info['mc_n1']     = n_1_mc.getValV()
    d_fit_info['mc_n1_err'] = n_1_mc.getError()
    d_fit_info['mc_n2']     = n_2_mc.getValV()
    d_fit_info['mc_n2_err'] = n_2_mc.getError()
    d_fit_info['mc_r_cb1_cb2']     = r_cb1_cb2_mc.getValV()
    d_fit_info['mc_r_cb1_cb2_err'] = r_cb1_cb2_mc.getError()
    d_fit_info['mc_fit_chi2'] = frame_mc.chiSquare()

    ## Data
    b_mass    = ROOT.RooRealVar("b_mass"   , "B mass [MeV]" , 5220., 5380.)
    mean      = ROOT.RooRealVar("mean"     , "mean"         , mean_mc     .getValV(), 5220., 5380.)
    sig_1     = ROOT.RooRealVar("sig_1"    , "sig_1"        , sig_1_mc    .getValV(), .1, 15.)
    alpha_1   = ROOT.RooRealVar("alpha_1"  , "alpha_1"      , alpha_1_mc  .getValV(), .1, 10.)
    n_1       = ROOT.RooRealVar("n_1"      , "n_1"          , n_1_mc      .getValV(), 0., 15.)
    r_s1_s2   = ROOT.RooRealVar("r_s1_s2"  , "r_s1_s2"      , r_s1_s2_mc  .getValV(), 0.01, 10)
    sig_2     = ROOT.RooFormulaVar("sig_2" , "sig_1*r_s1_s2", ROOT.RooArgList(sig_1, r_s1_s2))
    alpha_2   = ROOT.RooRealVar("alpha_2"  , "alpha_2"      , alpha_2_mc  .getValV(), -10., -.1)
    n_2       = ROOT.RooRealVar("n_2"      , "n_2"          , n_2_mc      .getValV(), 0., 15.)
    r_cb1_cb2 = ROOT.RooRealVar("r_cb1_cb2",  "r_cb1_cb2"   , r_cb1_cb2_mc.getValV(), 0.01  , 1.);
    sig_yield = ROOT.RooRealVar("sig_yield", "sig_yield"    , 0, (max_yield_est - bkg_yield_est)*1.5)
    # Background
    exp_c     = ROOT.RooRealVar("exp_c"        , "exp_c"    , 0., -.02, .02)
    bgr_yield = ROOT.RooRealVar("bgr_yield"    , "bgr_yield", bkg_yield_est*.5, max_yield_est)
    # Initialise fit model
    cb_1      = ROOT.RooCBShape("cb_1"     , "cb_1"     , b_mass, mean, sig_1, alpha_1, n_1)
    cb_2      = ROOT.RooCBShape("cb_2"     , "cb_2"     , b_mass, mean, sig_2, alpha_2, n_2)
    exp_bg    = ROOT.RooExponential("exp_bg", "exp_bg", b_mass, exp_c)
    model_sig = ROOT.RooAddPdf("model_sig", "model_sig", ROOT.RooArgList(cb_1, cb_2)       , ROOT.RooArgList(r_cb1_cb2))  
    model_tot = ROOT.RooAddPdf("model_tot", "model_tot", ROOT.RooArgList(model_sig, exp_bg), ROOT.RooArgList(sig_yield, bgr_yield))
    # Initialise dataset
    dataset = ROOT.RooDataSet("dataset","dataset from tree", t_fit_data, ROOT.RooArgSet(b_mass))
    # Perform fit - kTRUE vars determined from MC fit
    mean     .setConstant(ROOT.kFALSE)
    sig_1    .setConstant(ROOT.kFALSE)
    alpha_1  .setConstant(ROOT.kTRUE)
    n_1      .setConstant(ROOT.kTRUE)
    r_s1_s2  .setConstant(ROOT.kTRUE)
    alpha_2  .setConstant(ROOT.kTRUE)
    n_2      .setConstant(ROOT.kTRUE)
    r_cb1_cb2.setConstant(ROOT.kTRUE) 
    sig_yield.setConstant(ROOT.kFALSE)
    exp_c    .setConstant(ROOT.kFALSE)
    bgr_yield.setConstant(ROOT.kFALSE)
    model_tot.fitTo(dataset, ROOT.RooFit.Range(5220., 5380.))
    # Store fitted values
    f_mean      = ROOT.RooRealVar("f_mean"     , "f_mean"     , mean     .getValV())
    f_sig_1     = ROOT.RooRealVar("f_sig_1"    , "f_sig_1"    , sig_1    .getValV())
    f_alpha_1   = ROOT.RooRealVar("f_alpha_1"  , "f_alpha_1"  , alpha_1  .getValV())
    f_n_1       = ROOT.RooRealVar("f_n_1"      , "f_n_1"      , n_1      .getValV())
    f_sig_2     = ROOT.RooRealVar("f_sig_2"    , "f_sig_2"    , sig_2    .getValV())
    f_alpha_2   = ROOT.RooRealVar("f_alpha_2"  , "f_alpha_2"  , alpha_2  .getValV())
    f_n_2       = ROOT.RooRealVar("f_n_2"      , "f_n_2"      , n_1      .getValV())
    f_exp_c     = ROOT.RooRealVar("f_exp_c"    , "f_exp_c"    , exp_c    .getValV())
    f_sig_yield = ROOT.RooRealVar("f_sig_yield", "f_sig_yield", sig_yield.getValV())
    f_bgr_yield = ROOT.RooRealVar("f_bgr_yield", "f_bgr_yield", bgr_yield.getValV())
    # Store fitted models
    f_cb_1      = ROOT.    RooCBShape("cb_1"     , "cb_1"     , b_mass, f_mean, f_sig_1, f_alpha_1, f_n_1)
    f_cb_2      = ROOT.    RooCBShape("cb_2"     , "cb_2"     , b_mass, f_mean, f_sig_2, f_alpha_2, f_n_2)
    f_exp_bg    = ROOT.RooExponential("exp_bg"   , "exp_bg"   , b_mass, f_exp_c)
    f_model_sig = ROOT.     RooAddPdf("model_sig", "model_sig", ROOT.RooArgList(f_cb_1, f_cb_2), ROOT.RooArgList(r_cb1_cb2))
    f_model_tot = ROOT.     RooAddPdf("model_tot", "model_tot", ROOT.RooArgList(f_model_sig, f_exp_bg), ROOT.RooArgList(f_sig_yield, f_bgr_yield))
    # Plot
    # Frame for fit
    frame = b_mass.frame()
    # Frame for pulls
    frame_pull = b_mass.frame()
    # Add data and fit to frame
    dataset.plotOn(frame)
    f_model_tot.plotOn(frame)
    # Plot on split canvas
    c = ROOT.TCanvas("data_fit", "data_fit", 400, 500)
    c.Divide(1, 2, 0, 0)
    # Plot data and fit
    c.cd(2)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetPad(.01,.01,.95,.77)
    frame.SetTitle("Fitted Data Bmass")
    frame.SetMaximum(frame.GetMaximum()*1.1)
    frame.GetYaxis().SetTitleOffset(1.6)
    frame.Draw()
    # Plot pulls
    c.cd(1)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetPad(.01,.76,.95,.97)
    # Determine pulls and format
    h_pull = frame.pullHist()
    h_pull.SetFillColor(15)
    h_pull.SetFillStyle(3144)
    # Add pulls to frame
    frame_pull.addPlotable(h_pull,'L3')
    frame_pull.GetYaxis().SetNdivisions(505)
    frame_pull.GetYaxis().SetLabelSize(0.20)
    frame_pull.SetTitle("")
    frame_pull.Draw()
    # Save plot
    c.SaveAs(out_path+'/'+s_info+'_FittedMassDistribution_Data.pdf')
    # Store fit variables
    d_fit_info['data_mean']     = mean.getValV()
    d_fit_info['data_mean_err'] = mean.getError()
    d_fit_info['data_sig1']     = sig_1.getValV() 
    d_fit_info['data_sig1_err'] = sig_1.getError()
    d_fit_info['data_r_s1_s2']     = r_s1_s2.getValV() 
    d_fit_info['data_r_s1_s2_err'] = r_s1_s2.getError()
    d_fit_info['data_alpha1']     = alpha_1.getValV()
    d_fit_info['data_alpha1_err'] = alpha_1.getError()
    d_fit_info['data_alpha2']     = alpha_2.getValV()
    d_fit_info['data_alpha2_err'] = alpha_2.getError()
    d_fit_info['data_n1']     = n_1.getValV()
    d_fit_info['data_n1_err'] = n_1.getError()
    d_fit_info['data_n2']     = n_2.getValV()
    d_fit_info['data_n2_err'] = n_2.getError()
    d_fit_info['data_r_cb1_cb2']     = r_cb1_cb2.getValV()
    d_fit_info['data_r_cb1_cb2_err'] = r_cb1_cb2.getError()
    d_fit_info['data_expc']     = exp_c.getValV()
    d_fit_info['data_expc_err'] = exp_c.getError()
    d_fit_info['data_sig_yield']     = sig_yield.getValV()
    d_fit_info['data_sig_yield_err'] = sig_yield.getError()
    d_fit_info['data_bgr_yield']     = bgr_yield.getValV()
    d_fit_info['data_bgr_yield_err'] = bgr_yield.getError()
    d_fit_info['data_fit_chi2'] = frame.chiSquare()

    # Return fit info dictionary
    return d_fit_info

def train_model(max_depth, learning_rate, n_estimators, fit_vars, df_train_in, df_test_in=None, s0=None, bkg_cut=False, determine_probas=False):
    # Copy dataframes so originals are unchanged
    df_train = df_train_in.copy()
    if determine_probas:
        print('      *** Determining probabilities ***')
        # Initialise classifier
        xgbc = xgb.XGBClassifier(max_depth=max_depth, learning_rate=learning_rate, n_estimators=n_estimators)
        # Reclassify as binary for training
        df_train['training_class'] = df_train['class'] != 0
        # Train model
        print('      *** Training model ***')
        xgbc.fit(df_train[fit_vars], df_train['training_class'])
        # Assign a probability to each training event
        print('      *** Fitting model ***')
        df_train['probability'] = xgbc.predict_proba(df_train[fit_vars])[:,1] # returns (1-p, p)
    # Optimise the probability cut
    sig_effs_mcj = []  # record signal efficiencies - on MC X(3823) only
    sig_effs_mcx = []  # record signal efficiencies - on MC J(2S) only
    sig_effs_all = []  # record signal efficiencies
    bgr_rejs     = []  # record background rejections
    cut_scores   = []  # record cut optimisation metric
    # Determine cut metric for a range of cuts
    cuts = np.linspace(.0, .99, 200)
    for prob_threshold in cuts:
        # Determine how many predictions are correct
        signal_efficiency_mcj = float(df_train[(df_train['probability'] > prob_threshold) & (df_train['class'] == 1)].shape[0]) / float(df_train[df_train['class'] == 1].shape[0])
        signal_efficiency_mcx = float(df_train[(df_train['probability'] > prob_threshold) & (df_train['class'] == 2)].shape[0]) / float(df_train[df_train['class'] == 2].shape[0])
        signal_efficiency_all = float(df_train[(df_train['probability'] > prob_threshold) & (df_train['class'] != 0)].shape[0]) / float(df_train[df_train['class'] != 0].shape[0])
        background_rejection  = float(df_train[(df_train['probability'] < prob_threshold) & (df_train['class'] == 0)].shape[0]) / float(df_train[df_train['class'] == 0].shape[0])
        # Store scores
        sig_effs_all.append(signal_efficiency_all)
        sig_effs_mcj.append(signal_efficiency_mcj)
        sig_effs_mcx.append(signal_efficiency_mcx)
        bgr_rejs    .append(background_rejection)
        # Optimize cut
        eff = signal_efficiency_all
        a   = 5. # expected significance
        # Background events, scaled to 40MeV window about B peak, considering only those in X(3823) region
        B = df_train[((df_train['probability'] > prob_threshold)) & ((df_train['scaledmass'] > 5350.) & (df_train['scaledmass'] < 5400.)) & ((df_train['mjpipi'] > 3773) & (df_train['mjpipi'] < 3873))].shape[0] * .8
        if s0 is not None:
            cut_scores.append((s0 * eff) / sqrt((s0 * eff) + B))
        else:
            cut_scores.append(eff / ((a / 2) + sqrt(B)))
    # Find optimal cut
    if bkg_cut:
        cut_index = np.argmax(np.array(bgr_rejs)>.99)
        print("Background used: {:.3f}".format(bgr_rejs[cut_index]))
        prob_threshold = cuts[cut_index]
    else:
        cut_index = np.argmax(cut_scores)
        prob_threshold  = cuts[cut_index]
    # Cross validate if test set supplied
    cv_accuracy = None
    if df_test_in is not None:
        # Apply model to test set
        df_test  = df_test_in.copy()
        if determine_probas:
            df_test['probability'] = xgbc.predict_proba(df_test[fit_vars])[:,1]
        # Determine accuracy
        cv_accuracy = float(df_test[((df_test['probability'] > prob_threshold) & (df_test['class'] != 0)) | ((df_test['probability'] < prob_threshold) & (df_test['class'] == 0))].shape[0]) / float(df_test.shape[0])
    # Store training objects in output dictionary
    d_training = {'sig_effs_all': sig_effs_all,
                  'sig_effs_mcj': sig_effs_mcj,
                  'sig_effs_mcx': sig_effs_mcx,
                  'bgr_rejs'    : bgr_rejs,
                  'cuts'        : cuts,
                  'cut_scores'  : cut_scores,
                  'cut_index'   : cut_index,
                  'best_cut'    : prob_threshold,
                  'cv_accuracy' : cv_accuracy}
    if determine_probas:
        d_training['xgb_trained'] = xgbc 
        d_training['probability'] = df_train['probability']
    # Return dictionary holding everything we've learned
    return d_training

def main(args):
    ### Prepare data for fitting
    # Ensure output directory exists
    out_path = args.out_path+'/'+args.info
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # Open config file
    f_config = open(args.config)
    global D_CONFIGS
    D_CONFIGS = yaml.load(f_config)

    ### Estimate signal yield
    s0 = None
    if args.find_s0:
        a_mc_x = root_numpy.root2array(args.mc_x, treename = args.tree_name, branches='scaledmass')
        a_data_x = root_numpy.root2array(args.data, treename = args.tree_name, branches=['scaledmass', 'mjpipi'])
        a_data_x = a_data_x[(a_data_x['mjpipi'] > 3862) & (a_data_x['mjpipi'] < 3882)]['scaledmass']
        d_sig_est = fit_doubleCB(a_mc_x, a_data_x, out_path, s_info='signal_yield_est', )
        print("*** Expected X(3872) signal yield: %d ***"%(d_sig_est['data_sig_yield']))
        print("*** Expected X(3823) signal yield: %d ***"%(float(d_sig_est['data_sig_yield'])/20.))
        s0 = float(d_sig_est['data_sig_yield'])/20.

    d_run_info = {}
    d_roc_plot = {}
    d_run_info['sig_est_params'] = d_sig_est

    for run in list(D_CONFIGS.keys()):
        print('*** Performing run %s ***'%(run))
        d_run_info[run] = {}
        d_roc_plot[run] = {}
        ### Prepare data
        load_branches = list(D_CONFIGS[run]['fit_vars'])
        load_branches.extend(['scaledmass', 'mjpipi'])
        a_mc_j   = root_numpy.root2array(args.mc_j  , treename = args.tree_name, branches = load_branches)
        a_mc_x   = root_numpy.root2array(args.mc_x  , treename = args.tree_name, branches = load_branches)
        a_backgr = root_numpy.root2array(args.backgr, treename = args.tree_name, branches = load_branches)
        a_data   = root_numpy.root2array(args.data  , treename = args.tree_name, branches = load_branches)
        # Convert to dataframes
        df_mc_j   = pd.DataFrame(a_mc_j)
        df_mc_x   = pd.DataFrame(a_mc_x)
        df_backgr = pd.DataFrame(a_backgr)
        df_data   = pd.DataFrame(a_data)
        # Add probabilities if already determined
        if not args.det_prob:
            df_mc_j  ['probability'] = root_numpy.root2array(args.mc_j  , treename = args.tree_name, branches = '_'.join((args.info, run, args.prob_loc)))
            df_mc_x  ['probability'] = root_numpy.root2array(args.mc_x  , treename = args.tree_name, branches = '_'.join((args.info, run, args.prob_loc)))
            df_backgr['probability'] = root_numpy.root2array(args.backgr, treename = args.tree_name, branches = '_'.join((args.info, run, args.prob_loc)))
            df_data  ['probability'] = root_numpy.root2array(args.data  , treename = args.tree_name, branches = '_'.join((args.info, run, args.prob_loc)))
        # Add classification to training set
        df_mc_j['class']   = 1
        df_mc_x['class']   = 2
        df_backgr['class'] = 0
        # Combine into single training set
        df_train = pd.concat([df_mc_j, df_mc_x, df_backgr])
        # Print summary stats
        print('   *** Data loaded ***')
        print('   *** Training events: %d ***'%(df_train.shape[0]))
        print('   ***     Data events: %d ***'%(df_data .shape[0]))
        
        
        ### Cross validate k times
        if args.cross_val:
            print('   *** Cross validating ***')
            # Initiate K-fold object
            k  = 10
            skf = StratifiedKFold(df_train['class'], n_folds=k, shuffle=True, random_state=0)
            # Store cv accuracies
            cv_accuracies = []
            for train_index, test_index in skf:
                df_cv_train = df_train.iloc[train_index]
                df_cv_test  = df_train.iloc[test_index]
                # Train cross validaiton model
                d_cv_model = train_model(D_CONFIGS[run]['max_depth'], D_CONFIGS[run]['learning_rate'], D_CONFIGS[run]['n_estimators'], D_CONFIGS[run]['fit_vars'], df_cv_train, df_cv_test, s0=s0, bkg_cut=args.bck_cut, determine_probas=True)
                cv_accuracies.append(d_cv_model['cv_accuracy'])
            # Determine average accuracy
            cv_accuracy = np.mean(cv_accuracies)
            d_run_info[run]['cv_accuracy'] = cv_accuracy
            print('   *** CV accuracy: %1.2f ***'%(cv_accuracy))
        else:
            cv_accuracy = None
        
        ### Train actual model
        print('   *** Training final model ***')
        d_model = train_model(D_CONFIGS[run]['max_depth'], D_CONFIGS[run]['learning_rate'], D_CONFIGS[run]['n_estimators'], D_CONFIGS[run]['fit_vars'], df_train, s0=s0, bkg_cut=args.bck_cut, determine_probas=args.det_prob)
        # Access some important objects
        if args.det_prob:
            xgbc = d_model['xgb_trained']
        prob_threshold = np.asscalar(d_model['best_cut'])
        d_run_info[run]['optimal_cut'] = prob_threshold
        d_run_info[run]['all_signal_efficiency'] = d_model['sig_effs_all'][d_model['cut_index']]
        d_run_info[run]['mcj_signal_efficiency'] = d_model['sig_effs_mcj'][d_model['cut_index']]
        d_run_info[run]['mcx_signal_efficiency'] = d_model['sig_effs_mcx'][d_model['cut_index']]
        d_roc_plot[run]['sig_effs'] = d_model['sig_effs_all']
        d_roc_plot[run]['bgr_rejs'] = d_model['bgr_rejs']

        # Print some summary stats
        print('   ***           Optimal cut: %1.2f ***'%(prob_threshold))
        print('   ***     Signal efficiency: %1.2f ***'%(d_run_info[run]['all_signal_efficiency']))
        print('   *** MCJ Signal efficiency: %1.2f ***'%(d_run_info[run]['mcj_signal_efficiency']))
        print('   *** MCX Signal efficiency: %1.2f ***'%(d_run_info[run]['mcx_signal_efficiency']))
        
        ### Apply model to data
        if args.det_prob:
            df_train['probability'] = d_model['probability']
            df_data['probability'] = xgbc.predict_proba(df_data[D_CONFIGS[run]['fit_vars']])[:,1]
        df_data['class'] = df_data['probability'] > prob_threshold

        ### Plot cut optimisation
        print('   *** Plotting cut optimisation ***')
        fig = plt.figure()
        plt.plot(d_model['cuts'], d_model['cut_scores'])
        plt.ylabel("Cut Score")
        plt.xlabel("Probability Threshold")
        plt.xlim(0.,1.)
        plt.title("Cut Score "+args.info+" "+run)
        plt.tight_layout(pad=2.0)
        fig.savefig(out_path+'/'+run+'_cut_score.pdf')
        plt.close()

        
        ### Plot mass histogram for optimal cut
        print('   *** Plotting Mass histograms ***')
        # Initialise canvas
        c_name = 'B_Mass_Distribution '+args.info+'_'+run
        c = ROOT.TCanvas(c_name, c_name, 600, 400)
        c.cd()
        # Select required quantity
        a_raw = df_data['scaledmass'].as_matrix()
        a_cut = df_data[df_data['class']==1]['scaledmass'].as_matrix()
        # Create and format histograms
        h_raw = ROOT.TH1F(c_name+'_No_Cut' , c_name+'_No_Cut;B Mass [MeV/#it{c}^{2}];candidates/18[MeV/#it{c}^{2}]', 100, 5220., 5400.)
        h_cut = ROOT.TH1F(c_name+'_XGB_Cut', c_name+'_XGB_Cut;B Mass [MeV/#it{c}^{2}];candidates/18[MeV/#it{c}^{2}]', 100, 5220., 5400.)
        # Fill histograms
        map(h_raw.Fill, a_raw)
        map(h_cut.Fill, a_cut)
        # Normalise
        ## Make it pretty
        h_raw.SetTitle('B Mass Distribution '+args.info+' '+run)
        # Format for each case of x-axis
        h_raw.GetYaxis().SetTitleOffset(1.6)
        y_max = 1.1*max(h_raw.GetBinContent(h_raw.GetMaximumBin()), h_cut.GetBinContent(h_cut.GetMaximumBin()))
        y_min = 0.9*min(h_raw.GetBinContent(h_raw.GetMinimumBin()), h_cut.GetBinContent(h_cut.GetMinimumBin()))
        h_raw.GetYaxis().SetRangeUser(y_min, y_max)
        # Format plotting style
        h_raw.SetLineColor(ROOT.kRed)
        h_raw.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
        h_cut.SetLineColor(ROOT.kBlue)
        h_cut.SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
        # Remove stats boxes
        h_raw.SetStats(False)
        h_cut.SetStats(False)
        # Print
        h_raw.Draw('HIST')
        h_cut.Draw('HISTsame')
        # Create legend
        leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
        leg.AddEntry(h_raw, 'Uncut Data', 'L')
        leg.AddEntry(h_cut, 'XGBoost cut: {:.3f}'.format(prob_threshold), 'L')
        leg.SetLineColor(0)
        leg.SetLineStyle(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw('same')
        # Save
        c.SaveAs(out_path +'/'+run+'_Mass_histogram_B_XGBcut.pdf')


        ### BDT answer
        print('   *** Plotting classification probabilities ***')
        # Initialise canvas
        c_name = 'BDT_Predicted_Probability '+args.info+'_'+run
        c = ROOT.TCanvas(c_name, c_name, 600, 400)
        c.cd()
        # Select required quantity
        a_train_sig_prob = df_train['probability'][df_train['class'] != 0].as_matrix()
        a_train_bkg_prob = df_train['probability'][df_train['class'] == 0].as_matrix()
        a_data_prob      = df_data['probability'].as_matrix()
        # Create and format histograms
        h_train_sig_prob = ROOT.TH1F(c_name+'_Sig_Prob' , c_name+'_Sig_Prob;Probability;Candidates', 100, 0., 1.)
        h_train_bkg_prob = ROOT.TH1F(c_name+'_Bkg_Prob' , c_name+'_Bkg_Prob;Probability;Candidates', 100, 0., 1.)
        h_data_prob      = ROOT.TH1F(c_name+'_Data_Prob', c_name+'_Data_Prob;Probability;Candidates', 100, 0., 1.)
        # Fill histograms
        map(h_train_sig_prob.Fill, a_train_sig_prob)
        map(h_train_bkg_prob.Fill, a_train_bkg_prob)
        map(h_data_prob.Fill, a_data_prob)
        # Normalise
        h_train_sig_prob.Scale(1./h_train_sig_prob.Integral())
        h_train_bkg_prob.Scale(1./h_train_bkg_prob.Integral())
        h_data_prob     .Scale(1./h_data_prob     .Integral())
        ## Make it pretty
        h_train_sig_prob.SetTitle('Event Probability Distribution '+args.info+' '+run)
        # Format for each case of x-axis
        h_train_sig_prob.GetYaxis().SetTitleOffset(1.6)
        y_max = 1.1*max(max(h_train_sig_prob.GetBinContent(h_train_sig_prob.GetMaximumBin()), h_train_bkg_prob.GetBinContent(h_train_bkg_prob.GetMaximumBin())), h_data_prob.GetBinContent(h_data_prob.GetMaximumBin()))
        y_min = 0.9*min(min(h_train_sig_prob.GetBinContent(h_train_sig_prob.GetMinimumBin()), h_train_bkg_prob.GetBinContent(h_train_bkg_prob.GetMinimumBin())), h_data_prob.GetBinContent(h_data_prob.GetMinimumBin()))
        h_train_sig_prob.GetYaxis().SetRangeUser(y_min, y_max)
        # Format plotting style
        h_train_sig_prob.SetLineColor(ROOT.kRed)
        h_train_sig_prob.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
        h_train_bkg_prob.SetLineColor(ROOT.kBlue)
        h_train_bkg_prob.SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
        h_data_prob.SetLineColor(ROOT.kGreen)
        h_data_prob.SetFillColorAlpha(ROOT.kGreen - 10, 0.7)
        # Remove stats boxes
        h_train_sig_prob.SetStats(False)
        h_train_bkg_prob.SetStats(False)
        h_data_prob.SetStats(False)
        # Print
        h_train_sig_prob.Draw('HIST')
        h_train_bkg_prob.Draw('HISTsame')
        h_data_prob     .Draw('HISTsame')
        # Create legend
        leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
        leg.AddEntry(h_train_sig_prob, 'Training Signal Events'    , 'L')
        leg.AddEntry(h_train_bkg_prob, 'Training Background Events', 'L')
        leg.AddEntry(h_data_prob     , 'Data Events'               , 'L')
        leg.SetLineColor(0)
        leg.SetLineStyle(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw('same')
        # Save
        c.SaveAs(out_path+'/'+run+'_BDT_answer.pdf')

        ### Plot variable distributions
        print('   *** Plotting training variable distributions ***')
        for var in D_CONFIGS[run]['fit_vars']:
            # Initialise canvas
            c_name = var+'_Distribution_'+args.info+'_'+run
            c = ROOT.TCanvas(c_name, c_name, 600, 400)
            c.cd()
            # Select required quantity
            a_plt_sig = df_train[var][df_train['probability'] >= prob_threshold].as_matrix()
            a_plt_bkg = df_train[var][df_train['probability'] <  prob_threshold].as_matrix()
            # Create and format histograms
            x_max = max(max(a_plt_sig), max(a_plt_bkg))
            x_min = min(min(a_plt_sig), min(a_plt_bkg))
            h_plt_sig = ROOT.TH1F(c_name+'_Sig', c_name+'_Sig;'+var+';candidates', 100, x_min, x_max)
            h_plt_bkg = ROOT.TH1F(c_name+'_Bkg', c_name+'_Bkg;'+var+';candidates', 100, x_min, x_max)
            # Fill histograms
            map(h_plt_sig.Fill, a_plt_sig)
            map(h_plt_bkg.Fill, a_plt_bkg)
            ## Make it pretty
            h_plt_sig.SetTitle(var+' Distribution '+args.info+' '+run)
            # Format for each case of x-axis
            h_plt_sig.GetYaxis().SetTitleOffset(1.6)
            y_max = 1.1*max(h_plt_sig.GetBinContent(h_plt_sig.GetMaximumBin()), h_plt_bkg.GetBinContent(h_plt_bkg.GetMaximumBin()))
            h_plt_sig.GetYaxis().SetRangeUser(0, y_max)
            h_plt_sig.GetXaxis().SetRangeUser(x_min, x_max)
            # Format plotting style
            h_plt_sig.SetLineColor(ROOT.kRed)
            h_plt_sig.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
            h_plt_bkg.SetLineColor(ROOT.kBlue)
            h_plt_bkg.SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
            # Remove stats boxes
            h_plt_sig.SetStats(False)
            h_plt_bkg.SetStats(False)
            # Print
            h_plt_sig.Draw('HIST')
            h_plt_bkg.Draw('HISTsame')
            # Create legend
            leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
            leg.AddEntry(h_plt_sig, 'Training Events Identified as Signal', 'L')
            leg.AddEntry(h_plt_bkg, 'Training Events Identified as Background', 'L')
            leg.SetLineColor(0)
            leg.SetLineStyle(0)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.Draw('same')
            # Save
            c.SaveAs(out_path+'/'+run+'_'+var+'_distribution.pdf')

        ### Plot comparison to MC data
        print('   *** Plotting comparison to psi(2S) MC ***')
        df_comp = df_data[((df_data['mjpipi'] < 3696) & (df_data['mjpipi'] > 3676)) & ((df_data['scaledmass'] < 5299) & (df_data['scaledmass'] > 5259))]
        df_side = df_data[((df_data['mjpipi'] < 3696) & (df_data['mjpipi'] > 3676)) & ((df_data['scaledmass'] < 5400) & (df_data['scaledmass'] > 5350))]
        for var in D_CONFIGS[run]['fit_vars']:
            # Initialise canvas
            c_name = var+'_MC_#psi(2S)_Comparison_'+args.info+'_'+run
            c = ROOT.TCanvas(c_name, c_name, 600, 400)
            c.cd()
            # Select required quantity
            a_mc_j = df_mc_j[var].as_matrix()
            a_comp = df_comp[var].as_matrix()
            a_side = df_side[var].as_matrix()
            # Create and format histograms
            x_max = max(max(a_mc_j), max(a_comp))
            x_min = min(min(a_mc_j), min(a_comp))
            h_mc_j = ROOT.TH1F(c_name+'_mc_#psi', c_name+'_mc_#psi;'+var+';candidates', 100, x_min, x_max)
            h_comp = ROOT.TH1F(c_name+'_data'   , c_name+'_data;'+var+';candidates'   , 100, x_min, x_max)
            h_side = ROOT.TH1F(c_name+'_side'   , c_name+'_side;'+var+';candidates'   , 100, x_min, x_max)
            # Fill histograms
            map(h_mc_j.Fill, a_mc_j)
            map(h_comp.Fill, a_comp)
            map(h_side.Fill, a_side)
            # Background reduce
            h_comp.Add(h_side, -1)
            # Normalise
            h_mc_j.Scale(1./h_mc_j.Integral())
            h_comp.Scale(1./h_comp.Integral())
            ## Make it pretty
            h_mc_j.SetTitle(var+' Data vs MC J(2S) Distribution '+args.info+' '+run)
            # Format for each case of x-axis
            h_mc_j.GetYaxis().SetTitleOffset(1.6)
            y_max = 1.1*max(h_mc_j.GetBinContent(h_mc_j.GetMaximumBin()), h_comp.GetBinContent(h_comp.GetMaximumBin()))
            y_min = 0.9*min(h_mc_j.GetBinContent(h_mc_j.GetMinimumBin()), h_comp.GetBinContent(h_comp.GetMinimumBin()))
            h_mc_j.GetYaxis().SetRangeUser(y_min, y_max)
            # Format plotting style
            h_mc_j.SetLineColor(ROOT.kRed)
            h_mc_j.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
            h_comp.SetLineColor(ROOT.kBlue)
            h_comp.SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
            # Remove stats boxes
            h_mc_j.SetStats(False)
            h_comp.SetStats(False)
            # Print
            h_mc_j.Draw('HIST')
            h_comp.Draw('HISTsame')
            # Create legend
            leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
            leg.AddEntry(h_mc_j, '#psi(2S) Monte-Carlo', 'L')
            leg.AddEntry(h_comp, 'Background Reduced Data in #psi(2S) Region', 'L')
            leg.SetLineColor(0)
            leg.SetLineStyle(0)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.Draw('same')
            # Save
            c.SaveAs(out_path+'/'+run+'_'+var+'_XGBvsMC_distribution.pdf')


        ## Perform fit to XGB cut data
        # Filter dataframes
        a_cut_mc   = df_train[df_train['class'] != 0]['scaledmass'].as_matrix()
        a_cut_data = df_data[df_data['class'] == 1]['scaledmass'].as_matrix()
        # Fit
        d_cut_fit = fit_doubleCB(a_cut_mc, a_cut_data, out_path, s_info=run+'_cut_data_plot')
        # Store params
        d_run_info[run]['cut_fit_params'] = d_cut_fit
        # Fit in X region only
        a_mc_x = root_numpy.root2array(args.mc_x, treename = args.tree_name, branches='scaledmass')
        a_data_x = df_data[(df_data['class'] == 1) & ((df_data['mjpipi'] > 3862) & (df_data['mjpipi'] < 3882))]['scaledmass'].as_matrix()
        d_sig_est = fit_doubleCB(a_mc_x, a_data_x, out_path, s_info=run+'x_signal_yield_est')
        d_run_info[run]['x_reg_fit_params'] = d_sig_est


        ## Update root files with probabilitities
        if args.det_prob:
            # Extract probabilities
            a_mc_j_probs = np.array(df_train[df_train['class'] == 1]['probability'].as_matrix(), dtype = [('_'.join((args.info, run, args.prob_loc)), 'f8')])
            a_mc_x_probs = np.array(df_train[df_train['class'] == 2]['probability'].as_matrix(), dtype = [('_'.join((args.info, run, args.prob_loc)), 'f8')])
            a_bgrd_probs = np.array(df_train[df_train['class'] == 0]['probability'].as_matrix(), dtype = [('_'.join((args.info, run, args.prob_loc)), 'f8')])
            a_data_probs = np.array(                         df_data['probability'].as_matrix(), dtype = [('_'.join((args.info, run, args.prob_loc)), 'f8')])
            # Update root files
            root_numpy.array2root(a_mc_j_probs, args.mc_j  , args.tree_name)
            root_numpy.array2root(a_mc_x_probs, args.mc_x  , args.tree_name)
            root_numpy.array2root(a_bgrd_probs, args.backgr, args.tree_name)
            root_numpy.array2root(a_data_probs, args.data  , args.tree_name)

        

    print('*** Plotting ROC curve ***')
    ### Plot ROC curve
    fig = plt.figure()
    for run in list(D_CONFIGS.keys()):
        plt.plot(d_roc_plot[run]['bgr_rejs'], d_run_info[run]['sig_effs'], label=run)
    plt.legend(loc=3)
    plt.ylabel("Background Rejection")
    plt.xlabel("Signal Efficiency")
    plt.xlim(0.,1.)
    plt.ylim(0.,1.)
    plt.title("ROC Curve "+args.info)
    plt.tight_layout(pad=2.0)
    fig.savefig(out_path+'/ROC_curve.pdf')
    plt.close()

    print('*** Dumping run information ***')
    with open(out_path + '/' + args.out_dict, 'w') as outfile:
        yaml.dump(d_run_info, outfile, default_flow_style=False)


## Command line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "--add description--")
    parser.add_argument("mc_x"       ,            default = None                                   	   , help = "input MC X(3823) files")
    parser.add_argument("mc_j"       ,            default = None                                          , help = "input MC J(2S) files")
    parser.add_argument("backgr"     ,            default = None                                          , help = "input root background files")
    parser.add_argument("data"       ,            default = None                                          , help = "input data files")
    parser.add_argument("--out_path" , "-o"     , default = '/Disk/ds-sopa-group/PPE/users/s1305440/project_stuff/data/results'      , help = "directory for storing results")
    parser.add_argument("--config"   , "-c"     , default = '/Disk/ds-sopa-group/PPE/users/s1305440/project_stuff/config/xgboost.yml', help =  "config file")
    parser.add_argument("--out_dict" , "-q"     , default = 'run_info.yml'                                             , help = "dictionary summarising run information")
    parser.add_argument("--tree_name", "-t"     , default = "DecayTree"                                                , help = "input tree name")
    parser.add_argument("--info"     , "-i"     , default = ""                                                         , help = "string specifying information about the run")
    parser.add_argument("--find_s0"  , "-s"     , action = 'store_true'                                                , help = "if specified, a fit will be made to estimate s0 for X(3823)")
    parser.add_argument("--bck_cut"  , "-b"     , action = 'store_true'                                                , help = "if specified, the optimal cut will be where background rejection reaches 99%")
    parser.add_argument("--cross_val", "-v"     , action = 'store_true')
    parser.add_argument("--prob_loc" , "-p"     , default = ""                                                         , help = "if not default, specify the branch name where probabilities are stored")
    parser.add_argument("--det_prob" , "-d"     , action = 'store_true')
    args = parser.parse_args()
    main(args)
