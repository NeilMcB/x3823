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

def fit_doubleCB(a_mc_x, a_data, out_path, s_info=''):
    # Initialise dictionary for storing fit info
    d_fit_info = {}
    # Format arrays
    a_fit_mc   = a_mc_x.astype(dtype=[('b_mass_mc', np.float)])
    a_fit_data = a_data.astype(dtype=[('b_mass'   , np.float)])
    # Estimate sig/bkg yields
    max_yield_est = a_fit_data.shape[0]                         # max possible signal yield
    bkg_yield_est = len(a_fit_data[a_fit_data['b_mass'] > 5350.]) * 4.     # estimate background from sideband
    # Create tree for fitting
    t_fit_data = root_numpy.array2tree(a_fit_data)
    t_fit_mc   = root_numpy.array2tree(a_fit_mc)
    
    ## Monte Carlo
    # Inialise parameters for fit
    b_mass_mc     = ROOT.RooRealVar("b_mass_mc"   , "B mass MC [MeV]", 5200. , 5400.)
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
    model_sig_mc.fitTo(dataset_mc, ROOT.RooFit.Range(5220., 5400.))
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
    c.SaveAs(out_path+s_info+'_FittedMassDistribution_MonteCarlo.pdf')
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
    b_mass    = ROOT.RooRealVar("b_mass"   , "B mass [MeV]" , 5220., 5400.)
    mean      = ROOT.RooRealVar("mean"     , "mean"         , mean_mc     .getValV(), 5220., 5400.)
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
    c.SaveAs(out_path+s_info+'_FittedMassDistribution_Data.pdf')
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


def main(args):
	### Perpare data for processing

    # Ensure output directory exists
    out_path = args.data_dir+'results/'
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    ## Load files
    l_fit_vars = ['logDIRA', 'log_bplus_IPCHI2_OWNPV', 'bplus_LOKI_DTF_CHI2NDOF', 'log_bplus_FDCHI2_OWNPV', 'bplus_ETA', 
                  'log_1_IPCHI2_OWNPV', 'log_2_IPCHI2_OWNPV', 'log_3_IPCHI2_OWNPV', 'log_4_IPCHI2_OWNPV', 'log_5_IPCHI2_OWNPV', 
                  'mu_PT_max', 'mu_PT_min']
    l_mass_vars = ['scaledmass', 'mjpipi']
    l_load_branches = l_fit_vars + l_mass_vars + ['prob']
    # Load files into arrays
    print('*** Loading Data ***')
    a_mc_x = root_numpy.root2array(args.data_dir+'mc_x_proba.root', treename = args.tree_name, branches = l_load_branches)
    a_mc_p = root_numpy.root2array(args.data_dir+'mc_p_proba.root', treename = args.tree_name, branches = l_load_branches)
    a_side = root_numpy.root2array(args.data_dir+'side_proba.root', treename = args.tree_name, branches = l_load_branches)
    a_data = root_numpy.root2array(args.data_dir+'data_proba.root', treename = args.tree_name, branches = l_load_branches)

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
    # Combine into training set
    df_train = pd.concat([df_mc_x, df_mc_p, df_side])
    # Print summary stats
    print('   *** Data loaded ***')
    print('   *** Training events: %d ***'%(df_train.shape[0]))
    print('   ***     Data events: %d ***'%(df_data .shape[0]))

    # Dictionaries for storing information on each run
    d_info = {}

    ### Estimate signal yield - all data
    d_sig_est_alldata = fit_doubleCB(pd.concat([df_mc_x, df_mc_p])['scaledmass'].as_matrix(), df_data['scaledmass'].as_matrix(), out_path, s_info='alldata_signal_est')
    ### Estimate signal yield - X region
    if args.s0 is not None:
        s0 = int(args.s0)
    else:
        df_data_p = df_data[(df_data['mjpipi'] > 3676) & (df_data['mjpipi'] < 3696)]
        df_data_x = df_data[(df_data['mjpipi'] > 3862) & (df_data['mjpipi'] < 3882)]
        d_sig_est_p = fit_doubleCB(df_mc_p['scaledmass'].as_matrix(), df_data_p['scaledmass'].as_matrix(), out_path, s_info='psi(2S)_s0_est', )
        d_sig_est_x = fit_doubleCB(df_mc_x['scaledmass'].as_matrix(), df_data_x['scaledmass'].as_matrix(), out_path, s_info='x(3872)_s0_est', )
        print("*** Expected psi(2S) signal yield: %d ***"%(d_sig_est_p['data_sig_yield']))
        print("*** Expected X(3872) signal yield: %d ***"%(d_sig_est_x['data_sig_yield']))
        print("*** Expected X(3823) signal yield: %d ***"%(float(d_sig_est_x['data_sig_yield'])/20.))
        s0 = float(d_sig_est_x['data_sig_yield'])/20.
        d_info['sig_est_x_reg'] = d_sig_est_x
    
    d_info = {}

    d_info['raw_fit_params'] = d_sig_est_alldata

    out_path_plots = out_path+'plots/'
    if not os.path.exists(out_path_plots):
        os.makedirs(out_path_plots)

    if args.opt_cut is None:
        ### Find optimal cut
        print('   *** Determining optimal cut ***')
        # Optimise the probability cut
        sig_effs_mcp = []  # record signal efficiencies - on MC psi(2S) only
        sig_effs_mcx = []  # record signal efficiencies - on MC X(3823) only
        sig_effs_all = []  # record signal efficiencies
        bgr_rejs     = []  # record background rejections
        cut_scores   = []  # record cut optimisation metric
        # Determine cut metric for a range of cuts
        cuts = np.linspace(.0, 1., 200, endpoint=False)
        for prob_threshold in cuts:
            # Determine how many predictions are correct
            signal_efficiency_mcx = float(df_train[(df_train['prob'] > prob_threshold) & (df_train['cat']   == 'mc_x')].shape[0]) / float(df_train[df_train['cat']   == 'mc_x'].shape[0])
            signal_efficiency_mcp = float(df_train[(df_train['prob'] > prob_threshold) & (df_train['cat']   == 'mc_p')].shape[0]) / float(df_train[df_train['cat']   == 'mc_p'].shape[0])
            signal_efficiency_all = float(df_train[(df_train['prob'] > prob_threshold) & (df_train['class'] == 1     )].shape[0]) / float(df_train[df_train['class'] == 1     ].shape[0])
            background_rejection  = float(df_train[(df_train['prob'] < prob_threshold) & (df_train['class'] == 0     )].shape[0]) / float(df_train[df_train['class'] == 0     ].shape[0])
            # Store scores
            sig_effs_all.append(signal_efficiency_all)
            sig_effs_mcp.append(signal_efficiency_mcp)
            sig_effs_mcx.append(signal_efficiency_mcx)
            bgr_rejs    .append(background_rejection)
            # Optimize cut
            eff = signal_efficiency_all
            a   = 5. # expected significance
            # Background events, scaled to 40MeV window about B peak, considering only those in X(3823) region
            B = df_train[((df_train['prob'] > prob_threshold)) & ((df_train['scaledmass'] > 5400.) & (df_train['scaledmass'] < 5450.)) & ((df_train['mjpipi'] > 3773) & (df_train['mjpipi'] < 3873))].shape[0] * .8
            cut_scores.append((s0 * eff) / sqrt((s0 * eff) + B))
        # Find optimal cut
        cut_index = np.argmax(cut_scores)
        prob_threshold  = cuts[cut_index]

        ### Store some parameters of interest
        d_info['optimal_cut'] = np.asscalar(prob_threshold)
        d_info['all_signal_efficiency'] = sig_effs_all[cut_index]
        d_info['mcp_signal_efficiency'] = sig_effs_mcp[cut_index]
        d_info['mcx_signal_efficiency'] = sig_effs_mcx[cut_index]
        ### Print some summary stats
        print('   ***           Optimal cut: %1.2f ***'%(prob_threshold))
        print('   ***     Signal efficiency: %1.2f ***'%(d_info['all_signal_efficiency']))
        print('   *** MCP Signal efficiency: %1.2f ***'%(d_info['mcp_signal_efficiency']))
        print('   *** MCX Signal efficiency: %1.2f ***'%(d_info['mcx_signal_efficiency']))
    else:
        prob_threshold = float(args.opt_cut)

    ### Apply model to data
    df_data['class'] = df_data['prob'] > prob_threshold

    ### Plot cut optimisation
    print('   *** Plotting cut optimisation ***')
    fig = plt.figure()
    plt.plot(cuts, cut_scores)
    plt.ylabel("Cut Score")
    plt.xlabel("Probability Threshold")
    plt.xlim(0.,1.)
    plt.title("Cut Score")
    plt.tight_layout(pad=2.0)
    fig.savefig(out_path_plots+'cut_score.pdf')
    plt.close()

    ### Plot mass histogram for optimal cut
    print('   *** Plotting Mass histograms ***')
    # Initialise canvas
    c_name = 'B_Mass_Distribution'
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
    h_raw.SetTitle('B Mass Distribution')
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
    h_raw.Draw('E0HIST')
    h_cut.Draw('E0HISTsame')
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
    c.SaveAs(out_path_plots+'Mass_histogram_B_XGBcut.pdf')

    ### BDT answer
    print('   *** Plotting classification probabilities ***')
    # Initialise canvas
    c_name = 'BDT_Predicted_Probability'
    c = ROOT.TCanvas(c_name, c_name, 600, 400)
    c.cd()
    # Select required quantity
    a_train_sig_prob = df_train['prob'][df_train['class'] == 1].as_matrix()
    a_train_bkg_prob = df_train['prob'][df_train['class'] == 0].as_matrix()
    a_data_prob      = df_data['prob'].as_matrix()
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
    h_train_sig_prob.SetTitle('Event Probability Distribution')
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
    h_train_sig_prob.Draw('E0HIST')
    h_train_bkg_prob.Draw('E0HISTsame')
    h_data_prob     .Draw('E0HISTsame')
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
    c.SaveAs(out_path_plots+'BDT_answer.pdf')

    ### Plot variable distributions
    print('   *** Plotting training variable distributions ***')
    out_path_var = out_path_plots+'dist_vars/'
    if not os.path.exists(out_path_var):
        os.makedirs(out_path_var)
    for var in l_fit_vars + ['prob']:
        # Initialise canvas
        c_name = var+'_Distribution'
        c = ROOT.TCanvas(c_name, c_name, 600, 400)
        c.cd()
        # Select required quantity
        a_plt_sig = df_train[var][df_train['prob'] >= prob_threshold].as_matrix()
        a_plt_bkg = df_train[var][df_train['prob'] <  prob_threshold].as_matrix()
        
        # Create and format histograms
        x_max = max(max(a_plt_sig), max(a_plt_bkg))
        x_min = min(min(a_plt_sig), min(a_plt_bkg))
        h_plt_sig = ROOT.TH1F(c_name+'_Sig', c_name+'_Sig;'+var+';candidates', 100, x_min, x_max)
        h_plt_bkg = ROOT.TH1F(c_name+'_Bkg', c_name+'_Bkg;'+var+';candidates', 100, x_min, x_max)
        # Fill histograms
        map(h_plt_sig.Fill, a_plt_sig)
        map(h_plt_bkg.Fill, a_plt_bkg)
        ## Make it pretty
        h_plt_sig.SetTitle(var+' Distribution')
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
        h_plt_sig.Draw('E0HIST')
        h_plt_bkg.Draw('E0HISTsame')
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
        c.SaveAs(out_path_var+var+'.pdf')

        ### Plot comparison to MC data
        out_path_mcp_data = out_path_plots+'mcp_v_data/'
        if not os.path.exists(out_path_mcp_data):
            os.makedirs(out_path_mcp_data)
        print('   *** Plotting comparison to psi(2S) MC ***')
        df_data_comp = df_data[((df_data['mjpipi'] < 3696) & (df_data['mjpipi'] > 3676)) & ((df_data['scaledmass'] < 5299) & (df_data['scaledmass'] > 5259))]
        df_side_comp = df_side[ (df_side['mjpipi'] < 3696) & (df_side['mjpipi'] > 3676)]
        for var in l_fit_vars + ['prob']:
            # Initialise canvas
            c_name = var+'_MC_#psi(2S)_Comparison'
            c = ROOT.TCanvas(c_name, c_name, 600, 400)
            c.cd()
            # Select required quantity
            a_mc_p = df_mc_p[var].as_matrix()
            a_data_comp = df_data_comp[var].as_matrix()
            a_side_comp = df_side_comp[var].as_matrix()
            
            # Create and format histograms
            x_max = max(max(a_mc_p), max(a_data_comp))
            x_min = min(min(a_mc_p), min(a_data_comp))
            h_mc_p = ROOT.TH1F(c_name+'_mc_#psi', c_name+'_mc_#psi;'+var+';candidates', 100, x_min, x_max)
            h_comp = ROOT.TH1F(c_name+'_data'   , c_name+'_data;'+var+';candidates'   , 100, x_min, x_max)
            h_side = ROOT.TH1F(c_name+'_side'   , c_name+'_side;'+var+';candidates'   , 100, x_min, x_max)
            # Fill histograms
            map(h_mc_p.Fill, a_mc_p)
            map(h_comp.Fill, a_data_comp)
            map(h_side.Fill, a_side_comp)
            # Background reduce
            h_comp.Add(h_side, -1)
            # Normalise
            h_mc_p.Scale(1./h_mc_p.Integral())
            h_comp.Scale(1./h_comp.Integral())
            h_side.Scale(1./h_side.Integral())
            ## Make it pretty
            h_mc_p.SetTitle(var+' Data vs MC J(2S) Distribution')
            # Format for each case of x-axis
            h_mc_p.GetYaxis().SetTitleOffset(1.6)
            y_max = 1.1*max((h_mc_p.GetBinContent(h_mc_p.GetMaximumBin()), h_comp.GetBinContent(h_comp.GetMaximumBin()), h_side.GetBinContent(h_side.GetMaximumBin())))
            y_min = 0.9*min((h_mc_p.GetBinContent(h_mc_p.GetMinimumBin()), h_comp.GetBinContent(h_comp.GetMinimumBin()), h_side.GetBinContent(h_side.GetMinimumBin())))
            h_mc_p.GetYaxis().SetRangeUser(y_min, y_max)
            # Format plotting style
            h_mc_p.SetLineColor(ROOT.kRed)
            h_mc_p.SetFillColorAlpha(ROOT  .kRed - 10, 0.7)
            h_comp.SetLineColor(ROOT.kBlue)
            h_comp.SetFillColorAlpha(ROOT .kBlue - 10, 0.7)
            h_side.SetLineColor(ROOT.kGreen)
            h_side.SetFillColorAlpha(ROOT.kGreen - 10, 0.7)
            # Remove stats boxes
            h_mc_p.SetStats(False)
            h_comp.SetStats(False)
            h_side.SetStats(False)
            # Print
            h_mc_p.Draw('E0HIST')
            h_comp.Draw('E0HISTsame')
            h_side.Draw('E0HISTsame')
            # Create legend
            leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
            leg.AddEntry(h_mc_p, '#psi(2S) Monte-Carlo', 'L')
            leg.AddEntry(h_comp, 'Background Reduced Data in #psi(2S) Region', 'L')
            leg.AddEntry(h_side, 'Background Data in #psi(2S) Region', 'L')
            leg.SetLineColor(0)
            leg.SetLineStyle(0)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.Draw('same')
            # Save
            c.SaveAs(out_path_mcp_data+var+'.pdf')

        ### Plot comparison to MC data
        out_path_mcx_data = out_path_plots+'mcx_v_data/'
        if not os.path.exists(out_path_mcx_data):
            os.makedirs(out_path_mcx_data)
        print('   *** Plotting comparison to X(3872) MC ***')
        df_data_comp = df_data[((df_data['mjpipi'] < 3882) & (df_data['mjpipi'] > 3862)) & ((df_data['scaledmass'] < 5299) & (df_data['scaledmass'] > 5259))]
        df_side_comp = df_side[ (df_side['mjpipi'] < 3882) & (df_side['mjpipi'] > 3862)]
        for var in l_fit_vars + ['prob']:
            # Initialise canvas
            c_name = var+'_MC_X(3872)_Comparison'
            c = ROOT.TCanvas(c_name, c_name, 600, 400)
            c.cd()
            # Select required quantity
            a_mc_p = df_mc_x[var].as_matrix()
            a_data_comp = df_data_comp[var].as_matrix()
            a_side_comp = df_side_comp[var].as_matrix()
            # Scale DIRA and IPCHI2
            
            # Create and format histograms
            x_max = max(max(a_mc_p), max(a_data_comp))
            x_min = min(min(a_mc_p), min(a_data_comp))
            h_mc_p = ROOT.TH1F(c_name+'_mc_#psi', c_name+'_mc_#psi;'+var+';candidates', 100, x_min, x_max)
            h_comp = ROOT.TH1F(c_name+'_data'   , c_name+'_data;'+var+';candidates'   , 100, x_min, x_max)
            h_side = ROOT.TH1F(c_name+'_side'   , c_name+'_side;'+var+';candidates'   , 100, x_min, x_max)
            # Fill histograms
            map(h_mc_p.Fill, a_mc_p)
            map(h_comp.Fill, a_data_comp)
            map(h_side.Fill, a_side_comp)
            # Background reduce
            h_comp.Add(h_side, -1)
            # Normalise
            h_mc_p.Scale(1./h_mc_p.Integral())
            h_comp.Scale(1./h_comp.Integral())
            h_side.Scale(1./h_side.Integral())
            ## Make it pretty
            h_mc_p.SetTitle(var+' Data vs MC X(3872) Distribution')
            # Format for each case of x-axis
            h_mc_p.GetYaxis().SetTitleOffset(1.6)
            y_max = 1.1*max((h_mc_p.GetBinContent(h_mc_p.GetMaximumBin()), h_comp.GetBinContent(h_comp.GetMaximumBin()), h_side.GetBinContent(h_side.GetMaximumBin())))
            y_min = 0.9*min((h_mc_p.GetBinContent(h_mc_p.GetMinimumBin()), h_comp.GetBinContent(h_comp.GetMinimumBin()), h_side.GetBinContent(h_side.GetMinimumBin())))
            h_mc_p.GetYaxis().SetRangeUser(y_min, y_max)
            # Format plotting style
            h_mc_p.SetLineColor(ROOT.kRed)
            h_mc_p.SetFillColorAlpha(ROOT  .kRed - 10, 0.7)
            h_comp.SetLineColor(ROOT.kBlue)
            h_comp.SetFillColorAlpha(ROOT .kBlue - 10, 0.7)
            h_side.SetLineColor(ROOT.kGreen)
            h_side.SetFillColorAlpha(ROOT.kGreen - 10, 0.7)
            # Remove stats boxes
            h_mc_p.SetStats(False)
            h_comp.SetStats(False)
            h_side.SetStats(False)
            # Print
            h_mc_p.Draw('E0HIST')
            h_comp.Draw('E0HISTsame')
            h_side.Draw('E0HISTsame')
            # Create legend
            leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
            leg.AddEntry(h_mc_p, 'X(3872) Monte-Carlo', 'L')
            leg.AddEntry(h_comp, 'Background Reduced Data in X(3872) Region', 'L')
            leg.AddEntry(h_side, 'Background Data in X(3872) Region', 'L')
            leg.SetLineColor(0)
            leg.SetLineStyle(0)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.Draw('same')
            # Save
            c.SaveAs(out_path_mcx_data+var+'.pdf')

        ### Plot comparison of MC data
        out_path_mc_mc = out_path_plots+'mcp_v_mcx/'
        if not os.path.exists(out_path_mc_mc):
            os.makedirs(out_path_mc_mc)
        print('   *** Plotting comparison of psi(2S) MC and X(3972) MC ***')
        for var in l_fit_vars + ['prob']:
            # Initialise canvas
            c_name = var+'_MC_Comparison'
            c = ROOT.TCanvas(c_name, c_name, 600, 400)
            c.cd()
            # Select required quantity
            a_mc_p = df_mc_p[var].as_matrix()
            a_mc_x = df_mc_x[var].as_matrix()
            # Scale DIRA and IPCHI2
            
            # Create and format histograms
            x_max = max(max(a_mc_p), max(a_mc_x))
            x_min = min(min(a_mc_p), min(a_mc_x))
            h_mc_p = ROOT.TH1F(c_name+'_mc_#psi', c_name+'_mc_#psi;'+var+';candidates', 100, x_min, x_max)
            h_mc_x = ROOT.TH1F(c_name+'_mc_X'   , c_name+'_mc_X;'+var+';candidates'   , 100, x_min, x_max)
            # Fill histograms
            map(h_mc_p.Fill, a_mc_p)
            map(h_mc_x.Fill, a_mc_x)
            # Normalise
            h_mc_p.Scale(1./h_mc_p.Integral())
            h_mc_x.Scale(1./h_mc_x.Integral())
            ## Make it pretty
            h_mc_p.SetTitle(var+' MC X(3872) vs MC #psi(2S) Distribution')
            # Format for each case of x-axis
            h_mc_p.GetYaxis().SetTitleOffset(1.6)
            y_max = 1.1*max(h_mc_p.GetBinContent(h_mc_p.GetMaximumBin()), h_mc_x.GetBinContent(h_mc_x.GetMaximumBin()))
            y_min = 0.9*min(h_mc_p.GetBinContent(h_mc_p.GetMinimumBin()), h_mc_x.GetBinContent(h_mc_x.GetMinimumBin()))
            h_mc_p.GetYaxis().SetRangeUser(y_min, y_max)
            # Format plotting style
            h_mc_p.SetLineColor(ROOT.kRed)
            h_mc_p.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
            h_mc_x.SetLineColor(ROOT.kBlue)
            h_mc_x.SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
            # Remove stats boxes
            h_mc_p.SetStats(False)
            h_mc_x.SetStats(False)
            # Print
            h_mc_p.Draw('E0HIST')
            h_mc_x.Draw('E0HISTsame')
            # Create legend
            leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
            leg.AddEntry(h_mc_p, '#psi(2S) Monte-Carlo', 'L')
            leg.AddEntry(h_mc_x, 'X(3872) Monte-Carlo', 'L')
            leg.SetLineColor(0)
            leg.SetLineStyle(0)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.Draw('same')
            # Save
            c.SaveAs(out_path_mc_mc+var+'.pdf')

        ## Perform fit to XGB cut data
        # Filter dataframes
        a_cut_mc   = df_train[df_train['class'] == 1]['scaledmass'].as_matrix()
        a_cut_data = df_data[df_data['class'] == 1]['scaledmass'].as_matrix()
        # Fit
        d_cut_fit = fit_doubleCB(a_cut_mc, a_cut_data, out_path_plots, s_info='cut_data_plot')
        # Store params
        d_info['cut_fit_params'] = d_cut_fit
        # Fit in X region only
        #a_mc_x = df_train[df_train['cat']=='mc_x']['scaledmass'].as_matrix()
        #a_data_x = df_data[(df_data['class'] == 1) & ((df_data['mjpipi'] > 3862) & (df_data['mjpipi'] < 3882))]['scaledmass'].as_matrix()
        #d_sig_est = fit_doubleCB(a_mc_x, a_data_x, out_path_plots, s_info='x_signal_yield_est')
        #d_info['x_reg_fit_params'] = d_sig_est
        # Estimate signal efficiency in the data
        est_data_eff = float(d_cut_fit['data_sig_yield'])/d_sig_est_alldata['data_sig_yield']
        d_info['estimated_data_efficiency'] = est_data_eff
        print(   '*** Estimated fitted signal efficiency: {:.3f} ***'.format(est_data_eff))

    print('*** Dumping run information ***')
    with open(out_path + args.out_dict, 'w') as outfile:
        yaml.dump(d_info, outfile, default_flow_style=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "--add description--")
    parser.add_argument("data_dir"   ,       default = None           , help = "directory containing trained root files")
    parser.add_argument("--out_path" , "-o", default = 'results/'     , help = "directory for storing results")
    parser.add_argument("--out_dict" , "-q", default = 'run_info.yml' , help = "dictionary summarising run information")
    parser.add_argument("--tree_name", "-t", default = "DecayTree"    , help = "input tree name")
    parser.add_argument("--s0"       , "-s", default = None           , help = "if not specified, a fit will be made to estimate s0 for X(3823)")
    parser.add_argument("--opt_cut"  , "-r", default = None           , help = "if specified, skip optimisation and take this as the cut factor")
    args = parser.parse_args()
    main(args)