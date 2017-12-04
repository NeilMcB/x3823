# -*- coding: utf-8 -*-
import root_numpy
import ROOT
import numpy as np
import pandas as pd
from shutil import copyfile
ROOT.gROOT.SetBatch(True)

# Load data B mass branches
print('*** Loading Data ***')
data_loc = '/home/s1305440/PPE_disk/project_stuff/data/data_Qcut.root'
a_data = root_numpy.root2array(data_loc, treename = 'DecayTree', branches = ['scaledmass', 'mppp', 'mjprp', 'mjpk'])
# Load RapidSim B mass branches
Bu2Jpsipipipi_loc = '/home/s1305440/PPE_disk/project_stuff/RapidSim/validation/Bu2Jpsipipipi_tree.root'
a_Bu2Jpsipipipi   = root_numpy.root2array(Bu2Jpsipipipi_loc, treename = 'DecayTree', branches = ['Bp_0_M_pip_12Kp', 'Bp_0_M'])
Bu2JpsipipiK_loc  = '/home/s1305440/PPE_disk/project_stuff/RapidSim/validation/Bu2JpsipipiK_tree.root'
a_Bu2JpsipipiK    = root_numpy.root2array(Bu2JpsipipiK_loc , treename = 'DecayTree', branches = ['Bp_0_M_Kp_02pip', 'Bp_0_M'])
B02psi2skpi_loc   = '/home/s1305440/PPE_disk/project_stuff/RapidSim/validation/B02psi2skpi_tree.root'
a_B02psi2skpi     = root_numpy.root2array(B02psi2skpi_loc  , treename = 'DecayTree', branches = ['m_pim_1_drop', 'B0_0_M'])
Bs2psi2Sphi_loc   = '/home/s1305440/PPE_disk/project_stuff/RapidSim/validation/Bs2psi2Sphi_tree.root'
a_Bs2psi2Sphi     = root_numpy.root2array(Bs2psi2Sphi_loc  , treename = 'DecayTree', branches = ['m_Km_0_drop', 'Bsst0_0_M'])

print('   *** Plotting Bs -> psi(2S) phi ***')
# Initialise canvas
c_name = 'Bs2psi2Sphi'
c = ROOT.TCanvas(c_name, c_name, 600, 400)
c.cd()
# Select required quantity
a_true = a_data['scaledmass']
a_sim  = a_Bs2psi2Sphi['m_Km_0_drop'] * 1000. # Convert to MeV
# Create and format histograms
h_data = ROOT.TH1F(c_name+'_data', c_name+'_data;B Mass [MeV/#it{c}^{2}];candidates', 200, 4200., 5500.)
h_sim  = ROOT.TH1F(c_name+'_sim' , c_name+'_sim ;B Mass [MeV/#it{c}^{2}];candidates', 200, 4200., 5500.)
# Fill histograms
map(h_data.Fill, a_true)
map(h_sim .Fill, a_sim)
# Normalise
h_data.Scale(1./h_data.Integral())
h_sim .Scale(1./h_sim .Integral())
## Make it pretty
h_data.SetTitle('Rapid Sim B Mass Distribution')
# Format for each case of x-axis
h_data.GetYaxis().SetTitleOffset(1.6)
y_max = 1.1*max(h_data.GetBinContent(h_data.GetMaximumBin()), h_sim.GetBinContent(h_sim.GetMaximumBin()))
y_min = 0.9*min(h_data.GetBinContent(h_data.GetMinimumBin()), h_sim.GetBinContent(h_sim.GetMinimumBin()))
h_data.GetYaxis().SetRangeUser(y_min, y_max)
# Format plotting style
h_data.SetLineColor(ROOT.kRed)
h_data.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
h_sim .SetLineColor(ROOT.kBlue)
h_sim .SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
# Remove stats boxes
h_data.SetStats(False)
h_sim .SetStats(False)
# Print
h_data.Draw('HIST')
h_sim .Draw('HISTsame')
# Create legend
leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
leg.AddEntry(h_data, 'B^{+} #rightarrow j/#psi #pi^{+}#pi^{-}K^{+}', 'L')
leg.AddEntry(h_sim , 'B^{*}_{S} #rightarrow #psi(2S)#phi', 'L')
leg.SetLineColor(0)
leg.SetLineStyle(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw('same')
# Save
c.SaveAs('/home/s1305440/PPE_disk/project_stuff/data/rapdiSimComp_Bs2psi2Sphi.pdf')




print('   *** Plotting B0 -> psi(2S) K pi ***')
# Initialise canvas
c_name = 'B02psi2Skpi'
c = ROOT.TCanvas(c_name, c_name, 600, 400)
c.cd()
# Select required quantity
a_true = a_data['mjpk'][a_data['scaledmass'] < 5250.]
a_sim  = a_B02psi2skpi['m_pim_1_drop'] * 1000. # Convert to MeV
# Create and format histogram
h_data = ROOT.TH1F(c_name+'_data', c_name+'_data;B Mass [MeV/#it{c}^{2}];candidates', 200, 3650., 5150.)
h_sim  = ROOT.TH1F(c_name+'_sim' , c_name+'_sim ;B Mass [MeV/#it{c}^{2}];candidates', 200, 3650., 5150.)
# Fill histograms
map(h_data.Fill, a_true)
map(h_sim .Fill, a_sim)
# Normalise
h_data.Scale(1./h_data.Integral())
h_sim .Scale(1./h_sim .Integral())
## Make it pretty
h_data.SetTitle('Rapid Sim B Mass Distribution')
# Format for each case of x-axis
h_data.GetYaxis().SetTitleOffset(1.6)
y_max = 1.1*max(h_data.GetBinContent(h_data.GetMaximumBin()), h_sim.GetBinContent(h_sim.GetMaximumBin()))
y_min = 0.9*min(h_data.GetBinContent(h_data.GetMinimumBin()), h_sim.GetBinContent(h_sim.GetMinimumBin()))
h_data.GetYaxis().SetRangeUser(y_min, y_max)
# Format plotting style
h_data.SetLineColor(ROOT.kRed)
h_data.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
h_sim .SetLineColor(ROOT.kBlue)
h_sim .SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
# Remove stats boxes
h_data.SetStats(False)
h_sim .SetStats(False)
# Print
h_data.Draw('HIST')
h_sim .Draw('HISTsame')
# Create legend
leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
leg.AddEntry(h_data, 'B^{+} #rightarrow j/#psi #pi^{+}#pi^{-}K^{+}', 'L')
leg.AddEntry(h_sim , 'B^{0} #rightarrow #psi(2S)K^{+}#pi^{-}', 'L')
leg.SetLineColor(0)
leg.SetLineStyle(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw('same')
# Save
c.SaveAs('/home/s1305440/PPE_disk/project_stuff/data/rapdiSimComp_B02psi2SpiK.pdf')


print('   *** Plotting B+ -> J/psi pi pi K, K->pi misID ***')
# Initialise canvas
c_name = 'Bu2JpsipipiK_pi2k'
c = ROOT.TCanvas(c_name, c_name, 600, 400)
c.cd()
# Select required quantity
a_true = a_data['mppp']
a_sim  = a_Bu2JpsipipiK['Bp_0_M_Kp_02pip'] * 1000. # Convert to MeV
# Create and format histogram
h_data = ROOT.TH1F(c_name+'_data', c_name+'_data;B Mass [MeV/#it{c}^{2}];candidates', 200, 3650., 5400.)
h_sim  = ROOT.TH1F(c_name+'_sim' , c_name+'_sim ;B Mass [MeV/#it{c}^{2}];candidates', 200, 3650., 5400.)
# Fill histograms
map(h_data.Fill, a_true)
map(h_sim .Fill, a_sim)
# Normalise
h_data.Scale(1./h_data.Integral())
h_sim .Scale(1./h_sim .Integral())
## Make it pretty
h_data.SetTitle('Rapid Sim B Mass Distribution')
# Format for each case of x-axis
h_data.GetYaxis().SetTitleOffset(1.6)
y_max = 1.1*max(h_data.GetBinContent(h_data.GetMaximumBin()), h_sim.GetBinContent(h_sim.GetMaximumBin()))
y_min = 0.9*min(h_data.GetBinContent(h_data.GetMinimumBin()), h_sim.GetBinContent(h_sim.GetMinimumBin()))
h_data.GetYaxis().SetRangeUser(y_min, y_max)
# Format plotting style
h_data.SetLineColor(ROOT.kRed)
h_data.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
h_sim .SetLineColor(ROOT.kBlue)
h_sim .SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
# Remove stats boxes
h_data.SetStats(False)
h_sim .SetStats(False)
# Print
h_data.Draw('HIST')
h_sim .Draw('HISTsame')
# Create legend
leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
leg.AddEntry(h_data, 'B^{+} #rightarrow J/#psi #pi^{+}#pi^{-}K^{+}, K#rightarrow#pi misID data', 'L')
leg.AddEntry(h_sim , 'B^{+} #rightarrow J/#psi #pi^{+}#pi^{-}K^{+}, K#rightarrow#pi misID sim ', 'L')
leg.SetLineColor(0)
leg.SetLineStyle(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw('same')
# Save
c.SaveAs('/home/s1305440/PPE_disk/project_stuff/data/rapdiSimComp_Bu2JpsipipiK_k2pi.pdf')




print('   *** Plotting B+ -> J/psi pi pi pi ***')
# Initialise canvas
c_name = 'Bu2Jpsipipipi'
c = ROOT.TCanvas(c_name, c_name, 600, 400)
c.cd()
# Select required quantity
a_true = a_data['mppp']
a_sim  = a_Bu2JpsipipiK['Bp_0_M'] * 1000. # Convert to MeV
# Create and format histogram
h_data = ROOT.TH1F(c_name+'_data', c_name+'_data;B Mass [MeV/#it{c}^{2}];candidates', 200, 3650., 5400.)
h_sim  = ROOT.TH1F(c_name+'_sim' , c_name+'_sim ;B Mass [MeV/#it{c}^{2}];candidates', 200, 3650., 5400.)
# Fill histograms
map(h_data.Fill, a_true)
map(h_sim .Fill, a_sim)
# Normalise
h_data.Scale(1./h_data.Integral())
h_sim .Scale(1./h_sim .Integral())
## Make it pretty
h_data.SetTitle('Rapid Sim B Mass Distribution')
# Format for each case of x-axis
h_data.GetYaxis().SetTitleOffset(1.6)
y_max = 1.1*max(h_data.GetBinContent(h_data.GetMaximumBin()), h_sim.GetBinContent(h_sim.GetMaximumBin()))
y_min = 0.9*min(h_data.GetBinContent(h_data.GetMinimumBin()), h_sim.GetBinContent(h_sim.GetMinimumBin()))
h_data.GetYaxis().SetRangeUser(y_min, y_max)
# Format plotting style
h_data.SetLineColor(ROOT.kRed)
h_data.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
h_sim .SetLineColor(ROOT.kBlue)
h_sim .SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
# Remove stats boxes
h_data.SetStats(False)
h_sim .SetStats(False)
# Print
h_data.Draw('HIST')
h_sim .Draw('HISTsame')
# Create legend
leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
leg.AddEntry(h_data, 'B^{+} #rightarrow J/#psi #pi^{+}#pi^{-}K^{+}, K#rightarrow#pi misID', 'L')
leg.AddEntry(h_sim , 'B^{+} #rightarrow J/#psi #pi^{+}#pi^{-}#pi^{+}', 'L')
leg.SetLineColor(0)
leg.SetLineStyle(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw('same')
# Save
c.SaveAs('/home/s1305440/PPE_disk/project_stuff/data/rapdiSimComp_Bu2Jpsipipipi.pdf')




print('   *** Plotting B+ -> J/psi pi pi pi, K misID ***')
# Initialise canvas
c_name = 'Bu2Jpsipipipi_pi2k'
c = ROOT.TCanvas(c_name, c_name, 600, 400)
c.cd()
# Select required quantity
a_true = a_data['scaledmass']
a_sim  = a_Bu2Jpsipipipi['Bp_0_M_pip_12Kp'] * 1000. # Convert to MeV
# Create and format histogram
h_data = ROOT.TH1F(c_name+'_data', c_name+'_data;B Mass [MeV/#it{c}^{2}];candidates', 200, 5200., 5600.)
h_sim  = ROOT.TH1F(c_name+'_sim' , c_name+'_sim ;B Mass [MeV/#it{c}^{2}];candidates', 200, 5200., 5600.)
# Fill histograms
map(h_data.Fill, a_true)
map(h_sim .Fill, a_sim)
# Normalise
h_data.Scale(1./h_data.Integral())
h_sim .Scale(1./h_sim .Integral())
## Make it pretty
h_data.SetTitle('Rapid Sim B Mass Distribution')
# Format for each case of x-axis
h_data.GetYaxis().SetTitleOffset(1.6)
y_max = 1.1*max(h_data.GetBinContent(h_data.GetMaximumBin()), h_sim.GetBinContent(h_sim.GetMaximumBin()))
y_min = 0.9*min(h_data.GetBinContent(h_data.GetMinimumBin()), h_sim.GetBinContent(h_sim.GetMinimumBin()))
h_data.GetYaxis().SetRangeUser(y_min, y_max)
# Format plotting style
h_data.SetLineColor(ROOT.kRed)
h_data.SetFillColorAlpha(ROOT .kRed - 10, 0.7)
h_sim .SetLineColor(ROOT.kBlue)
h_sim .SetFillColorAlpha(ROOT.kBlue - 10, 0.7)
# Remove stats boxes
h_data.SetStats(False)
h_sim .SetStats(False)
# Print
h_data.Draw('HIST')
h_sim .Draw('HISTsame')
# Create legend
leg = ROOT.TLegend(0.6,0.75,0.9,0.9)
leg.AddEntry(h_data, 'B^{+} #rightarrow J/#psi #pi^{+}#pi^{-}K^{+}', 'L')
leg.AddEntry(h_sim , 'B^{+} #rightarrow J/#psi #pi^{+}#pi^{-}#pi^{+}, #pi#rightarrow K misID', 'L')
leg.SetLineColor(0)
leg.SetLineStyle(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.Draw('same')
# Save
c.SaveAs('/home/s1305440/PPE_disk/project_stuff/data/rapdiSimComp_Bu2Jpsipipipi_pi2k.pdf')