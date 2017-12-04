#include <string>
#include <iostream>
#include "TROOT.h"

// Arrays for fitting values
void interpRes(){
	double m[3]    = {3376., 3686.1  , 3872.2  };
	double em[3]   = {   0.,    0.   ,    0.   };
	double eres[3] = {   0.,    0.011,    0.013};
	double res[3]  = {   0.,    2.057,    2.721};
	// Fill graph object
	TGraphErrors* gr = new TGraphErrors(3, m, res, em, eres);
	// Initiate function for interpolating
	TF1* fun = new TF1("fun","pow([0]*TMath::Abs((x-3376)),0.5)", 0.1, 100e3);
	// 

	// Make it pretty
	fun->SetLineStyle(2);
	fun->SetLineColor(2);
	fun->SetParameters(1,333);
	// Interpolating
	gr->Fit(fun, "","",3373., 4000);
	// Print results
	std::cout << "Best fit: "    << fun->Eval(3823.1) << std::endl;
	std::cout << "Interpolate: " <<  gr->Eval(3823.1) << std::endl;


	TCanvas *c1 = new TCanvas("c1","pretty plot",200,10,700,500);

	gr->SetTitle("Energy Resolution Interpolation");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("A*E");

	fun->Draw("Lsame");

	c1->Update();

	c1->SaveAs("fitted_interp.pdf");
}