#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include "GlobalVariables.h"
#include "FitFunctions.h"
#include <iostream>
#include <fstream>

// -----=====-----=====----- //
// PARAMETERS                //
// -----=====-----=====----- //
static const Int_t npar =  4;
static Int_t fixparams[npar] = {0};
static Double_t params[npar] = {0};
static Double_t errors[npar] = {0};
static TString parname[npar];
static TString initparam = "param/initial_dwave.txt";	
static TString outparam  = "param/result_dwave.txt";	
static TString outname   = "fig/fit_result_dwave.pdf";
static TString out_root_file_name= "root/fit_result_dwave.root";

static Double_t mass_ll;
static Double_t mass_ul;
static Int_t mass_nbin;
static Double_t momtrans_ll;
static Double_t momtrans_ul;
static Int_t momtrans_nbin;

// -----=====-----=====----- //
// CONSTANT VALUES           //
// -----=====-----=====----- //
const Double_t analysis_efficiency  = 0.5;
const Double_t selection_efficiency = 0.8;
const Double_t facility_efficiency  = 0.9;
const Double_t week = 4.;
const Double_t luminosity_for_aweek = 1.9e3; // /ub
const Double_t cs_bs = 9.29; // ub

// -----=====-----=====----- //
// FITTING                   //
// -----=====-----=====----- //
void Fitting(TString name){


	TF2* f_fit = new TF2("f_fit",kpp,mass_ll,mass_ul,momtrans_ll,momtrans_ul,npar);
	
	for(int i=0; i<npar; i++){
		f_fit->SetParName(i,parname[i].Data());
		f_fit->SetParameter(i,params[i]);
		f_fit->SetParLimits(i,params[i]*0.01,params[i]*5);
		if(fixparams[i]){
			f_fit->FixParameter(i,params[i]);
		}
	}
	f_fit->SetNpx(mass_nbin);
	f_fit->SetNpy(momtrans_nbin);

	// Fitting
	data_mass_momtrans->Fit("f_fit","VLN0");
	for(int ipar=0; ipar<npar; ipar++){
		params[ipar] = (TVirtualFitter::GetFitter())->GetParameter(ipar);
		errors[ipar] = (TVirtualFitter::GetFitter())->GetParError(ipar);
	}

  // Drawing
  TCanvas* canvas = new TCanvas("canvas_fit","canvas_fit",2100,1800);
	canvas->cd();
	canvas->Divide(2,2);
  for(Int_t i_pad=1; i_pad<=4; ++i_pad){
    TPad* pad = (TPad*)canvas->cd(i_pad);
    pad->SetGrid();
    pad->SetTicks(); 
    pad->SetLeftMargin(0.16);
    pad->SetRightMargin(0.03);
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.16);         
  }

  canvas->cd(1);
  data_mass_momtrans->Draw("col");

  canvas->cd(2);
  TH2F* fit_mass_momtrans = (TH2F*)f_fit->CreateHistogram();
  fit_mass_momtrans->Draw("col");

  canvas->cd(3);
  data_mass->Draw("hist c");
  TH1F* fit_mass = (TH1F*)fit_mass_momtrans->ProjectionX();
  fit_mass->SetLineColor(kBlack);
  fit_mass->SetLineWidth(2);
  fit_mass->SetLineStyle(2);
  fit_mass->SetFillStyle(0);
  fit_mass->Draw("hist l same");

  canvas->cd(4);
  data_momtrans->Draw("hist c");
  TH1F* fit_momtrans = (TH1F*)fit_mass_momtrans->ProjectionY();
  fit_momtrans->SetLineColor(kBlack);
  fit_momtrans->SetLineWidth(2);
  fit_momtrans->SetLineStyle(2);
  fit_momtrans->SetFillStyle(0);
  fit_momtrans->Draw("hist l same");

  canvas->Print(outname.Data());

  TFile* out_root_file = new TFile(out_root_file_name.Data(),"RECREATE");
  out_root_file->cd();
  TH2F* th2_copy = new TH2F();
  th2_copy = (TH2F*)fit_mass_momtrans->Clone("lp_mass_vs_momtrans_fit");
  out_root_file->Write();
  out_root_file->Close();
}

// -----=====-----=====----- //
// MAIN FUNCTION             //
// -----=====-----=====----- //
int main(int argc, char* argv[])
{	
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetTitleOffset(1.0,"XYZ");

	if(argc>1){
		initparam = argv[1];
		if(argc>2){
			outparam = argv[2];
			if(argc>3){
				outname = argv[3];
			}
		}
	}

	// -----=====-----=====----- //
	// INITIAL PARAMETERS        //
	// -----=====-----=====----- //
	std::ifstream infile(initparam.Data());
	if(infile.fail()){
		std::cerr << "Failed to open file." << std::endl;
		std::cerr << initparam.Data() << std::endl;
		return 0;
	}
	{
		std::string str;
		int ipar=0;
		while(getline(infile,str)){
			char tmpchar[20];
			sscanf(str.c_str(),"%d, %lf : %s",&fixparams[ipar],&params[ipar],tmpchar);
			parname[ipar] = tmpchar;
			std::cout << parname[ipar].Data() << " >>> " << fixparams[ipar] << " , " << params[ipar] << std::endl;
			ipar++;
		}
	}
	// -----=====-----=====----- //
	// DATA                      //
	// -----=====-----=====----- //
  TH1F* h1_scaling_lpn_bs = (TH1F*)file_data->Get("beam_momentum_lab_gene");
  Double_t scaling_lpn_bs = h1_scaling_lpn_bs->GetEntries();
	data_mass_momtrans = (TH2F*)file_data->Get("lp_mass_vs_momtrans_kfit_lp");
  data_mass_momtrans->Scale(luminosity_for_aweek*week*cs_bs*analysis_efficiency*selection_efficiency*facility_efficiency/scaling_lpn_bs);
	data_mass = (TH1F*)data_mass_momtrans->ProjectionX("data_mass");
	data_momtrans = (TH1F*)data_mass_momtrans->ProjectionY("data_momtrans");
  mass_ll = data_mass->GetBinCenter(data_mass->GetXaxis()->GetFirst()) - data_mass->GetBinWidth(data_mass->GetXaxis()->GetFirst())/2.;
  mass_ul = data_mass->GetBinCenter(data_mass->GetXaxis()->GetLast()) + data_mass->GetBinWidth(data_mass->GetXaxis()->GetLast())/2.;
  mass_nbin = data_mass->GetNbinsX();
  momtrans_ll = data_momtrans->GetBinCenter(data_momtrans->GetXaxis()->GetFirst()) - data_momtrans->GetBinWidth(data_momtrans->GetXaxis()->GetFirst())/2.;
  momtrans_ul = data_momtrans->GetBinCenter(data_momtrans->GetXaxis()->GetLast()) + data_momtrans->GetBinWidth(data_momtrans->GetXaxis()->GetLast())/2.;
  momtrans_nbin = data_momtrans->GetNbinsX();
	// -----=====-----=====----- //
	// SIMULATED ACCEPTANCE      //
	// -----=====-----=====----- //
	generated_mass_momtrans = (TH2F*)file_acc->Get("lp_mass_vs_momtrans_gene");
	generated_mass_momtrans->Scale(1.0/generated_mass_momtrans->GetEntries());
	generated_mass = (TH1F*)generated_mass_momtrans->ProjectionX("generated_mass");
	generated_momtrans = (TH1F*)generated_mass_momtrans->ProjectionY("generated_momtrans");
	phasespace_mass_momtrans = (TH2F*)file_acc->Get("lp_mass_vs_momtrans_acce_lp");
	phasespace_mass_momtrans->Scale(1.0/generated_mass_momtrans->GetEntries());
	phasespace_mass = (TH1F*)phasespace_mass_momtrans->ProjectionX("phasespace_mass");
	phasespace_momtrans = (TH1F*)phasespace_mass_momtrans->ProjectionY("phasespace_momtrans");
	// -----=====-----=====----- //

  TCanvas* canvas = new TCanvas("canvas","canvas",3200,900);
	canvas->cd();
	canvas->Divide(3,1);
  for(Int_t i_pad=1; i_pad<=3; ++i_pad){
    TPad* pad = (TPad*)canvas->cd(i_pad);
    pad->SetGrid();
    pad->SetTicks(); 
    pad->SetLeftMargin(0.16);
    pad->SetRightMargin(0.03);
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.16);         
  }

	// data
	canvas->cd(1);
	data_mass_momtrans->Draw("col");
	// phasespace
	canvas->cd(2);
	generated_mass_momtrans->Draw("col");
	// phasespace
	canvas->cd(3);
	phasespace_mass_momtrans->Draw("col");
	canvas->Print("fig/fit_1.pdf");

	Fitting(outname.Data());

}
