#ifndef FITFUNCTIONS_HH
#define FITFUNCTIONS_HH
//#define PWAVE
#define DWAVE

// -----=====-----=====----- //
// DATA                      //
// -----=====-----=====----- //
static TFile* file_data = new TFile("root/long_cdc_lpn_bs_wcap.root","READ");
static TH2F* data_mass_momtrans;
static TH1F* data_mass;
static TH1F* data_momtrans;
// -----=====-----=====----- //
// SIMULATED ACCEPTANCE      //
// -----=====-----=====----- //
static TFile* file_acc = new TFile("root/long_cdc_lpn_wcap.root","READ");
static TH2F* generated_mass_momtrans;
static TH1F* generated_mass;
static TH1F* generated_momtrans;
static TH2F* phasespace_mass_momtrans;
static TH1F* phasespace_mass;
static TH1F* phasespace_momtrans;
// -----=====-----=====----- //
// OTHER FUNCTIONS           //
// -----=====-----=====----- //
Double_t eval_q(Double_t m_lp, Double_t cos){
	if(cos<-1.0||1.0<cos){
		return -1;
	}
	Float_t sin = sqrt(1-cos*cos);
	Double_t p_k = 1.003;
	Double_t E_k = sqrt(kpMass*kpMass+p_k*p_k);
	Double_t E_cm = sqrt( pow(E_k+ThreeHeMass,2) - pow(p_k,2) );
	if(m_lp<lMass+pMass||E_cm-nMass<m_lp){
		return -1;
	}
	Double_t beta = p_k/(E_k+ThreeHeMass);
	Double_t gamma = (E_k+ThreeHeMass)/E_cm;
	Double_t p_n = (1.0/2.0/E_cm) * sqrt(pow(m_lp,4) - 2.0*(pow(E_cm,2)+pow(nMass,2))*pow(m_lp,2) + pow((pow(E_cm,2)-pow(nMass,2)),2));
	Double_t E_lp_cm = sqrt(nMass*nMass+p_n*p_n);
	Double_t val = sqrt( pow(p_k - gamma*beta*E_lp_cm - gamma*p_n*cos,2) + pow(p_n*sin,2) );
	return val;
}
// -----=====-----=====----- //
// FUNCTIONS                 //
// -----=====-----=====----- //
// ----- Kpp shape ----- //
Double_t kpp(Double_t* x, Double_t* par){
  // x[0] : m, x[1] : q
  // par[0] : Kpp_Factor0
  // par[1] : Kpp_Mass
  // par[2] : Kpp_Width
  // par[3] : Kpp_Q
  Double_t val_1 = (par[2]/2.0) / ( pow(x[0]-par[1],2) + pow(par[2]/2.0,2) );
#if defined(PWAVE)
  Double_t val_2 = x[1]*x[1]*TMath::Exp( -pow(x[1],2)/pow(par[3],2) );
#elif defined(DWAVE)
  Double_t val_2 = x[1]*x[1]*TMath::Exp( -pow(x[1],2)/pow(par[3],2) );
#else
  Double_t val_2 = TMath::Exp( -pow(x[1],2)/pow(par[3],2) );
#endif
  // Set value
  Double_t val = par[0] * val_1 * val_2;
  val *= phasespace_mass_momtrans->GetBinContent(phasespace_mass->FindBin(x[0]),phasespace_momtrans->FindBin(x[1]));
  val *= phasespace_mass->GetBinWidth(1)*1000.0;
  val *= phasespace_momtrans->GetBinWidth(1)*1000.0;
  return val; 
}
// -----=====-----=====----- //
// DRAW FUNCTION             //
// -----=====-----=====----- //
void SetHist(TH1* hist){
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(1);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.0);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetLabelFont(132);
  hist->GetYaxis()->SetLabelFont(132);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleFont(132);
  hist->GetYaxis()->SetTitleFont(132);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->SetStats(false);
  hist->SetTitle("");
}
// -----=====-----=====----- //
// DRAW FUNCTION             //
// -----=====-----=====----- //
//   CalcConfidenceInterval(        hist,    fitfunc,    cov. matrix,        conf level)
void CalcConfidenceInterval(TH1* hfit, TF1* f, Double_t* matr,  Double_t cl=0.95){
  Int_t npar = f->GetNpar();
  Double_t *grad = new Double_t[npar];
  Double_t *sum_vector = new Double_t[npar];
  Double_t x[3];

  Int_t hxfirst = hfit->GetXaxis()->GetFirst();
  Int_t hxlast  = hfit->GetXaxis()->GetLast();
  Int_t hyfirst = hfit->GetYaxis()->GetFirst();
  Int_t hylast  = hfit->GetYaxis()->GetLast();
  Int_t hzfirst = hfit->GetZaxis()->GetFirst();
  Int_t hzlast  = hfit->GetZaxis()->GetLast();

  TAxis *xaxis  = hfit->GetXaxis();
  TAxis *yaxis  = hfit->GetYaxis();
  TAxis *zaxis  = hfit->GetZaxis();
  Double_t t = TMath::StudentQuantile(0.5 + cl/2, f->GetNDF());
  Double_t chidf = TMath::Sqrt(f->GetChisquare()/f->GetNDF());
  Double_t c=0;
  for (Int_t binz=hzfirst; binz<=hzlast; binz++) {
    x[2]=zaxis->GetBinCenter(binz);
    for (Int_t biny=hyfirst; biny<=hylast; biny++) {
      x[1]=yaxis->GetBinCenter(biny);
      for (Int_t binx=hxfirst; binx<=hxlast; binx++) {
        x[0]=xaxis->GetBinCenter(binx);
        f->GradientPar(x, grad);
        for (Int_t irow=0; irow<npar; irow++){
          sum_vector[irow]=0;
          for (Int_t icol=0; icol<npar; icol++)
            sum_vector[irow]+=matr[irow*npar+icol]*grad[icol];
        }
        c = 0;
        for (Int_t i=0; i<npar; i++){
          c+=grad[i]*sum_vector[i];
        }
        c=TMath::Sqrt(c);
        hfit->SetBinContent(binx, biny, binz, f->EvalPar(x));
        hfit->SetBinError(binx, biny, binz, c*t*chidf);
      }
    }
  }
  delete [] grad;
  delete [] sum_vector;
}
#endif
