// Analyzer.cc

#include "Analyzer.hh"

Analyzer::Analyzer()
{
  root_file_name_ = "tmp.root";
  CreateHistograms();
}

Analyzer::Analyzer(std::string file_name)
{
  root_file_name_ = file_name;
  CreateHistograms();
}

Analyzer::~Analyzer()
{
  root_file_->Write();
  root_file_->Close();
}

void Analyzer::CreateHistograms()
{
  root_file_ = new TFile(root_file_name_.data(),"recreate");
  root_file_->cd();


  new TH1F("lp_mass","invariant mass of lp;#font[12]{m}_{#Lambda#font[12]{p}} (GeV/#font[12]{c}^{2});Counts",100,2.,3.);
  new TH1F("lp_momtrans","momentum transfer of lp;#font[12]{q}_{#Lambda#font[12]{p}} (GeV/#font[12]{c});Counts",100,0.,2.);

  new TH2F("lp_mass_vs_lp_momtrans","mass vs. momentum transfer of lp;#font[12]{m}_{#Lambda#font[12]{p}} (GeV/#font[12]{c}^{2});#font[12]{q}_{#Lambda#font[12]{p}} (GeV/#font[12]{c});Counts",100,2.,3.,100,0.,2.);

  new TH1F("ppim_mass","invariant mass of ppim;#font[12]{m}_{#font[12]{p}#pi^{-}} (GeV/#font[12]{c}^{2});Counts",100,1.0,1.2);


  std::string ml_name[8] = {"","_ml_0","_ml_p1","_ml_m1",
    "_cos_select","_cos_select_ml_0","_cos_select_ml_p1","_cos_select_ml_m1"};

  for(Int_t i_ml=0; i_ml<8; ++i_ml){
    new TH1F(Form("reference_dot_lambda_spin%s",ml_name[i_ml].data()),"refernece dot lambda spin;#font[12]{n}_{ref}#upoint#font[12]{S}_{#Lambda};Counts",100,-1.001,1.001);
    new TH1F(Form("lambda_spin_dot_proton_spin%s",ml_name[i_ml].data()),"lambda spin dot proton spin;#font[12]{S}_{#Lambda}#upoint#font[12]{S}_{#font[12]{p}};Counts",100,-1.001,1.001);
    new TH1F(Form("reference_dot_lambda_momentum%s",ml_name[i_ml].data()),"refernece dot lambda momentum;#font[12]{n}_{ref}#upoint#font[12]{p}_{#Lambda}^{#Lambda#font[12]{p} rest};Counts",100,-1.001,1.001);
    new TH1F(Form("boost_lambda_dot_lambda_spin%s",ml_name[i_ml].data()),"lambda boost vector dot lambda spin;#font[12]{p}_{#Lambda}^{#Lambda#font[12]{p} rest}#upoint#font[12]{S}_{#Lambda};Counts",100,-1.001,1.001);
    new TH1F(Form("boost_lambda_dot_proton_from_lambda_momentum%s",ml_name[i_ml].data()),"lambda boost vector dot proton from lambda momentum;#font[12]{p}_{#Lambda}^{#Lambda#font[12]{p} rest}#upoint#font[12]{p}_{#font[12]{p} from #Lambda};Counts",100,-1.001,1.001);
    new TH1F(Form("lambda_spin_dot_proton_from_lambda_momentum%s",ml_name[i_ml].data()),"lambda spin dot proton from lambda momentum;#font[12]{S}_{#Lambda}#upoint#font[12]{p}_{#font[12]{p}}^{#Lambda rest};Counts",100,-1.001,1.001);
    new TH1F(Form("lambda_spin_dot_proton_from_lambda_spin%s",ml_name[i_ml].data()),"lambda spin dot proton from lambda momentum;#font[12]{S}_{#Lambda}#upoint#font[12]{S}_{#font[12]{p} from #Lambda};Counts",100,-1.001,1.001);
    new TH1F(Form("boost_lambda_dot_proton_from_lambda_spin%s",ml_name[i_ml].data()),"lambda boost vector dot proton from lambda spin;#font[12]{p}_{#Lambda}^{#Lambda#font[12]{p} rest}#upoint#font[12]{S}_{#font[12]{p} form #Lambda};Counts",100,-1.001,1.001);
    new TH1F(Form("proton_spin_dot_proton_from_lambda_spin%s",ml_name[i_ml].data()),"proton spin dot proton from lambda spin;#font[12]{S}_{#font[12]{p}}#upoint#font[12]{S}_{#font[12]{p} from #Lambda};Counts",100,-1.001,1.001);
    new TH1F(Form("proton_spin_dot_proton_from_lambda_momentum%s",ml_name[i_ml].data()),"proton spin dot proton from lambda momentum;#font[12]{S}_{#font[12]{p}}#upoint#font[12]{p}_{#font[12]{p} from #Lambda}^{#Lambda rest};Counts",100,-1.001,1.001);
  }
}

void Analyzer::Analysis(KppGenerator* kpp_generator, LambdaDecay* lambda_decay)
{
  TLorentzVector lv_lambda = kpp_generator->LvLambda();
  TLorentzVector lv_proton = kpp_generator->LvProton();
  TLorentzVector lv_neutron = kpp_generator->LvNeutron();

  TLorentzVector lv_lp = lv_lambda + lv_proton;

  TLorentzVector lv_proton_from_lambda = lambda_decay->LvProton();
  TLorentzVector lv_pim_from_lambda = lambda_decay->LvPim();

  TLorentzVector lv_ppim_from_lambda = lv_proton_from_lambda + lv_pim_from_lambda;

  TVector3 vec_boost_lp = lv_lp.BoostVector();
  TLorentzVector clv_lambda_lp = lv_lambda;
  clv_lambda_lp.Boost(-vec_boost_lp);
  TLorentzVector clv_proton_from_lambda_lp = lv_proton_from_lambda;
  clv_proton_from_lambda_lp.Boost(-vec_boost_lp);
  TVector3 vec_boost_lambda = clv_lambda_lp.BoostVector();
  TLorentzVector clv_proton_from_lambda_lambda = clv_proton_from_lambda_lp;
  clv_proton_from_lambda_lambda.Boost(-vec_boost_lambda);

  Fill("lp_mass",lv_lp.M());
  Fill("lp_momtrans",lv_lp.P());
  Fill("lp_mass_vs_lp_momtrans",lv_lp.M(),lv_lp.P());


  Int_t ml_decay = kpp_generator->MLDecay();
  TVector3 vec_reference_direction = kpp_generator->VecReference().Unit();
  TVector3 vec_lambda_spin_direction = kpp_generator->VecLambdaSpin().Unit();
  TVector3 vec_proton_spin_direction = kpp_generator->VecProtonSpin().Unit();
  TVector3 vec_proton_from_lambda_spin_direction = lambda_decay->VecProtonSpin().Unit();

  TVector3 vec_lambda_momentum_direction = clv_lambda_lp.Vect().Unit();
  TVector3 vec_proton_from_lambda_momentum_direction = clv_proton_from_lambda_lambda.Vect().Unit();

  Double_t reference_dot_lambda_spin = vec_reference_direction.Dot(vec_lambda_spin_direction);
  Double_t lambda_spin_dot_proton_spin = vec_lambda_spin_direction.Dot(vec_proton_spin_direction);
  Double_t reference_dot_lambda_momentum = vec_reference_direction.Dot(vec_lambda_momentum_direction);
  Double_t boost_lambda_dot_lambda_spin = vec_boost_lambda.Unit().Dot(vec_lambda_spin_direction);
  Double_t boost_lambda_dot_proton_from_lambda_momentum = vec_boost_lambda.Unit().Dot(vec_proton_from_lambda_momentum_direction);
  Double_t lambda_spin_dot_proton_from_lambda_momentum = vec_lambda_spin_direction.Dot(vec_proton_from_lambda_momentum_direction);
  Double_t lambda_spin_dot_proton_from_lambda_spin = vec_lambda_spin_direction.Dot(vec_proton_from_lambda_spin_direction);
  Double_t boost_lambda_dot_proton_from_lambda_spin = vec_boost_lambda.Unit().Dot(vec_proton_from_lambda_spin_direction);
  Double_t proton_spin_dot_proton_from_lambda_spin = vec_proton_spin_direction.Dot(vec_proton_from_lambda_spin_direction);
  Double_t proton_spin_dot_proton_from_lambda_momentum = vec_proton_spin_direction.Dot(vec_proton_from_lambda_momentum_direction);

  Bool_t fill_flags[8] = {false};
  fill_flags[0] = true;
  if(ml_decay==0){
    fill_flags[1] = true;
  }
  if(ml_decay==1){
    fill_flags[2] = true;
  }
  if(ml_decay==-1){
    fill_flags[3] = true;
  }

  Double_t cos_select_ll = -0.5;
  Double_t cos_select_ul =  0.5;

  if(cos_select_ll<boost_lambda_dot_proton_from_lambda_momentum&&boost_lambda_dot_proton_from_lambda_momentum<cos_select_ul){
    fill_flags[4] = true;
    if(ml_decay==0){
      fill_flags[5] = true;
    }
    if(ml_decay==1){
      fill_flags[6] = true;
    }
    if(ml_decay==-1){
      fill_flags[7] = true;
    }
  }


  std::string ml_name[8] = {"","_ml_0","_ml_p1","_ml_m1",
    "_cos_select","_cos_select_ml_0","_cos_select_ml_p1","_cos_select_ml_m1"};

  for(Int_t i_ml=0; i_ml<8; ++i_ml){
    if(fill_flags[i_ml]){
      Fill("reference_dot_lambda_spin"+ml_name[i_ml],reference_dot_lambda_spin);
      Fill("lambda_spin_dot_proton_spin"+ml_name[i_ml],lambda_spin_dot_proton_spin);
      Fill("reference_dot_lambda_momentum"+ml_name[i_ml],reference_dot_lambda_momentum);
      Fill("boost_lambda_dot_lambda_spin"+ml_name[i_ml],boost_lambda_dot_lambda_spin);
      Fill("boost_lambda_dot_proton_from_lambda_momentum"+ml_name[i_ml],boost_lambda_dot_proton_from_lambda_momentum);
      Fill("lambda_spin_dot_proton_from_lambda_momentum"+ml_name[i_ml],lambda_spin_dot_proton_from_lambda_momentum);
      Fill("lambda_spin_dot_proton_from_lambda_spin"+ml_name[i_ml],lambda_spin_dot_proton_from_lambda_spin);
      Fill("boost_lambda_dot_proton_from_lambda_spin"+ml_name[i_ml],boost_lambda_dot_proton_from_lambda_spin);
      Fill("proton_spin_dot_proton_from_lambda_spin"+ml_name[i_ml],proton_spin_dot_proton_from_lambda_spin);
      Fill("proton_spin_dot_proton_from_lambda_momentum"+ml_name[i_ml],proton_spin_dot_proton_from_lambda_momentum);
    }
  }

  return;
}

void Analyzer::Fill(std::string hist_name, Double_t val)
{
  TH1* hist = (TH1*)(root_file_->Get(hist_name.data()));
  if(hist){
    hist->Fill(val);
  }
  return;
}

void Analyzer::Fill(std::string hist_name, Double_t val1, Double_t val2)
{
  TH1* hist = (TH1*)(root_file_->Get(hist_name.data()));
  if(hist){
    hist->Fill(val1,val2);
  }
  return;
}

void Analyzer::Draw(std::string hist_name, std::string option="")
{
  TH1* hist = (TH1*)(root_file_->Get(hist_name.data()));
  if(hist){
    hist->SetLineColor(kBlack);
    hist->Draw(option.data());
    hist->SetMinimum(0.);
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(132);
    hist->GetXaxis()->SetLabelFont(132);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleFont(132);
    hist->GetYaxis()->SetLabelFont(132);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.05);
  }
  return;
}

void Analyzer::PrintHistograms(std::string pdf_name)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* title = new TCanvas("tltle","title",3200,1800);
  TLatex* text = new TLatex();
  text->SetTextAlign(22);
  text->SetTextFont(132);
  text->SetTextSize(0.05);

  TCanvas* canvas = new TCanvas("canvas","canvas",3200,1800);
  canvas->Divide(3,2);
  for(Int_t i_pad=1; i_pad<=6; ++i_pad){
    TPad* pad = (TPad*)canvas->cd(i_pad);
    pad->SetLeftMargin(0.16);
    pad->SetRightMargin(0.03);
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.16);         
  }
  canvas->Print(std::string(pdf_name+"[").data());


  // ---------------------------------------
  canvas->cd(1);
  Draw("lp_mass");
  canvas->cd(2);
  Draw("lp_momtrans");
  canvas->cd(3);
  Draw("lp_mass_vs_lp_momtrans","col");
  canvas->cd(4);
  Draw("ppim_mass");
  canvas->Print(pdf_name.data());
  for(Int_t i_pad=1; i_pad<=6; ++i_pad){
    canvas->cd(i_pad);
    gPad->Clear();
  }
  // ---------------------------------------

  std::string ml_name[8] = {"","_ml_0","_ml_p1","_ml_m1",
    "_cos_select","_cos_select_ml_0","_cos_select_ml_p1","_cos_select_ml_m1"};
  for(Int_t i_ml=0; i_ml<8; ++i_ml){
    title->cd();
    text->DrawLatex(0.5,0.5,ml_name[i_ml].data());
    title->Print(pdf_name.data());
    title->Clear();
    // ---------------------------------------
    canvas->cd(1);
    Draw("reference_dot_lambda_spin"+ml_name[i_ml]);
    canvas->cd(2);
    Draw("lambda_spin_dot_proton_spin"+ml_name[i_ml]);
    canvas->cd(3);
    Draw("reference_dot_lambda_momentum"+ml_name[i_ml]);
    canvas->cd(4);
    Draw("boost_lambda_dot_lambda_spin"+ml_name[i_ml]);
    canvas->cd(5);
    Draw("boost_lambda_dot_proton_from_lambda_momentum"+ml_name[i_ml]);
    canvas->cd(6);
    Draw("lambda_spin_dot_proton_from_lambda_momentum"+ml_name[i_ml]);
    canvas->Print(pdf_name.data());
    for(Int_t i_pad=1; i_pad<=6; ++i_pad){
      canvas->cd(i_pad);
      gPad->Clear();
    }
    // ---------------------------------------

    // ---------------------------------------
    canvas->cd(1);
    Draw("lambda_spin_dot_proton_from_lambda_spin"+ml_name[i_ml]);
    canvas->cd(2);
    Draw("boost_lambda_dot_proton_from_lambda_spin"+ml_name[i_ml]);
    canvas->cd(3);
    Draw("proton_spin_dot_proton_from_lambda_spin"+ml_name[i_ml]);
    canvas->cd(4);
    Draw("proton_spin_dot_proton_from_lambda_momentum"+ml_name[i_ml]);
    canvas->Print(pdf_name.data());
    for(Int_t i_pad=1; i_pad<=6; ++i_pad){
      canvas->cd(i_pad);
      gPad->Clear();
    }
    // ---------------------------------------
  }

  canvas->Print(std::string(pdf_name+"]").data());
  return;
}
