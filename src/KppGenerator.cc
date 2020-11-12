// KppGenerator.cc

#include "KppGenerator.hh"

KppGenerator::KppGenerator()
{
  kpp_spin_ = 0;
  kpp_parity_ = -1;

  lv_target_ = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  lv_target_.SetVectM(TVector3(0.,0.,0.),kThreeHeMass);

	f_kpp_mass_ = new TF1("f_kpp_mass","([1]*[1]/4.) / ( (x-[0])*(x-[0]) + [1]*[1]/4)",2.,3.);
  f_kpp_mass_->SetParameter(0,2.327518e+00);
  f_kpp_mass_->SetParameter(1,1.003277e-01);
	f_kpp_momtrans_ = new TF1("f_kpp_momtrans","exp(-x*x/[0]/[0])",0.,2.);
  f_kpp_momtrans_->SetParameter(0,3.825931e-01);

  random_ = new TRandom3(0);
}

KppGenerator::KppGenerator(Int_t kpp_spin, Int_t kpp_parity)
{
  kpp_spin_ = kpp_spin;
  kpp_parity_ = kpp_parity;

  lv_target_ = TLorentzVector(0.0, 0.0, 0.0, 0.0);
  lv_target_.SetVectM(TVector3(0.,0.,0.),kThreeHeMass);

	f_kpp_mass_ = new TF1("f_kpp_mass","([1]*[1]/4.) / ( (x-[0])*(x-[0]) + [1]*[1]/4)",2.,3.);
  f_kpp_mass_->SetParameter(0,2.327518e+00);
  f_kpp_mass_->SetParameter(1,1.003277e-01);
	f_kpp_momtrans_ = new TF1("f_kpp_momtrans","exp(-x*x/[0]/[0])",0.,2.);
  f_kpp_momtrans_->SetParameter(0,3.825931e-01);

  random_ = new TRandom3(0);
}

KppGenerator::~KppGenerator()
{
  delete f_kpp_mass_;
  delete f_kpp_momtrans_;
  delete random_;
}

void KppGenerator::Generate(TLorentzVector& lv_beam)
{
  TLorentzVector lv_cm = lv_beam + lv_target_;
  generator_.SetDecay(lv_cm, 3, cm_masses_);

  while(true){
    Double_t weight_gen = generator_.Generate();
    Double_t weight_ran = random_->Uniform()*kMaxWeightThreeBody;
    if(weight_gen>weight_ran){
      TLorentzVector lv_lambda = *(generator_.GetDecay(0));
      TLorentzVector lv_proton = *(generator_.GetDecay(1));
      Double_t mass_lp = (lv_lambda+lv_proton).M();
      Double_t momtrans_lp = (lv_lambda+lv_proton).P();
      Double_t val_func = f_kpp_mass_->Eval(mass_lp);
      val_func *= f_kpp_momtrans_->Eval(momtrans_lp);
      Double_t val_ran = random_->Uniform();
      if(val_func>val_ran){ break; }
    }
  }

  lv_kpp_ = *(generator_.GetDecay(0)) + *(generator_.GetDecay(1));
  lv_neutron_ = *(generator_.GetDecay(2));

  KppZeroMinusDecay();
}

void KppGenerator::KppZeroMinusDecay()
{
  // L-direction
  Double_t cos_theta_l = random_->Uniform(-1.,1.);
  Double_t sin_theta_l = sqrt(1. - cos_theta_l*cos_theta_l);
  Double_t phi_l = random_->Uniform(-TMath::Pi(),TMath::Pi());
  vec_reference_direction_ = TVector3(sin_theta_l*cos(phi_l),sin_theta_l*sin(phi_l),cos_theta_l);

  // ml (of P-wave decay)
  Double_t random_ml = random_->Uniform(0.,3.);
  if(random_ml<1.){
    ml_ = -1;
  }
  else if(random_ml<2.){
    ml_ = 0;
  }
  else{
    ml_ = 1;
  }

  // momentum directions (P-wave decay)
  Double_t cos_theta_decay = -999.;
  if(ml_==-1){
      while(true){
        Double_t val_cos_theta = random_->Uniform(-1.,1.);
        Double_t val_random = random_->Uniform();
        if(1.-val_cos_theta*val_cos_theta>val_random){
          cos_theta_decay = val_cos_theta;
          break;
        }
      }
  }
  if(ml_==0){
      while(true){
        Double_t val_cos_theta = random_->Uniform(-1.,1.);
        Double_t val_random = random_->Uniform();
        if(val_cos_theta*val_cos_theta>val_random){
          cos_theta_decay = val_cos_theta;
          break;
        }
      }
  }
  if(ml_==1){
      while(true){
        Double_t val_cos_theta = random_->Uniform(-1.,1.);
        Double_t val_random = random_->Uniform();
        if(1.-val_cos_theta*val_cos_theta>val_random){
          cos_theta_decay = val_cos_theta;
          break;
        }
      }
  }
  Double_t phi_decay = random_->Uniform(-TMath::Pi(),TMath::Pi());
  TVector3 vec_decay_direction = vec_reference_direction_;
  vec_decay_direction.RotateZ(-phi_l);
  vec_decay_direction.RotateY(TMath::ACos(cos_theta_decay));
  vec_decay_direction.RotateZ(phi_l);
  vec_decay_direction.Rotate(phi_decay,vec_reference_direction_);

  Double_t momentum_decay = sqrt( (lv_kpp_.M2()-(kLambdaMass+kProtonMass)*(kLambdaMass+kProtonMass))*(lv_kpp_.M2()-(kLambdaMass-kProtonMass)*(kLambdaMass-kProtonMass)) )/2./lv_kpp_.M();

  TVector3 vec_lambda = momentum_decay*vec_decay_direction;
  TVector3 vec_proton = -vec_lambda;

  lv_lambda_.SetVectM(vec_lambda,kLambdaMass);
  lv_proton_.SetVectM(vec_proton,kProtonMass);

  TVector3 vec_boost_kpp = lv_kpp_.BoostVector();
  lv_lambda_.Boost(vec_boost_kpp);
  lv_proton_.Boost(vec_boost_kpp);

  // spin directions
  if(ml_==-1){
    vec_lambda_spin_direction_ = vec_reference_direction_;
    vec_proton_spin_direction_ = vec_reference_direction_;
  }
  if(ml_==0){
    if(random_->Uniform()<0.5){
      vec_lambda_spin_direction_ = vec_reference_direction_;
    }
    else{
      vec_lambda_spin_direction_ = -vec_reference_direction_;
    }
    vec_proton_spin_direction_ = -vec_lambda_spin_direction_;
  }
  if(ml_==1){
    vec_lambda_spin_direction_ = -vec_reference_direction_;
    vec_proton_spin_direction_ = -vec_reference_direction_;
  }

  return;
}

void KppGenerator::KppZeroPlusDecay()
{
  // reference direction (isotropic)
  Double_t cos_theta_ref = random_->Uniform(-1.,1.);
  Double_t sin_theta_ref = sqrt(1. - cos_theta_ref*cos_theta_ref);
  Double_t phi_ref = random_->Uniform(-TMath::Pi(),TMath::Pi());
  vec_reference_direction_ = TVector3(sin_theta_ref*cos(phi_ref),sin_theta_ref*sin(phi_ref),cos_theta_ref);

  // ml is always 0 (S-wave decay)
  ml_ = 0;

  // momentum directions (isotropic due to S-wave decay)
  Double_t cos_theta_decay = random_->Uniform(-1.,1.);
  Double_t sin_theta_decay = sqrt(1. - cos_theta_decay*cos_theta_decay);
  Double_t phi_decay = random_->Uniform(-TMath::Pi(),TMath::Pi());
  TVector3 vec_decay_direction = TVector3(sin_theta_decay*cos(phi_decay),sin_theta_decay*sin(phi_decay),cos_theta_decay);

  Double_t momentum_decay = sqrt( (lv_kpp_.M2()-(kLambdaMass+kProtonMass)*(kLambdaMass+kProtonMass))*(lv_kpp_.M2()-(kLambdaMass-kProtonMass)*(kLambdaMass-kProtonMass)) )/2./lv_kpp_.M();

  TVector3 vec_lambda = momentum_decay*vec_decay_direction;
  TVector3 vec_proton = -vec_lambda;

  lv_lambda_.SetVectM(vec_lambda,kLambdaMass);
  lv_proton_.SetVectM(vec_proton,kProtonMass);

  TVector3 vec_boost_kpp = lv_kpp_.BoostVector();
  lv_lambda_.Boost(vec_boost_kpp);
  lv_proton_.Boost(vec_boost_kpp);

  // spin directions (always antiparallel)
  if(random_->Uniform()<0.5){
    vec_lambda_spin_direction_ = vec_reference_direction_;
  }
  else{
    vec_lambda_spin_direction_ = -vec_reference_direction_;
  }
  vec_proton_spin_direction_ = -vec_lambda_spin_direction_;

  return;
}

void KppGenerator::KppOnePlusDecay()
{
  // reference direction (isotropic)
  Double_t cos_theta_ref = random_->Uniform(-1.,1.);
  Double_t sin_theta_ref = sqrt(1. - cos_theta_ref*cos_theta_ref);
  Double_t phi_ref = random_->Uniform(-TMath::Pi(),TMath::Pi());
  vec_reference_direction_ = TVector3(sin_theta_ref*cos(phi_ref),sin_theta_ref*sin(phi_ref),cos_theta_ref);

  // ml (of J)
  Double_t random_ml = random_->Uniform(0.,3.);
  if(random_ml<1.){
    ml_ = -1;
  }
  else if(random_ml<2.){
    ml_ = 0;
  }
  else{
    ml_ = 1;
  }

  // momentum directions (isotropic due to S-wave decay)
  Double_t cos_theta_decay = random_->Uniform(-1.,1.);
  Double_t sin_theta_decay = sqrt(1. - cos_theta_decay*cos_theta_decay);
  Double_t phi_decay = random_->Uniform(-TMath::Pi(),TMath::Pi());
  TVector3 vec_decay_direction = TVector3(sin_theta_decay*cos(phi_decay),sin_theta_decay*sin(phi_decay),cos_theta_decay);

  Double_t momentum_decay = sqrt( (lv_kpp_.M2()-(kLambdaMass+kProtonMass)*(kLambdaMass+kProtonMass))*(lv_kpp_.M2()-(kLambdaMass-kProtonMass)*(kLambdaMass-kProtonMass)) )/2./lv_kpp_.M();

  TVector3 vec_lambda = momentum_decay*vec_decay_direction;
  TVector3 vec_proton = -vec_lambda;

  lv_lambda_.SetVectM(vec_lambda,kLambdaMass);
  lv_proton_.SetVectM(vec_proton,kProtonMass);

  TVector3 vec_boost_kpp = lv_kpp_.BoostVector();
  lv_lambda_.Boost(vec_boost_kpp);
  lv_proton_.Boost(vec_boost_kpp);

  // spin directions
  if(ml_==-1){
    vec_lambda_spin_direction_ = -vec_reference_direction_;
    vec_proton_spin_direction_ = -vec_reference_direction_;
  }
  if(ml_==0){
    if(random_->Uniform()<0.5){
      vec_lambda_spin_direction_ = vec_reference_direction_;
    }
    else{
      vec_lambda_spin_direction_ = -vec_reference_direction_;
    }
    vec_proton_spin_direction_ = -vec_lambda_spin_direction_;
  }
  if(ml_==1){
    vec_lambda_spin_direction_ = vec_reference_direction_;
    vec_proton_spin_direction_ = vec_reference_direction_;
  }

  return;
}

void KppGenerator::Print()
{
  std::cout << "KppGenerator::Print()------------------ " << std::endl;
  std::cout << "kpp_spin : " << kpp_spin_ << std::endl;
  std::cout << "kpp_parity : " << kpp_parity_ << std::endl;
  std::cout << "ml : " << ml_ << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "lv_kpp : " << lv_kpp_[0] << ", " << lv_kpp_[1] << ", " << lv_kpp_[2] << ", " << lv_kpp_[3] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "lv_lambda : " << lv_lambda_[0] << ", " << lv_lambda_[1] << ", " << lv_lambda_[2] << ", " << lv_lambda_[3] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "lv_proton : " << lv_proton_[0] << ", " << lv_proton_[1] << ", " << lv_proton_[2] << ", " << lv_proton_[3] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "lv_neutron : " << lv_neutron_[0] << ", " << lv_neutron_[1] << ", " << lv_neutron_[2] << ", " << lv_neutron_[3] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "vec_reference_direction : " << vec_reference_direction_[0] << ", " << vec_reference_direction_[1] << ", " << vec_reference_direction_[2] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "vec_lambda_spin_direction : " << vec_lambda_spin_direction_[0] << ", " << vec_lambda_spin_direction_[1] << ", " << vec_lambda_spin_direction_[2] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "vec_proton_spin_direction : " << vec_proton_spin_direction_[0] << ", " << vec_proton_spin_direction_[1] << ", " << vec_proton_spin_direction_[2] << std::endl;
  std::cout << "--------------------------------------- " << std::endl;
}
