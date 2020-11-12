// LambdaDecay.cc

#include "LambdaDecay.hh"

LambdaDecay::LambdaDecay()
{
  alpha_ = 0.732;
  beta_ = sqrt(1.-alpha_*alpha_)*sin(6.5/180.*TMath::Pi());
  gamma_ = 0.76;

  random_ = new TRandom3(0);
}

LambdaDecay::~LambdaDecay()
{
  delete random_;
}

void LambdaDecay::Generate(TLorentzVector& lv_lambda,TVector3& vec_lambda_spin_direction)
{
  lv_lambda_ = lv_lambda;
  vec_lambda_spin_direction_ = vec_lambda_spin_direction;
  TVector3 vec_decay_direction(0.,0.,0.);

  // momentum directions
  while(true){
    Double_t cos_theta_decay = random_->Uniform(-1.,1.);
    Double_t sin_theta_decay = sqrt(1. - cos_theta_decay*cos_theta_decay);
    Double_t phi_decay = random_->Uniform(-TMath::Pi(),TMath::Pi());
    vec_decay_direction = TVector3(sin_theta_decay*cos(phi_decay),sin_theta_decay*sin(phi_decay),cos_theta_decay);

    Double_t cos_theta_spin = random_->Uniform(-1.,1.);
    Double_t sin_theta_spin = sqrt(1. - cos_theta_spin*cos_theta_spin);
    Double_t phi_spin = random_->Uniform(-TMath::Pi(),TMath::Pi());
    vec_proton_spin_direction_ = TVector3(sin_theta_spin*cos(phi_spin),sin_theta_spin*sin(phi_spin),cos_theta_spin);

    Double_t val_ratio = 1. + gamma_*vec_proton_spin_direction_.Dot(vec_lambda_spin_direction_)
      + (1.-gamma_)*vec_proton_spin_direction_.Dot(vec_decay_direction)*vec_lambda_spin_direction_.Dot(vec_decay_direction)
      + alpha_*(vec_proton_spin_direction_.Dot(vec_decay_direction)+vec_lambda_spin_direction_.Dot(vec_decay_direction))
      + beta_*vec_decay_direction.Dot(vec_proton_spin_direction_.Cross(vec_lambda_spin_direction_));
    Double_t val_random = random_->Uniform(0.,kMaxRatioDecay);
    if(val_ratio>val_random){
      break;
    }
  }

  Double_t momentum_decay = sqrt( (lv_lambda_.M2()-(kProtonMass+kPiMass)*(kProtonMass+kPiMass))*(lv_lambda_.M2()-(kProtonMass-kPiMass)*(kProtonMass-kPiMass)) )/2./lv_lambda_.M();
  TVector3 vec_proton = momentum_decay*vec_decay_direction;
  TVector3 vec_pim = -vec_proton;

  lv_proton_.SetVectM(vec_proton,kProtonMass);
  lv_pim_.SetVectM(vec_pim,kPiMass);

  TVector3 vec_boost_lambda = lv_lambda_.BoostVector();
  lv_proton_.Boost(vec_boost_lambda);
  lv_pim_.Boost(vec_boost_lambda);

  return;
}

void LambdaDecay::Print()
{
  std::cout << "LambdaDecay::Print()------------------ " << std::endl;
  std::cout << "alpha : " << alpha_ << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "lv_lambda : " << lv_lambda_[0] << ", " << lv_lambda_[1] << ", " << lv_lambda_[2] << ", " << lv_lambda_[3] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "lv_proton : " << lv_proton_[0] << ", " << lv_proton_[1] << ", " << lv_proton_[2] << ", " << lv_proton_[3] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "lv_pim : " << lv_pim_[0] << ", " << lv_pim_[1] << ", " << lv_pim_[2] << ", " << lv_pim_[3] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "vec_lambda_spin_direction : " << vec_lambda_spin_direction_[0] << ", " << vec_lambda_spin_direction_[1] << ", " << vec_lambda_spin_direction_[2] << std::endl;
  std::cout << "--------------------------------------- " << std::endl;
}
