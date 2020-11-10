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
  Double_t phi_spin = vec_lambda_spin_direction_.Phi();

  // momentum directions
  Double_t cos_theta_decay = -999.;
  while(true){
    Double_t val_cos_theta = random_->Uniform(-1.,1.);
    Double_t val_random = random_->Uniform(0.,1.+alpha_);
    if(1.+alpha_*val_cos_theta>val_random){
      cos_theta_decay = val_cos_theta;
      break;
    }
  }
  Double_t phi_decay = random_->Uniform(-TMath::Pi(),TMath::Pi());
  TVector3 vec_decay_direction = vec_lambda_spin_direction_;
  vec_decay_direction.RotateZ(-phi_spin);
  vec_decay_direction.RotateY(TMath::ACos(cos_theta_decay));
  vec_decay_direction.RotateZ(phi_spin);
  vec_decay_direction.Rotate(phi_decay,vec_lambda_spin_direction_);
  
  Double_t momentum_decay = sqrt( (lv_lambda_.M2()-(kProtonMass+kPiMass)*(kProtonMass+kPiMass))*(lv_lambda_.M2()-(kProtonMass-kPiMass)*(kProtonMass-kPiMass)) )/2./lv_lambda_.M();
  TVector3 vec_proton = momentum_decay*vec_decay_direction;
  TVector3 vec_pim = -vec_proton;

  lv_proton_ = TLorentzVector(vec_proton,kProtonMass);
  lv_pim_ = TLorentzVector(vec_pim,kPiMass);

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
