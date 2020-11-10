#include "KppGenerator.hh"
#include "LambdaDecay.hh"
#include "KppGeneratorCommon.hh"
#include <iostream>
#include <fstream>

// -----=====-----=====----- //
// MAIN FUNCTION             //
// -----=====-----=====----- //
int main(int argc, char* argv[])
{	
  Int_t kpp_spin = 0; 
  Int_t kpp_parity = -1; 
  KppGenerator* kpp_generator = new KppGenerator(kpp_spin,kpp_parity);
  LambdaDecay* lambda_decay = new LambdaDecay();

  Double_t beam_momentum = 1.; // GeV/c
  TLorentzVector lv_beam(0.,0.,beam_momentum,kKaonMass);

  Int_t total_events = 1;
  for(Int_t i_event = 0; i_event<total_events; ++i_event){
    kpp_generator->Generate(lv_beam);
    kpp_generator->Print();
    TLorentzVector lv_lambda = kpp_generator->LvLambda();
    TVector3 vec_lambda_spin_direction = kpp_generator->VecLambdaSpin();
    lambda_decay->Generate(lv_lambda,vec_lambda_spin_direction);
    lambda_decay->Print();
  }

  return 0;
}
