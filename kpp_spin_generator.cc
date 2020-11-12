#include "KppGenerator.hh"
#include "LambdaDecay.hh"
#include "Analyzer.hh"
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
  Analyzer* analyzer = new Analyzer("root/kpp_generator.root");

  Double_t beam_momentum = 1.; // GeV/c
  TLorentzVector lv_beam(0.,0.,0.,0.);
  lv_beam.SetVectM(TVector3(0.,0.,beam_momentum),kKaonMass);

  Int_t total_events = 100000;
  for(Int_t i_event = 0; i_event<total_events; ++i_event){

    kpp_generator->Generate(lv_beam);
    TLorentzVector lv_lambda = kpp_generator->LvLambda();
    TVector3 vec_lambda_spin_direction = kpp_generator->VecLambdaSpin();
    lambda_decay->Generate(lv_lambda,vec_lambda_spin_direction);

    analyzer->Analysis(kpp_generator,lambda_decay);

    Int_t print_progress = (Int_t)log10(i_event+1);
    print_progress = pow(10,print_progress);
    if((i_event+1)%print_progress==0){
      std::cout << i_event+1 << std::endl;
    }

  }

  analyzer->PrintHistograms("fig/kpp_genrator.pdf");

  return 0;
}
