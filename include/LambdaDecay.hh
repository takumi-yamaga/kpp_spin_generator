// LambdaDecay.hh

#ifndef LambdaDecay_hh
#define LambdaDecay_hh

#include <iostream>
#include <iomanip>

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "KppGeneratorCommon.hh"

class LambdaDecay
{
  public:
    LambdaDecay();
    ~LambdaDecay();

  public:
    void Generate(TLorentzVector& lv_lambda, TVector3& vec_lambda_spin_direction);
    void Print();

  private:
    TRandom3* random_;
    TLorentzVector lv_lambda_;
    TLorentzVector lv_proton_;
    TLorentzVector lv_pim_;

    TVector3 vec_lambda_spin_direction_;
    Double_t alpha_;
    Double_t beta_;
    Double_t gamma_;

};

#endif
