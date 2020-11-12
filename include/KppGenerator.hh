// KppGenerator.hh

#ifndef KppGenerator_hh
#define KppGenerator_hh

#include <iostream>
#include <iomanip>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TF1.h"

#include "KppGeneratorCommon.hh"

class KppGenerator
{
  public:
    KppGenerator();
    KppGenerator(Int_t kpp_spin, Int_t kpp_parity);
    ~KppGenerator();

  public:
    void Generate(TLorentzVector& lv_beam);
    void Print();

    TLorentzVector LvKpp(){ return lv_kpp_; }
    TLorentzVector LvLambda(){ return lv_lambda_; }
    TLorentzVector LvProton(){ return lv_proton_; }
    TLorentzVector LvNeutron(){ return lv_neutron_; }

    Int_t KppSpin(){ return kpp_spin_; }
    Int_t KppParity(){ return kpp_parity_; }
    Int_t ML(){ return ml_; }
    TVector3 VecReference(){ return vec_reference_direction_; }
    TVector3 VecLambdaSpin(){ return vec_lambda_spin_direction_; }
    TVector3 VecProtonSpin(){ return vec_proton_spin_direction_; }

  private:
    void KppZeroMinusDecay();
    void KppZeroPlusDecay();
    void KppOnePlusDecay();
    TRandom3* random_;

    TF1* f_kpp_mass_;
    TF1* f_kpp_momtrans_;
    
    TLorentzVector lv_target_;
    TGenPhaseSpace generator_;
    Double_t cm_masses_[3] = { kLambdaMass, kProtonMass, kNeutronMass};

    TLorentzVector lv_kpp_;
    TLorentzVector lv_lambda_;
    TLorentzVector lv_proton_;
    TLorentzVector lv_neutron_;

    Int_t kpp_spin_;
    Int_t kpp_parity_;
    Int_t ml_;
    TVector3 vec_reference_direction_;
    TVector3 vec_lambda_spin_direction_;
    TVector3 vec_proton_spin_direction_;
};

#endif
