// Analyzer.hh

#ifndef Analyzer_hh
#define Analyzer_hh

#include <string>

#include "KppGenerator.hh"
#include "LambdaDecay.hh"

#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLatex.h"

class Analyzer
{
  public:
    Analyzer();
    Analyzer(std::string file_name);
    ~Analyzer();

  public:
    void Analysis(KppGenerator* kpp_generator, LambdaDecay* lambda_decay);
    void PrintHistograms(std::string pdf_name);

  private:
    void CreateHistograms();
    void Fill(std::string hist_name, Double_t val);
    void Fill(std::string hist_name, Double_t val1,Double_t val2);
    void Draw(std::string hist_name, std::string option);

  private:
    std::string root_file_name_;
    TFile* root_file_;

};

#endif
