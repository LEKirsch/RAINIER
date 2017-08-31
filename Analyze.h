#include <TStyle.h>
#include "TCanvas.h"
#include "TLegend.h"
#include <vector>
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"

#include "th22mama.C"
/******************************************************************************/
/***************************** RAINIER Analysis *******************************/
/******************************************************************************/
/////////////////////////////// Populations ////////////////////////////////////
TF1 *fnECrit = new TF1("fnECrit", Form("%f", g_dECrit), -g_dPlotSpMax, 
  g_dPlotSpMax);
void AnalyzePop(int exim0 = g_nExIMean - 1, int real0 = g_nReal -1) {
  TCanvas *cPop = new TCanvas("cPop","cPop", 800, 650);
  cPop->SetLogz();
  g_ah2PopLvl[real0][exim0]->Draw("colz");
  fnECrit->Draw("Same");
} // Pop

////////////////////////////// Spin Populations ////////////////////////////////
void AnalyzeJPop(int exim0 = g_nExIMean - 1, int real0 = g_nReal -1) {
  TGraphErrors *agrJPop[g_nReal][g_nExIMean];
  for(int real=0; real<g_nReal; real++) {
    for(int exim=0; exim<g_nExIMean; exim++) {
      g_ahJPop[real][exim]->Scale( 1 / double(g_nEvent) );
      agrJPop[real][exim] = new TGraphErrors(g_nPopLvl);
      for(int lvl=0; lvl<g_nPopLvl; lvl++) {
        int nLvl = g_anPopLvl[lvl];
        int nLvlJ = g_adDisSp[nLvl];
        int nLvlPop = g_ahDisPop[real][exim]->GetBinContent(nLvl);

        agrJPop[real][exim]->SetPoint(lvl, nLvlJ + 0.5,
          nLvlPop / double(g_nEvent) );
        agrJPop[real][exim]->SetPointError(lvl, 0,
          sqrt(nLvlPop) / double(g_nEvent));
      } // lvl
      agrJPop[real][exim]->Write(Form("grExI%dJPop_%d",exim,real));
    } // exim
  } // real

  TCanvas *cJPop = new TCanvas("cJPop","cJPop",800,650);
  // plot one initial J distribution
  g_ahJPop[real0][exim0]->SetLineColor(kBlack);
  g_ahJPop[real0][exim0]->SetLineWidth(2);
  g_ahJPop[real0][exim0]->GetXaxis()->SetTitle("J (#hbar)");
  g_ahJPop[real0][exim0]->GetYaxis()->SetTitle("Population");
  g_ahJPop[real0][exim0]->GetXaxis()->SetNdivisions(10,0,0,kFALSE);
  g_ahJPop[real0][exim0]->Draw();
  TLegend *legJ = new TLegend(0.89,0.51,0.7,0.89);
  //legJ->AddEntry(g_ahJPop[real0][exim0],"J_{I} Dist","L");
  int anColor[] = {kRed, kGreen, kBlue, kMagenta};

  for(int real=0; real<g_nReal; real++) {
    agrJPop[real][exim0]->SetLineColor(anColor[real]);
    agrJPop[real][exim0]->SetLineWidth(2);
    agrJPop[real][exim0]->Draw("same LPE");
    legJ->AddEntry(agrJPop[real][exim0],
      Form("Real %d",real),"LPE");
  }
  legJ->Draw("SAME");
} // J Pop

//////////////////////////// TSC ///////////////////////////////////////////////
void AnalyzeTSC(int nDisEx, int exim0 = g_nExIMean - 1) {
  TCanvas *cTSC = new TCanvas("cTSC","cTSC", 800, 650);
  cTSC->SetLogy();
  for(int real=0; real<g_nReal; real++) {
    g_ahTSC[real][exim0][nDisEx]->SetLineColor(real+2);
    g_ahTSC[real][exim0][nDisEx]->Draw("same");
  } // real
} // TSC

/////////////////////////// Primary 2+ TSC /////////////////////////////////////
TGraphErrors *agrPrim2[g_nReal][g_nExIMean];
void ScalePrim2(int exim, double dScale = 0.75) {
  for(int real=0; real<g_nReal; real++) {
    for(int prim2=0; prim2<g_nPrim2; prim2++) {
      double dEg, dY;
      agrPrim2[real][exim]->GetPoint(prim2, dEg, dY);
      double dEgErr = agrPrim2[real][exim]->GetErrorX(prim2);
      double dYErr  = agrPrim2[real][exim]->GetErrorY(prim2);
      agrPrim2[real][exim]->SetPoint(prim2, dEg, dY * dScale);
      agrPrim2[real][exim]->SetPointError(prim2, dEgErr, dYErr * dScale);
    } // prim2
  } // real
} // scale Prim2

void AnalyzePrim2() { // Primary to select levels
  int nMaxPrim2 = 0;
  for(int real=0; real<g_nReal; real++) {
    for(int exim=0; exim<g_nExIMean; exim++) {
      agrPrim2[real][exim] = new TGraphErrors(g_nPrim2);
      for(int prim2=0; prim2<g_nPrim2; prim2++) {
        int nPrim2Int = g_ahPrim2[real][exim][prim2]->Integral();
        if(nPrim2Int > nMaxPrim2) nMaxPrim2 = nPrim2Int;
        double dEi = g_adExIMean[exim];
        double dEf = g_adDisEne[g_anPrim2[prim2]];
        double dEg = dEi - dEf;
        double dGSF = nPrim2Int / pow(dEg,3);
          // dont have a convienent scale; dont know tot width or
          // cross sec for abs mag; diff Ei have diff vals
          // going to normalize later anyways
        agrPrim2[real][exim]->SetPoint(prim2, dEg, dGSF);
        agrPrim2[real][exim]->SetPointError(prim2, g_dExISpread,
          sqrt(nPrim2Int) / pow(dEg,3));
          // need to put error on dEg spread and density
      } // prim2
      agrPrim2[real][exim]->Write(Form("grExI%dPrim2_%d",exim,real));
    } // exim
  } // real
    
  TCanvas *cGrPrim2 = new TCanvas("cGrPrim2","cGrPrim2",800,650);
  cGrPrim2->SetLogy();
  TF1 *fE1M1E2 = new TF1("fE1M1E2", // assume a const temp
    "GetStrE1(5.0,x)/x**3 + GetStrM1(x)/x**3 + GetStrE2(x)/x**5",0,g_dExIMax);
  TH1D *hfEmpty = new TH1D("hfEmpty","hfEmpty",1000,0,20.0);
  hfEmpty->GetXaxis()->SetTitleSize(0.055);
  hfEmpty->GetXaxis()->SetTitleFont(132);
  hfEmpty->GetXaxis()->SetTitleOffset(0.8);
  hfEmpty->GetXaxis()->CenterTitle();

  hfEmpty->GetYaxis()->SetTitleSize(0.055);
  hfEmpty->GetYaxis()->SetTitleFont(132);     
  hfEmpty->GetYaxis()->SetTitleOffset(0.8);
  hfEmpty->GetYaxis()->CenterTitle();

  hfEmpty->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
  hfEmpty->GetYaxis()->SetTitle("f(E_{#gamma}) (arb)");
  hfEmpty->GetXaxis()->SetRangeUser(0.0,8.0);
  hfEmpty->GetYaxis()->SetRangeUser(1e-9,1e-7);
  hfEmpty->Draw();

  TLegend *legPrim = new TLegend(0.377,0.606,0.642,0.877);
  legPrim->AddEntry(fE1M1E2,"Input f(E_{#gamma})","L");
  for(int real=0; real<g_nReal; real++) {
    for(int exim=0; exim<g_nExIMean; exim++) {
      ///// Normalize /////
      double dMeanX = agrPrim2[real][exim]->GetMean(1); 
      double dMeanY = agrPrim2[real][exim]->GetMean(2);
      double dfEg = fE1M1E2->Eval(dMeanX);
      ScalePrim2(exim, dfEg / dMeanY); // not perfect, have to scale a bit more
      
      agrPrim2[real][exim]->SetLineColor(exim+2);
      agrPrim2[real][exim]->SetMarkerColor(exim+2);
      agrPrim2[real][exim]->SetMarkerStyle(20+exim);
      agrPrim2[real][exim]->SetMarkerSize(1.5);
      agrPrim2[real][exim]->SetLineWidth(2);
      agrPrim2[real][exim]->Draw("same LPE");
      double dExI = g_adExIMean[exim];
      legPrim->AddEntry(agrPrim2[real][exim],Form("E_{x,I} = %2.0f MeV",
        dExI),"P");
    } // exim
  } // real
  fE1M1E2->SetLineColor(kBlack);
  fE1M1E2->Draw("same"); // want on top
  legPrim->Draw();

  // Adjustments
  //ScalePrim2(0,2.0);
  //ScalePrim2(1,1.5);
  //ScalePrim2(2,1.25);
  //ScalePrim2(3,1.125);

} // Prim2

//////////////////////////// Initial Pop Analysis //////////////////////////////
void AnalyzePopDist() {
  TCanvas *cPopDist = new TCanvas("cPopDist", "cPopDist", 580, 650);
  cPopDist->SetLogz();
  g_h2PopDist->GetXaxis()->SetTitleSize(0.055);
  g_h2PopDist->GetXaxis()->SetTitleFont(132);
  g_h2PopDist->GetXaxis()->SetTitleOffset(0.8);
  g_h2PopDist->GetXaxis()->CenterTitle();
  g_h2PopDist->GetXaxis()->SetTitle("J_{I} (#hbar)");

  g_h2PopDist->GetYaxis()->SetTitleSize(0.055);
  g_h2PopDist->GetYaxis()->SetTitleFont(132);
  g_h2PopDist->GetYaxis()->SetTitleOffset(0.65);
  g_h2PopDist->GetYaxis()->CenterTitle();
  g_h2PopDist->GetYaxis()->SetTitle("E_{x,I} (MeV)");

  g_h2PopDist->GetZaxis()->SetTitleSize(0.045);
  g_h2PopDist->GetZaxis()->SetTitleFont(132);
  g_h2PopDist->GetZaxis()->SetTitleOffset(-0.4);
  g_h2PopDist->GetZaxis()->CenterTitle();
  g_h2PopDist->GetZaxis()->SetTitle("Relative Counts");

  g_h2PopDist->Draw("COLZ");
  fnECrit->Draw("SAME");
}

void AnalyzePopI(int exim0 = g_nExIMean - 1, int real0 = g_nReal -1) {
  TCanvas *cPopI = new TCanvas("cPopI","cPopI", 800,650);
  cPopI->SetLogz();
  g_h2PopI->GetXaxis()->SetTitleSize(0.055);
  g_h2PopI->GetXaxis()->SetTitleFont(132);
  g_h2PopI->GetXaxis()->SetTitleOffset(0.8);
  g_h2PopI->GetXaxis()->CenterTitle();

  g_h2PopI->GetYaxis()->SetTitleSize(0.055);
  g_h2PopI->GetYaxis()->SetTitleFont(132);
  g_h2PopI->GetYaxis()->SetTitleOffset(0.7);
  g_h2PopI->GetYaxis()->CenterTitle();

  g_h2PopI->GetZaxis()->SetTitleSize(0.045);
  g_h2PopI->GetZaxis()->SetTitleFont(132);
  g_h2PopI->GetZaxis()->SetTitleOffset(-0.4);
  g_h2PopI->GetZaxis()->CenterTitle();
  g_h2PopI->GetZaxis()->SetTitle("Counts");

  g_h2PopI->GetXaxis()->SetTitle("J_{I} (#hbar)");
  g_h2PopI->GetYaxis()->SetTitle("E_{x,I} (MeV)");

  g_h2PopI->Draw("COLZ");
  fnECrit->Draw("SAME");
} // initial pop

//////////////////////////// Oslo Analysis /////////////////////////////////////
void AnalyzeOslo(int exim0 = g_nExIMean - 1, int real0 = g_nReal -1) {
  TCanvas *cExEg = new TCanvas("cExEg","cExEg", 800, 650);
  cExEg->SetLogz();
  g_ah2ExEg[real0][exim0]->Draw("colz");
  th22mama(g_ah2ExEg[real0][exim0],"ExEg.m");

  TCanvas *c1Gen= new TCanvas("c1Gen","c1Gen", 800, 650);
  c1Gen->SetLogz();
  g_ah21Gen[real0][exim0]->Draw("colz");
  th22mama(g_ah21Gen[real0][exim0],"1Gen.m");
} // Oslo

//////////////////////////// Gamma Analysis ////////////////////////////////////
void AnalyzeGamma(int exim0 = g_nExIMean - 1, int real0 = g_nReal -1) {
  TCanvas *cGSpec= new TCanvas("cGSpec","cGSpec", 800, 650);
  cGSpec->SetLogy();
  g_ahGSpec[real0][exim0]->Draw("colz");
} // gammas

/////////////////////////// Feeding Analysis ///////////////////////////////////
void AnalyzeFeed(int exim0 = g_nExIMean - 1, int real0 = g_nReal -1) {
  TCanvas *cFeed = new TCanvas("cFeed","cFeed", 800, 650);
  cFeed->SetLogz();
  g_ah2FeedTime[real0][exim0]->Draw("colz");
  TGraphErrors *agrLvlFeed[g_nReal][g_nDisLvlMax];
  TH1D *hLvlTimeProj[g_nReal][g_nExIMean][g_nDisLvlMax];

  double dFeedTimeMax = 220;
  TF1 *expo2 = new TF1("expo2", "[0]*exp(-x/[1]) + [2]*exp(-x/[3])",
    0.001, dFeedTimeMax);
  expo2->SetLineColor(kBlack);

  TCanvas *cFit = new TCanvas("cFit","cFit", 100,100);
  for(int real=0; real<g_nReal; real++) {
    for(int lvl=0; lvl<g_nDisLvlMax; lvl++) {
      agrLvlFeed[real][lvl] = new TGraphErrors(g_nExIMean);

      for(int exim=0; exim<g_nExIMean; exim++) {
        double dExIMean = g_adExIMean[exim];

        hLvlTimeProj[real][exim][lvl]
          = (TH1D*)g_ah2FeedTime[real][exim]->ProjectionY(
            Form("hEx%dFeed%d_%d",exim,lvl,real), lvl, lvl+1);

        double dMaxCount = hLvlTimeProj[real][exim][lvl]->GetMaximum();
        expo2->SetParameter(0,dMaxCount/1.3);
        expo2->SetParameter(1,1.0);
        expo2->SetParameter(2,dMaxCount/3.0);
        expo2->SetParameter(3,20.0);
        hLvlTimeProj[real][exim][lvl]->Fit(expo2,"q","goff");
        hLvlTimeProj[real][exim][lvl]->GetYaxis()->SetTitleSize(0.055);
        hLvlTimeProj[real][exim][lvl]->GetYaxis()->SetTitleFont(132);
        hLvlTimeProj[real][exim][lvl]->GetYaxis()->SetTitleOffset(0.8);
        hLvlTimeProj[real][exim][lvl]->GetYaxis()->CenterTitle();
        hLvlTimeProj[real][exim][lvl]->GetYaxis()->SetTitle("Counts");
        double a0  = expo2->GetParameter(0); // dStep1Mag
        double da0 = expo2->GetParError(0);  // dStep1MagErr
        double a1  = expo2->GetParameter(1); // dStep1LifeT
        double da1 = expo2->GetParError(1);  // dStep1LifeTErr

        double a2  = expo2->GetParameter(2); // dStep2Mag
        double da2 = expo2->GetParError(2);  // dStep2MagErr
        double a3  = expo2->GetParameter(3); // dStep2LifeT
        double da3 = expo2->GetParError(3);  // dStep2LifeTErr

        double dFeedTimeMean = (a0*pow(a1,2) + a2*pow(a3,2))/(a0*a1 + a2*a3);
        // Mathematica is the best way to get this:
        // tmean[a0_, a1_, a2_, a3_] := (a0*a1^2 + a2*a3^2)/(a0*a1 + a2*a3); 
        // CForm[ FullSimplify[
        //  D[tmean[a0, a1, a2, a3], a0]^2*da0^2 + 
        //  D[tmean[a0, a1, a2, a3], a1]^2*da1^2 + 
        //  D[tmean[a0, a1, a2, a3], a2]^2*da2^2 + 
        //  D[tmean[a0, a1, a2, a3], a3]^2*da3^2]]
        // then Power -> pow, and tack on a sqrt
        double dFeedTimeMeanErr = sqrt(
          (pow(a1,2)*pow(a2,2)*pow(a1 - a3,2)*pow(a3,2)*pow(da0,2) +
          pow(a0,2)*pow(a0*pow(a1,2) + a2*(2*a1 - a3)*a3,2)*pow(da1,2) +
          pow(a0,2)*pow(a1,2)*pow(a1 - a3,2)*pow(a3,2)*pow(da2,2) +
          pow(a2,2)*pow(a0*a1*(a1 - 2*a3) - a2*pow(a3,2),2)*pow(da3,2))
          /pow(a0*a1 + a2*a3,4));

        // T_1/2 = log(2) * lifetime
        agrLvlFeed[real][lvl]->SetPoint(exim, dExIMean, dFeedTimeMean);
        agrLvlFeed[real][lvl]->SetPointError(
          exim,g_dExISpread,dFeedTimeMeanErr);
      } // initial excitation mean
      agrLvlFeed[real][lvl]->Write(Form("grLvlFeed%d_%d",lvl,real));
    } // lvl 
  } // realization
  cFit->Close();

  TCanvas *cFeedMean = new TCanvas("cFeedMean","cFeedMean", 800,650);
  // TGraphErrors give bad plot control, empty histogram is better in ROOT
  TH1D *hPlotEmpty = new TH1D("hPlotEmpty","Level Feed Times",
    300, 0, g_adExIMean[g_nExIMean-1] * 1.1);
  hPlotEmpty->GetXaxis()->SetTitleSize(0.055);
  hPlotEmpty->GetXaxis()->SetTitleFont(132);
  hPlotEmpty->GetXaxis()->SetTitleOffset(0.8);
  hPlotEmpty->GetXaxis()->CenterTitle();
  hPlotEmpty->GetXaxis()->SetRangeUser(3.0, 10.5);

  hPlotEmpty->GetYaxis()->SetTitleSize(0.055);
  hPlotEmpty->GetYaxis()->SetTitleFont(132);
  hPlotEmpty->GetYaxis()->SetTitleOffset(0.8);
  hPlotEmpty->GetYaxis()->CenterTitle();

  hPlotEmpty->GetXaxis()->SetTitle(
    "Initial Excitation Energy #bar{E}_{x,I} (MeV)");
  hPlotEmpty->GetYaxis()->SetTitle("Avg. Feeding Time #bar{t} (fs)");
  hPlotEmpty->GetYaxis()->SetRangeUser(0.0, 77.0);
  hPlotEmpty->Draw();

  TLegend *legFeed = new TLegend(0.9,0.7,0.7,0.9,"Fed Lvls");
  int anLvlCheck[] = {3,8,10}; // only measure these feed times experimentally
  const int nLvlCheck = sizeof(anLvlCheck) / sizeof(int);
  // plot last real; declared earlier

  for(int chk=0; chk<nLvlCheck; chk++) {
    int lvl = anLvlCheck[chk];
    agrLvlFeed[real0][lvl]->SetLineColor(chk+1);
    agrLvlFeed[real0][lvl]->SetMarkerStyle(21+chk);
    agrLvlFeed[real0][lvl]->SetMarkerSize(1.5);
    agrLvlFeed[real0][lvl]->SetMarkerColor(chk+1);
    agrLvlFeed[real0][lvl]->SetLineWidth(2);
    agrLvlFeed[real0][lvl]->Draw("same LPE");
    if(g_anDisPar[lvl] == 1) {
      legFeed->AddEntry(agrLvlFeed[real0][lvl],
        Form("%2.3f MeV %1.1f+", g_adDisEne[lvl], g_adDisSp[lvl]), "PE");
    } else {
      legFeed->AddEntry(agrLvlFeed[real0][lvl],
        Form("%2.3f MeV %1.1f-", g_adDisEne[lvl], g_adDisSp[lvl]), "PE");
    } // parity
  } // lvl
  legFeed->Draw("Same");

  TCanvas *cExFeed = new TCanvas("cExFeed","cExFeed",800,650);
  cExFeed->SetLogy();
  TLegend *legExFeed = new TLegend(0.9,0.52,0.7,0.9,"#bar{E}_{x,I} (MeV)");
  int nLvl = 3;
  for(int exim=0; exim<g_nExIMean; exim++) {
    hLvlTimeProj[real0][exim][nLvl]->GetXaxis()->SetTitleOffset(0.8);
    hLvlTimeProj[real0][exim][nLvl]->SetLineColor(exim+2);
    hLvlTimeProj[real0][exim][nLvl]->SetLineWidth(2);
    legExFeed->AddEntry(hLvlTimeProj[real0][exim][nLvl],
      Form("%2.1f MeV", g_adExIMean[exim] ), "L");
    hLvlTimeProj[real0][exim][nLvl]->Draw("same");
  }
  legExFeed->Draw("same");
} // Feed

void AnalyzeTotWid(int exim = g_nExIMean - 1) {
  // RAINIER
  TCanvas *cTotWid = new TCanvas("cTotWid","cTotWid",800,650);

  TH1D *htEmpty = new TH1D("htEmpty","htEmpty",g_nReal,0,g_nReal);
  htEmpty->GetXaxis()->SetTitleSize(0.055);
  htEmpty->GetXaxis()->SetTitleFont(132);
  htEmpty->GetXaxis()->SetTitleOffset(0.8);
  htEmpty->GetXaxis()->CenterTitle();

  htEmpty->GetYaxis()->SetTitleSize(0.055);
  htEmpty->GetYaxis()->SetTitleFont(132);
  htEmpty->GetYaxis()->SetTitleOffset(1.2);
  htEmpty->GetYaxis()->CenterTitle();

  htEmpty->GetXaxis()->SetTitle("Realization Number");
  htEmpty->GetYaxis()->SetTitle("#bar{#Gamma}_{Tot} (MeV)");
  htEmpty->GetXaxis()->SetRangeUser(0,g_nReal);
  htEmpty->GetYaxis()->SetRangeUser(10.15e-9,10.35e-9);
  htEmpty->Draw();

  g_grTotWidAvg[exim]->SetMarkerStyle(22);
  g_grTotWidAvg[exim]->SetLineWidth(2);
  g_grTotWidAvg[exim]->Draw("LP");

  // DICEBOX Width
  ifstream totFile;
  totFile.open("Bench/Re186/Re186/TOTWID.DAT");

  string sHeadLine;
  getline(totFile,sHeadLine);

  vector<double> vdTotWid;
  vector<double> vdTotWidErr;

  int nReal;
  double dTotWid, dTotWidErr;
  string sTotLine;
  while(getline(totFile,sTotLine)) {
    istringstream issTot(sTotLine);
    issTot >> nReal >> dTotWid >> dTotWidErr;
    vdTotWid.push_back(dTotWid);
    vdTotWidErr.push_back(dTotWidErr);
  }

  int nRealTot = vdTotWid.size();
  TGraph *grDICEBOXTotWid = new TGraph(nRealTot);
  for(int real=0; real<nRealTot; real++) {
    grDICEBOXTotWid->SetPoint(real,real,vdTotWid[real]);
  }
  grDICEBOXTotWid->SetLineColor(kRed);
  grDICEBOXTotWid->SetLineWidth(2);
  grDICEBOXTotWid->SetMarkerStyle(21);
  grDICEBOXTotWid->SetMarkerColor(kRed);
  grDICEBOXTotWid->Draw("same LP");
  TLegend *legTotWid = new TLegend(0.217,0.506,0.497,0.657);
  legTotWid->AddEntry(grDICEBOXTotWid,"DICEBOX","LP");
  legTotWid->AddEntry(g_grTotWidAvg[exim],"RAINIER","LP");
  legTotWid->Draw("SAME");

  double dRAINIERAvgTotWid,dDICEBOXAvgTotWid, dR;
  g_grTotWidAvg[exim]->GetPoint(g_nReal-1,dR,dRAINIERAvgTotWid);
  grDICEBOXTotWid->GetPoint(nRealTot-1,dR,dDICEBOXAvgTotWid);

  cout << "Average total width after " << g_nReal << " " << nRealTot
    << " realizations:" << endl;
  cout << "DICEBOX: " << dDICEBOXAvgTotWid << " MeV" << endl;
  cout << "RAINIER: " << dRAINIERAvgTotWid << " MeV" << endl;
  cout << "difference: " << (dDICEBOXAvgTotWid - dRAINIERAvgTotWid) 
    / dDICEBOXAvgTotWid * 100.0 << " \%" << endl;
} // benchmark width

void Bench() {
  TCanvas *cTotWid = new TCanvas("cTotWid","cTotWid",800,650);

  TH1D *hpEmpty = new TH1D("hpEmpty","hpEmpty",g_nReal,0,g_nReal);
  hpEmpty->GetXaxis()->SetTitleSize(0.055);
  hpEmpty->GetXaxis()->SetTitleFont(132);
  hpEmpty->GetXaxis()->SetTitleOffset(0.8);
  hpEmpty->GetXaxis()->CenterTitle();

  hpEmpty->GetYaxis()->SetTitleSize(0.055);
  hpEmpty->GetYaxis()->SetTitleFont(132);
  hpEmpty->GetYaxis()->SetTitleOffset(0.8);
  hpEmpty->GetYaxis()->CenterTitle();

  hpEmpty->GetXaxis()->SetTitle("Realization Number");
  hpEmpty->GetYaxis()->SetTitle("Average Level Population (\%)");
  hpEmpty->GetXaxis()->SetRangeUser(0,g_nReal-1);
  hpEmpty->GetYaxis()->SetRangeUser(0.0,100);
  hpEmpty->Draw();

  int nPlotLvlStart = 1;
  const int nPlotLvl = 6;
  // RAINIER
  TLegend *legPopRAINIER = new TLegend(0.7,0.88,0.8,0.6,"RAINIER");
  TGraph *grPopAvgRAINIER[nPlotLvl];
  for(int lvl=nPlotLvlStart; lvl<nPlotLvl; lvl++) {
    grPopAvgRAINIER[lvl] = new TGraph(g_nReal);
    for(int real=0; real<g_nReal; real++) {
      if(real != 0) {
        double dOldPopAvg, dR;
        grPopAvgRAINIER[lvl]->GetPoint(real-1, dR, dOldPopAvg);
        double dNewPopAvg = (dOldPopAvg * real 
          + g_ahDisPop[real][0]->GetBinContent(lvl+1) / double(g_nEvent) 
          * 100.0 ) / (real + 1);
        grPopAvgRAINIER[lvl]->SetPoint(real,real,dNewPopAvg);
      } else { // 1st real
        grPopAvgRAINIER[lvl]->SetPoint(real,real,
          g_ahDisPop[real][0]->GetBinContent(lvl+1) / double(g_nEvent) 
          * 100.0 );
      } // 1st real?
    } // real
    grPopAvgRAINIER[lvl]->SetLineColor(lvl-nPlotLvlStart+1);
    grPopAvgRAINIER[lvl]->SetLineWidth(2);
    grPopAvgRAINIER[lvl]->SetMarkerStyle(20+lvl);
    grPopAvgRAINIER[lvl]->SetMarkerColor(lvl-nPlotLvlStart+1);
    grPopAvgRAINIER[lvl]->Draw("same LP");
    legPopRAINIER->AddEntry(grPopAvgRAINIER[lvl], Form("Lvl %d",lvl), "LP");
  } // lvl
  legPopRAINIER->Draw("SAME");

  // DICEBOX
  TGraph *grPopAvgDICEBOX[nPlotLvl];
  TLegend *legPopDICEBOX= new TLegend(0.55,0.88,0.65,0.6,"DB");
  for(int lvl=nPlotLvlStart; lvl<nPlotLvl; lvl++) {
    ifstream popFile;
    popFile.open(Form("Bench/Re186/Re186/Pop%d.dat",lvl+1));
    //popFile.open(Form("Pop%d.dat",lvl+1));

    string sHeadLine;
    getline(popFile,sHeadLine);

    vector<double> vdPop;
    vector<double> vdPopErr;

    int nReal,nLvl;
    double dEnergy, dPop, dPopErr;
    string sPopLine;
    while(getline(popFile,sPopLine)) {
      istringstream issPop(sPopLine);
      issPop >> nReal >> nLvl >> dEnergy >> dPop >> dPopErr;
      vdPop.push_back(dPop);
      vdPopErr.push_back(dPopErr);
    } // getline

    int nRealPop = vdPop.size();
    grPopAvgDICEBOX[lvl] = new TGraph(nRealPop);
    for(int real=0; real<nRealPop; real++) {
      grPopAvgDICEBOX[lvl]->SetPoint(real,real,vdPop[real]*100.0);
      //grPopAvgDICEBOX[lvl]->SetPointError(real,0,vdPopErr[real]);
    } // real
    grPopAvgDICEBOX[lvl]->SetLineColor(nPlotLvl+lvl-2*nPlotLvlStart+1);
    grPopAvgDICEBOX[lvl]->SetLineWidth(2);
    grPopAvgDICEBOX[lvl]->SetMarkerStyle(20+nPlotLvl+lvl);
    grPopAvgDICEBOX[lvl]->SetMarkerColor(nPlotLvl+lvl-2*nPlotLvlStart+1);
    grPopAvgDICEBOX[lvl]->Draw("same LP");
    legPopDICEBOX->AddEntry(grPopAvgDICEBOX[lvl], Form("Lvl %d",lvl), "LP");

  } // lvl
  legPopDICEBOX->Draw("SAME");
  
} // AnalyzePopBench

TF1 *fnLDa    = new TF1("fnLDa",    "GetLDa(x)",0,10);
TF1 *fnSpCut  = new TF1("fnSpCut",  "sqrt(GetSpinCut2(x))",0,10);
TF1 *fnGSFE1  = new TF1("fnGSFE1",  "GetStrE1([0],x)/x**3",0,18);
TF1 *fnGSFM1  = new TF1("fnGSFM1",  "GetStrM1(x)/x**3",0,12);
TF1 *fnGSFE2  = new TF1("fnGSFE2",  "GetStrE2(x)/x**5",0,12);
TF1 *fnGSFTot = new TF1("fnGSFTot", "GetStrE1([0],x)/x**3 + GetStrM1(x)/x**3 + GetStrE2(x)/x**5",0,12);
