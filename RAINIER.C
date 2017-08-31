/********************* RAINIER l.kirsch  4/27/2017 start date *****************/
/********************* version 1.0.5:    8/31/2017 ****************************/
/********************* lekirsch@lbl.gov  **************************************/
//  ____________________________________________ 
// |* * * * * * * *|############################|
// | * * * * * * * |                            |
// |* * * * * * * *|############################|
// | * * * * * * * |                            |
// |* * * * * * * *|############################|
// | * * * * * * * |                            |
// |* * * * * * * *|############################|
// |~~~~~~~~~~~~~~~'                            |
// |############################################|
// |                                            |
// |############################################|
// |                                            |
// |############################################|
// '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
// Randomizer of Assorted Initial Nuclear Intensities and Emissions of Radiation
// to run in bash:
// $ root.exe .x RAINIER.C++

// doubles, ints, and bools marked with prescript "d", "n", and "b" respectively
// arrays marked with prescript "a"
// globals marked with "g_" are accessible after the run, so are the functions
// full lowercase named variables (i.e. no CamelCase) are index variables 
// precompiler commands "#" make things run fast when unnecessary items are out
// - better than 40 "if" statements for every call
// - also makes src code clear as to what contributes and what doesn't
// - dormant code is error prone

////////////////////////////////////////////////////////////////////////////////
////////////////////// Input Parameters ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////// Discrete Settings ///////////////////////////////////////
//#define bPrintLvl // print both discrete and constructed lvl schemes
//#define bLvlDICEBOX // if copy last few lines of DICEBOX input file
#define bLvlTALYS // using /home/lekirsch/talys/structure/levels/final/z???
const char zFile[] = "/home/lekirsch/talys/structure/levels/final/z060";
//const char zFile[] = "z060DICE";
const int g_nDisLvlMax = 13; // only trust level scheme to here
const int g_nDisLvlGamMax = 7; // dont have more than this many transistions

///////////////////// Constructed Level Scheme Settings ////////////////////////
///// Bins /////
#define bForceBinNum // else force bin spacing
const int g_nConSpbMax = 21; // constructed # spin bins, small for light ion rxn

#ifdef bForceBinNum
double g_dConESpac; // constructed E bin spacing
const int g_nConEBin = 400; // number of energy bins in constructed scheme
#else
const double g_dConESpac = 0.01; // MeV; wont matter if forcing bin number
int g_nConEBin;
#endif // bForceBinNum

///// Nucleus /////
const int g_nZ = 60; // proton number
const double g_dAMass = 144.0; // proton + neutron number
const bool g_bIsEvenA = true; // set even/odd A

///// Level Density, LD, model /////
// choose one, fill in corresponding parameters
#define bLD_BSFG
//#define bLD_CTM
#ifdef bLD_CTM
const double g_dTemp =  0.48473; // MeV
const double g_dE0   = -1.31817; // MeV
const double g_dDeuPair = 0.62834; // MeV; 
//const double g_dE0   = -1.004 + 0.5 * g_dDeuPair; // MeV; von egidy09 fit
#endif
#ifdef bLD_BSFG
const double g_dE1 = 0.968; // MeV, excitation energy shift
const double g_dDeuPair = 2.698; // MeV; can get from ROBIN: Pa_prime
//const double g_dE1 = g_dDeuPair * 0.5 - 0.381; // von Egidy fit
#endif
// deuteron pairing energy from mass table, related to backshift in BSFG or CTM
// used for effective energy in LD, spincut, GSF models

//#define bLDaEx // a(Ex) = aAsym * (1 + dW * (1 - exp(-Gam * Eff) / dEff) )
#ifdef bLDaEx
const double g_dLDaAsym   = 14.58; // MeV^-1; Asymptotic value, a(Ex->Inf)
const double g_dDampGam   = 0.0; // Damping Parameter
const double g_dShellDelW = 0.0; // MeV; M_exp - M_LDM ~ shell correction
#else // const a
const double g_dLDa = 14.58; // MeV^-1 aka "LD parameter a" 
#endif

///// Spin Cutoff /////
// choose one:
#define bJCut_VonEgidy05 // low-energy model
//#define bJCut_SingPart // single particle model
//#define bJCut_RigidSph // rigid sphere model
//#define bJCut_VonEgidy09 // empirical fit 
//#define bJCut_Other // so many models of spin cut... code your own

///// Pairity Dependence /////
// choose one:
#define bPar_Equipar // equal +/- states, usually a good approx
//#define bPar_Edep // exponential asymptotic dependence
#ifdef bPar_Edep // 0.5 (1 +/- 1 / (1 + exp( g_dParC * (dEx - g_dParD) ) ) )
const double g_dParC = 3.0; // MeV^-1
const double g_dParD = 0.0; // MeV
#endif

///// Level spacing distribution /////
#define bPoisson // good approx to lvl spacing, might have more severe fluct
//#define bWigner // more representative of nuc lvl spacing, but more t to init

/////////////////////// Gamma Strength Function, GSF, Settings /////////////////
///// Width Fluctuations Distribution (WFD)-reason statistical codes exist /////
// choose one:
#define bWFD_PTD // fastest version of the Porter Thomas Distribution, nu = 1
//#define bWFD_nu // set the chi^2 degrees of freedom, nu - slow
//#define bWFD_Off // no fluctuations - nearly TALYS with level spac fluct
#ifdef bWFD_nu
const double g_dNu = 0.5; // See Koehler PRL105,072502(2010): measured nu~0.5
#endif

// Parameters for 56Fe from TALYS defaults usually an acceptable start
///// fE1 /////
// choose one:
#define bE1_GenLor // General Lorentzian
//#define bE1_KMF // Kadmenskij Markushev Furman model
//#define bE1_KopChr // Kopecky Chrien model
//#define bE1_StdLor // standard Lorentzian
const double g_adSigE1[] = {317.00, 0.00}; // mb magnitude
const double g_adEneE1[] = { 15.05, 0.01}; // MeV centroid energy, non-zero
const double g_adGamE1[] = {  5.30, 0.00}; // MeV GDR width

///// fM1 /////
const double g_adSigM1[] = {0.37, 0.00}; // mb magnitude
const double g_adEneM1[] = {7.82, 0.01}; // MeV centroid energy, non-zero
const double g_adGamM1[] = {4.00, 0.00}; // MeV GDR width
//#define bM1StrUpbend // Oslo observed low energy upbend aka enhancement
#ifdef bM1StrUpbend // soft pole behavior: C * exp(-A * Eg)
const double g_dUpbendM1Const = 5e-8; // C
const double g_dUpbendM1Exp = 1.0; // (positive) A
#endif

///// fE2 /////
#define bE2_SingPart // single particle
//#define bE2_StdLor // standard Lorentzian, only one resonance
#ifdef bE2_SingPart
const double g_dSpSigE2 = 4e-11; // MeV^-5
#endif
#ifdef bE2_StdLor
const double g_dSigE2 =  0.39; // mb magnitude
const double g_dEneE2 = 11.04; // MeV centroid energy
const double g_dGamE2 =  3.88; // MeV GDR width
#endif

////////////////////// Run Settings ////////////////////////////////////////////
const int g_nReal = 1; // number of realizations of nuclear level scheme
const int g_nEvent = 1e5; // number of events per excitation and realization
const int g_nEvUpdate = 1e2; // print progress to screen at this interval

////////////////////// Internal Conversion Coefficient, ICC, Settings //////////
//#define bUseICC // ICC = 0.0 otherwise
// if issues, turn ICC off and get your briccs to work outside RAINIER first
const char g_sBrIccModel[] = "BrIccFO"; // Conversion data table
const int g_nBinICC = 100; // Energy bins of BrIcc - more takes lot init time
const double g_dICCMin = g_dConESpac / 2.0; // uses 1st Ebin ICC val below this
const double g_dICCMax = 1.0; // MeV; Uses last Ebin ICC value for higher E

////////////////////// Excitation Settings /////////////////////////////////////
// choose one, fill in corresponding params:
#define bExSingle  // single population input
//#define bExSelect // like Beta decay
//#define bExSpread  // ejectile detected input
//#define bExFullRxn // no ejectile detected input

#ifdef bExSingle
const double g_dExIMax = 7.8174; // MeV, Ei - "capture state energy"
const double g_dSpI    = 3.0; // hbar, Ji - "capture state spin"
const double g_dParI   = 0; // Pi - "capture state parity" 0=(-), 1=(+)

const double g_adExIMean[] = {0.0}; // unused in single
const double g_dExISpread = 0.0; // unused in single
const double g_dExRes = 0.0; // unused in single
#endif

#ifdef bExSelect
const double g_adExI[]  = {0, 0.88489, 1.78662, 2.69324, 2.89503, 3.0383, 3.3418, 3.47668, 3.59912, 3.78828, 3.8483, 3.9039, 4.0016, 4.06142, 4.2646, 4.30903, 4.46481, 4.5144, 4.5583, 4.5889, 4.7101, 4.7917, 4.8492, 5.0613}; // MeV
const double g_adSpI[]  = {0,  2,  4,  4,  6,  5,  3,  5,  6,  6,  5,  5,  6,  5,  6,  5,  7,  6,  5,  7,  7,  5,  5,  5}; // hbar
const double g_anParI[] = {1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0}; // Pi 
const double g_adBRI[]   = {0, 0, 0, 0, 0.014, 0.58, 0.0094, 0.061, 0.0103, 0.0325, 0.049, 0.0057, 0.0074, 0.0141, 0.0216, 0.0236, 0.0415, 0.0135, 0.0298, 0.014, 0.0104, 0.0057, 0.0519, 0.005}; // Branching Ratio: sums to 1.0

const int    g_nStateI = sizeof(g_adExI) / sizeof(double);
const double g_dExIMax = g_adExI[g_nStateI-1] + 0.25; // build above last val

const double g_adExIMean[] = {0.0}; // unused in select
const double g_dExISpread = 0.0; // unused in select
const double g_dExRes = 0.0; // unused in select
#endif

#ifdef bExSpread
const double g_dExIMax = 8.0; // MeV; constructed lvl scheme built up to this
// dont exceed with init excitations - gaus might sample higher than expected
const double g_adExIMean[] = {3.0, 4.0, 5.0, 6.0, 7.0}; // MeV
const double g_dSpIMean = 4.0; // poisson mean for initial spin state spread
const double g_dExISpread = 0.2 / 2.355; // MeV; std dev sigma = FWHM / 2.355 
const double g_dExRes = 0.2; // excitation resolution on g_ah2ExEg for analysis
#endif

#ifdef bExFullRxn // from a TALYS output file if available
const double g_dExIMax = 16.0; // MeV; 16 MeV beam cant excite beyond ~15.7
const char popFile[] = "Fe56Pop.dat"; // made from TALYS "outpopulation y"
// make sure to match # of discrete bins
const double g_dExRes = 0.2 / 2.355; // excitation resolution on g_ah2ExEg

const double g_adExIMean[] = {0.0}; // unsed in full rxn
const double g_dExISpread = 0.0; // unsed in full rxn
#endif

////////////////////// Analysis Settings ///////////////////////////////////////
const double g_dPlotSpMax = 10.0;
///// JPop Analysis /////
const int g_anPopLvl[] = {0}; // low-ly populated lvls, 0 = gs
const int g_nPopLvl = sizeof(g_anPopLvl) / sizeof(int);
///// Prim2 Analysis
const int g_anPrim2[] = {0}; // low-ly primary lvls, 0 = gs
const int g_nPrim2 = sizeof(g_anPrim2)/sizeof(int);
const int g_nEgBinTSC = 500;

/////////////////////////////// Parallel Settings //////////////////////////////
// should handle itself, email me if you get it to work on Mac or PC
#ifdef __CLING__
// cling in root6 won't parse omp.h. but man, root5 flies with 24 cores!
#else
#ifdef __linux__ // MacOS wont run omp.h by default, might exist workaround
    //#define bParallel // Parallel Option
    // ROOT hisograms not thread safe, but only miss ~1e-5 events
#endif // linux
#endif // cling

////////////////////////////////////////////////////////////////////////////////
////////////////////// End Input Parameters ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/////////////////////////////// Program Includes ///////////////////////////////
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <iomanip> 
#include "math.h"
#include "TTimeStamp.h"
using namespace std;
#include "TRandom2.h" 
// 2=Tausworthe is faster and smaller than 3=Mersenne Twister (MT19937)
// search and replace "TRandom2" to "TRandom3" to change PRNG
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include <TROOT.h>
// determine OS for briccs
#ifdef __linux__
char cbriccs[] = "briccs";
#endif
#ifdef __MACH__
char cbriccs[] = "briccsMac";
#endif
#ifdef _WIN32
char cbriccs[] = "BrIccS.exe"; // haven't done any windows testing yet
#endif

#ifdef bParallel
#ifndef CINT 
#include "omp.h" // for parallel on shared memory machine (not cluster yet)
#endif // cint
#endif // parallel

const int g_nExIMean = sizeof(g_adExIMean) / sizeof(double); // self adjusting
const int g_nParE1 = sizeof(g_adSigE1) / sizeof(double);
const int g_nParM1 = sizeof(g_adSigM1) / sizeof(double);

///////////////////////// Discrete Input File //////////////////////////////////
double g_adDisEne[g_nDisLvlMax]; // discrete lvl energy
double g_adDisSp [g_nDisLvlMax]; // discrete lvl spin
int    g_anDisPar[g_nDisLvlMax]; // discrete lvl parity
double g_adDisT12[g_nDisLvlMax]; // discrete lvl half-life
int    g_anDisGam[g_nDisLvlMax]; // number of gammas known
int    g_anDisGamToLvl[g_nDisLvlMax][g_nDisLvlGamMax]; // daughter lvls
double g_adDisGamBR [g_nDisLvlMax][g_nDisLvlGamMax]; // daughter branching ratio
double g_adDisGamICC[g_nDisLvlMax][g_nDisLvlGamMax]; // daughter Alpha ICC
// could do these dynamically, but takes time to code

double g_dECrit; // trust lvl scheme up to this E, determined by g_nDisLvlMax
void ReadTALYS_DisInputFile() {
  ifstream lvlFile;
  lvlFile.open(zFile);

  string sNucSearchLine;
  int nZ, nA, nLineTot, nLvlTot = 0;
  bool bFoundNuc = false;
  int nTries = 1e6, nTry = 0;
  while(!bFoundNuc && nTry < nTries) {
    nTry++;
    getline(lvlFile, sNucSearchLine); // Header
    istringstream issHeader(sNucSearchLine);
    issHeader >> nZ >> nA >> nLineTot >> nLvlTot;
    if(nA == g_dAMass && nZ == g_nZ) 
      {bFoundNuc = true; cout << "Nucleus: " << endl << sNucSearchLine << endl;}
  }
  if(nLvlTot < 2) { cerr << "err: No levels, check z file" << endl; cin.get(); }

  for(int lvl=0; lvl<g_nDisLvlMax; lvl++) {
    string sLvlLine;
    getline(lvlFile,sLvlLine);
    istringstream issLvl(sLvlLine);
    int nLvl, nLvlPar, nLvlGam;
    double dLvlEner, dLvlSp, dLvlT12;
    string sJP;
    issLvl >> nLvl >> dLvlEner >> dLvlSp >> nLvlPar >> nLvlGam >> dLvlT12
      >> sJP; 

    if(nLvlPar == -1) nLvlPar = 0; // 1=+, 0=- diff convention dicebox and talys
    if(nLvl != lvl || lvl > nLvlTot) cerr << "err: File mismatch" << endl; 
    if(int(dLvlT12) == dLvlT12) dLvlT12 = 999; // sometimes no halflife meas

    g_adDisEne[lvl] = dLvlEner;
    g_adDisSp[lvl]  = dLvlSp;
    g_anDisPar[lvl] = nLvlPar;
    g_anDisGam[lvl] = nLvlGam;
    g_adDisT12[lvl] = dLvlT12 * 1e15; // fs

    for(int gam=0; gam<nLvlGam; gam++) {
      string sGamLine;
      getline(lvlFile,sGamLine);
      istringstream issGam(sGamLine);
      int nGamToLvl;
      double dGamBR, dGamICC;
      issGam >> nGamToLvl >> dGamBR >> dGamICC; 

      g_anDisGamToLvl[lvl][gam] = nGamToLvl;
      g_adDisGamBR   [lvl][gam] = dGamBR;
      g_adDisGamICC  [lvl][gam] = dGamICC;
    } // gam
  } // lvl
  lvlFile.close();

  g_dECrit = g_adDisEne[g_nDisLvlMax-1];   
  #ifdef bForceBinNum
  g_dConESpac = (g_dExIMax - g_dECrit) / double(g_nConEBin);
  #else
  g_nConEBin = int((g_dExIMax - g_dECrit) / g_dConESpac) + 1; // 1 beyond
  #endif
} // Read Input

// DICEBOX-like input
void ReadDICEBOX_DisInputFile() {
  ifstream lvlFile;
  lvlFile.open(zFile);

  string sHeadLine;
  getline(lvlFile, sHeadLine); // header
  int nLvlTot = 0;
  istringstream issHeader(sHeadLine);
  issHeader >> nLvlTot;
  if(nLvlTot < 2) { cerr << "err: No levels, check z file" << endl; cin.get(); }
  if(nLvlTot-1 != g_nDisLvlMax){cerr << "err: Lvl mismatch" << endl; cin.get();}

  for(int lvl=0; lvl<g_nDisLvlMax; lvl++) {
    string sLvlLine;
    getline(lvlFile,sLvlLine);
    istringstream issLvl(sLvlLine);
    int nLvl, nLvlPar, nLvlGam;
    double dLvlEner, dLvlSp, dLvlT12;
    string sJP;
    issLvl >> dLvlEner >> dLvlSp >> nLvlPar >> nLvlGam;

    if(nLvlPar == 1) nLvlPar = 0; // 1=+, 0=- opposite convention for RAINIER
      else nLvlPar = 1;
    if(lvl >= nLvlTot) { cerr << "err: File mismatch" << endl; cin.get(); }
    dLvlT12 = 999; // no halflife meas

    g_adDisEne[lvl] = dLvlEner;
    g_adDisSp [lvl] = dLvlSp;
    g_anDisPar[lvl] = nLvlPar;
    g_anDisGam[lvl] = nLvlGam;
    g_adDisT12[lvl] = dLvlT12; // unused

    double dIntenTot = 0.0;
    for(int gam=0; gam<nLvlGam; gam++) {
      string sGamLine;
      getline(lvlFile,sGamLine);
      istringstream issGam(sGamLine);
      double dGamToLvlE, dInten, dIntenErr, dGamToSp, dGamToPar, dGamDel,
        dGamICC;
      issGam >> dGamToLvlE >> dInten >> dIntenErr >> dGamToSp >> dGamToPar 
        >> dGamDel >> dGamICC; 

      // search for gam - this is why I prefer talys input
      for(int lvlsearch=0; lvlsearch<lvl; lvlsearch++) {
        if( TMath::Abs(g_adDisEne[lvlsearch] - dGamToLvlE) < 0.001) { // 1 keV match
          g_anDisGamToLvl[lvl][gam] = lvlsearch;
        } // lvl E match
      } // lvl search

      g_adDisGamICC[lvl][gam] = dGamICC;
      g_adDisGamBR[lvl][gam] = dInten * (1 + dGamICC); // needs to be normalized
      dIntenTot += g_adDisGamBR[lvl][gam];
    } // gam

    // normalize - this is also why I prefer talys input
    for(int gamnorm=0; gamnorm<nLvlGam; gamnorm++) {
      g_adDisGamBR[lvl][gamnorm] /= dIntenTot;
    } // norm

    // decay isomers to g.s. - this is also also why I prefer talys input
    if( (nLvlGam == 0) && (lvl != 0) ) {
      g_anDisGam[lvl] = 1;
      g_anDisGamToLvl[lvl][0] = 0;
      g_adDisGamBR[lvl][0] = 1.0;
    }

  } // lvl
  lvlFile.close();

  g_dECrit = g_adDisEne[g_nDisLvlMax-1];
  #ifdef bForceBinNum
  g_dConESpac = (g_dExIMax - g_dECrit) / double(g_nConEBin);
  #else
  g_nConEBin = int((g_dExIMax - g_dECrit) / g_dConESpac) + 1; // 1 beyond
  #endif
} // Read Input

void PrintDisLvl() {
  cout << "****** Discrete ******" << endl;
  for(int lvl=0; lvl<g_nDisLvlMax; lvl++) {
    // levels
    cout << lvl << ":\t " << g_adDisEne[lvl] << "   " << g_adDisSp[lvl]
      << (g_anDisPar[lvl]==1?"+":"-")  // careful dicebox flips parity
      << " "; 
      if(g_adDisT12[lvl] < 1e9) cout << g_adDisT12[lvl] << " fs" << endl;
      else cout << "N/A" << endl; // lifetime not measured
    // gammas
    for(int gam=0; gam<g_anDisGam[lvl]; gam++) {
      cout << "   " << g_anDisGamToLvl[lvl][gam] << " : " 
        << g_adDisGamBR[lvl][gam] << endl;
    } // gam
  } // lvl
} // Print Discrete

TH2D *g_h2PopDist;
#ifdef bExFullRxn
/////////////////////// TALYS Rxn Population File //////////////////////////////
const int g_nExPopI = 71; // bins 0-70
const int g_nSpPopIBin = 10; // spins 0-9, havent tested with half-int yet
void ReadPopFile() {
  cout << "Reading Population File" << endl;
  // TALYS File:
  //   projectile p
  //   element fe
  //   mass 56
  //   energy 16
  //   outpopulation y
  //   ldmodel 2
  //   spincutmodel 1
  //   bins 65
  //   maxlevelstar 5
  // copy the "Population of Z= 26 N= 30 ( 56Fe) before decay" section to 
  // "Fe56Pop.dat" in this directory

  double adEx [g_nExPopI];
  double adPop[g_nExPopI][g_nSpPopIBin];

  ifstream filePop;
  filePop.open(popFile);
  string sLine;
  getline(filePop, sLine); // header
  getline(filePop, sLine); // header

  while( getline(filePop, sLine) ) {
    istringstream issPop(sLine);
    int nBin;
    double dPopTot, dEx;
    issPop >> nBin >> dEx >> dPopTot;
    adEx[nBin] = dEx; // iss cant input directly into array elements
    for(int s=0; s<g_nSpPopIBin; s++) {
      double dPop;
      issPop >> dPop;
      adPop[nBin][s] = dPop;
    } // read J     
  } // read E

  g_h2PopDist = new TH2D("h2PopDist","h2PopDist", 
    g_nSpPopIBin,0,g_nSpPopIBin, g_nExPopI-1, adEx);

  for(int binE=1; binE<g_nExPopI; binE++) {
    for(int binJ=1; binJ<=g_nSpPopIBin; binJ++) {
      g_h2PopDist->SetBinContent(binJ,binE,adPop[binE-1][binJ-1]);
    } // assign J
  } // assign E

} // ReadPopFile
#endif // talys input for Oslo Analysis

////////////////////////// Level Density ///////////////////////////////////////
double GetEff(double dEx) {
  #ifdef bLD_BSFG
  double dEff = dEx - g_dE1;
  #endif
  #ifdef bLD_CTM
  double dEff = dEx - g_dE0;
  #endif
  if(dEff < 0.0) dEff = 0.00000001;
  return dEff;
} // GetEff

double GetLDa(double dEx) { // TALYS 1.8 asymptotic dependence
  double dEff = GetEff(dEx);
  #ifdef bLDaEx
  return g_dLDaAsym * (1 + g_dShellDelW * (1 - exp(-g_dDampGam * dEff)) / dEff);
  #else
  return g_dLDa;
  #endif
} // GetLDa

double GetSpinCut2(double dEx) {
  double dEff = GetEff(dEx);
  double dLDa = GetLDa(dEx);

  #ifdef bJCut_VonEgidy05 // Von Egidy PRC72,044311(2005)
  double dSpinCut2 = 0.0146 * pow(g_dAMass, 5.0/3.0)
  * (1 + sqrt(1 + 4 * dLDa * dEff)) / (2 * dLDa);
  #endif

  #ifdef bJCut_SingPart // Gholami PRC75,044308(2007)
  double dSpinCut2 = 0.1461 * sqrt(dLDa * dEff) * pow(g_dAMass, 2.0/3.0);
  #endif

  #ifdef bJCut_RigidSph // Grimes PRC10(1974) 2373-2386
  double dSpinCut2 = 0.0145 * sqrt(dEff / dLDa) * pow(g_dAMass, 5.0/3.0);
  #endif

  #ifdef bJCut_VonEgidy09 // Von Egidy PRC80,054310(2009)
  // empirical fit to other data with only mass table parameters
  double dSpinCut2 = 0.391 * pow(g_dAMass, 0.675) * pow(dEx - 0.5 * g_dDeuPair,0.312);
  #endif

  #ifdef bJCut_Other // choose/add what you like, these are some I've found:
  // everyone and their mother seems to have a favorite spin cutoff model

  //double dSpinCut2 = pow(0.83 * pow(g_dAMass,0.26),2); // TALYS 1.8 global

  //double dSpinCut2 = pow(0.98 * pow(g_dAMass,0.29),2); // DICEBOX CTM

  double dSpinCut2 = 0.01389 * pow(g_dAMass, 5.0/3.0) / g_dLDaAsym
    * sqrt(GetLDa(dEx) * dEff); // TALYS 1.8 default

  //double dSpinCut2 = 0.0888*sqrt(dLDa * dEff) * pow(g_dAMass,2/3.0); //DICEBOX
  #endif
  return dSpinCut2;
}

double GetDensityTot(double dEx) {
  double dEff = GetEff(dEx);

  #ifdef bLD_CTM // Constant Temperture Function Model
  double dEnDen = exp(dEff / g_dTemp) / g_dTemp;
  #endif
  #ifdef bLD_BSFG // Back shifted Fermi Gas
  double dLDa = GetLDa(dEx);
  double dSpinCut2 = GetSpinCut2(dEx);
  double dEnDen = 1.0 / (12.0 * sqrt(2) * sqrt(dSpinCut2) * pow(dLDa, 0.25)
    * pow(dEff, 5.0 / 4.0) ) * exp(2.0 * sqrt(dLDa * dEff) );
  #endif
  return dEnDen;
}

double GetDensity(double dEx, double dSp, int nPar) {
  double dEnDen = GetDensityTot(dEx);
  double dSpinCut2 = GetSpinCut2(dEx);


  double dSpDen = (dSp + 0.5) * exp(-pow(dSp + 0.5, 2) / (2 * dSpinCut2))
   / dSpinCut2;
  
  #ifdef bPar_Equipar
  int nParity = nPar; // basically unused in this model
  double dParDen = 0.5;
  #endif
  #ifdef bPar_Edep
  // Al-Quraishi PRC67,015803(2003)
  double dParDen;
  double dExpTerm = 1.0 / (1.0 + exp( g_dParC * (dEx - g_dParD) ) );
  if(g_bIsEvenA) {
    if(nPar == 1){ // positive parity
      dParDen = 0.5 * (1.0 + dExpTerm ); // lot of positive states low E
    } else { // negative parity
      dParDen = 0.5 * (1.0 - dExpTerm ); // not many negative states low E
    } // +/-
  } else { // odd A
    if(nPar == 0){ // negative parity
      dParDen = 0.5 * (1.0 + dExpTerm );
    } else { // postive parity
      dParDen = 0.5 * (1.0 - dExpTerm );
    } // +/-
  }
  #endif

  double dDenEJP = dEnDen * dSpDen * dParDen;
  return dDenEJP;
} // plot with: TF2 *fDen2 = new TF2("fDen2","GetDensity(x,y,1)",0,16,0,9); fDen2->Draw("colz")

//////////////////// Build Nucleus /////////////////////////////////////////////
// spins marked with postscript "b" refer to "bin": necessary for half-int spins
// not to be confused with boolean type prescript "b"
// b bin:    0    1    2    3    4    ...
// sp int:   0    1    2    3    4    ...
// half-int: 0.5  1.5  2.5  3.5  4.5  ...

int g_nConLvlTot; // # lvl in constructed scheme. includes multiple in same bin
int g_nConMaxLvlBin; // the largest number of levels in an EJP bin

// * = memory to be dynamically allocated
double *g_adConExCen; // centroid energies of constructed bins
int *g_anConLvl; // number of levels in an EJP bin
int *g_anConCumul; // cumulative levels for random width seeding

int EJP(int ex, int spb, int par) { // index for Energy E, Spin J, Parity P
  return ex + spb * g_nConEBin + par * g_nConEBin * g_nConSpbMax;
} // EJP

double GetInBinE(int nReal, int nConEx, int nSpb, int nPar, int nLvlInBin) {
  TRandom2 ranInBinE(1 + // have seen issues with 0 seed
    nPar + 
    nConEx    * 2 + 
    nReal     * 2 * g_nConEBin +
    nSpb      * 2 * g_nConEBin * g_nReal + 
    nLvlInBin * 2 * g_nConEBin * g_nReal * g_nConSpbMax);
  double dInBinE = g_adConExCen[nConEx];
  dInBinE += (-g_dConESpac / 2.0 + ranInBinE.Uniform(g_dConESpac));
  return dInBinE;
} // GetInBinE

const double g_dPi = 3.14159265359;
void BuildConstructed(int nReal) {
  // get rid of previous allocations
  cout << "Constructing Lvl Scheme" << endl;
  delete[] g_adConExCen;
  delete[] g_anConLvl;
  delete[] g_anConCumul;
  g_adConExCen = new double[g_nConEBin]; // center energy of the bin
  g_anConLvl   = new int   [g_nConEBin * g_nConSpbMax * 2]; 
  g_anConCumul = new int   [g_nConEBin * g_nConSpbMax * 2]; 

  TRandom2 ranLvlGen(1 + nReal);
  g_nConLvlTot = g_nDisLvlMax; // include discrete lvls below
  g_nConMaxLvlBin = 0; // max of any bin

  #ifdef bPoisson
  for(int ex=0; ex<g_nConEBin; ex++) {
    g_adConExCen[ex] = g_dECrit + (ex + 0.5) * g_dConESpac; 
    // bottom of bin ex=0 equals g_dECrit: no gap or overlap with discrete
    // constructed in-bin-energies are discretized by gaussian seeds

    for(int spb=0; spb<g_nConSpbMax; spb++) {
      for(int par=0; par<2; par++) {
        
        double dSp; // integer or half integer determination for density
        if(g_bIsEvenA) dSp = spb; else dSp = spb + 0.5;

        double dAvgNumLev = g_dConESpac * GetDensity(g_adConExCen[ex],dSp,par);
        int nRanNumLvl = ranLvlGen.Poisson(dAvgNumLev); // integer poisson
        g_anConLvl[EJP(ex,spb,par)] = nRanNumLvl; 
        if(nRanNumLvl > g_nConMaxLvlBin) g_nConMaxLvlBin = nRanNumLvl;

        g_nConLvlTot += g_anConLvl[EJP(ex,spb,par)];
        g_anConCumul[EJP(ex,spb,par)] = g_nConLvlTot;
      } // par
    } // sp bin
  } // ex
  #endif

  #ifdef bWigner
  for(int ex=0; ex<g_nConEBin; ex++) { // assign level energies, same as above
    g_adConExCen[ex] = g_dECrit + (ex + 0.5) * g_dConESpac;
  } // ex

  for(int par=0; par<2; par++) {
    for(int spb=0; spb<g_nConSpbMax; spb++) {
      double dSp; // integer or half integer, for density
      if(g_bIsEvenA) dSp = spb; else dSp = spb + 0.5;

      ///// initialize E bins to 0 /////
      for(int ex=0; ex<g_nConEBin; ex++) { 
        g_anConLvl[EJP(ex,spb,par)] = 0;
      } // ex

      ///// expected cumulative constructed # of lvls /////
      // Each JP has independent average energy bin spacing according to
      // the density inverse which is itself a function of energy
      double adExpCumulCon[g_nConEBin] = {0.0};
      for(int ex=0; ex<g_nConEBin; ex++) {
        double dEx = g_adConExCen[ex];
        if(ex == 0)
          adExpCumulCon[ex] = GetDensity(dEx, dSp, par) * g_dConESpac;
        else // ex != 0
          adExpCumulCon[ex] = GetDensity(dEx, dSp, par) * g_dConESpac 
            + adExpCumulCon[ex-1];
      } // ex

      ///// Level assignment /////
      double dWigSampleSum = 0.0; 
      // expectation value of avg dist between neighboring levels is 1
      for(int ex=0; ex<g_nConEBin; ex++) {
        while( dWigSampleSum < adExpCumulCon[ex] ) {
          dWigSampleSum += 2.0 / sqrt(g_dPi) 
            * sqrt( -log( ranLvlGen.Uniform(1.0) ) );
          g_anConLvl[EJP(ex,spb,par)]++;
        } // WigSampleSum < ExpCumulCon
        g_nConLvlTot += g_anConLvl[EJP(ex,spb,par)];
        g_anConCumul[EJP(ex,spb,par)] = g_nConLvlTot;
      } // ex

    } // spb
  } // par
  #endif
} // BuildConstructed

void PrintConLvl() {
  int nSpbPrint = 9; // won't algin with double digit spins or bin content > 9
  cout << "****** Constructed ******" << endl;
  cout << "More levels exist at higher spins" << endl;
  cout << "Parity   ";
  for(int spb=0; spb<=nSpbPrint; spb++) { cout << "-" << " ";} cout << " ";
  for(int spb=0; spb<=nSpbPrint; spb++) { cout << "+" << " ";} cout << endl;
  cout << "Spin Bin ";
  for(int spb=nSpbPrint; spb>=0; spb--) { cout << spb << " ";} cout << " ";
  for(int spb=0; spb<=nSpbPrint; spb++) { cout << spb << " ";} cout << endl;

  cout << "E(MeV)   " << endl;
  for(int ex=0; ex<g_nConEBin; ex++) {
    cout << fixed << setprecision(3) << g_adConExCen[ex] << "    ";
    int par = 0;
    for(int spb=nSpbPrint; spb>=0; spb--) {
      cout << g_anConLvl[EJP(ex,spb,par)] << "|";
    } // sp bin
    cout << " ";
    par = 1;
    for(int spb=0; spb<=nSpbPrint; spb++) {
      cout << g_anConLvl[EJP(ex,spb,par)] << "|";
    } // sp bin
    // printing energy of every level in bin would be so many more lines 
    // but could be done with something like:
    // for(int lvl=0; lvl<g_anConLvl[EJP(ex,spb,par)]; lvl++) 
    //   cout << GetInBinE(ex,spb,par,lvl) << endl;
    cout << scientific << endl;
  } // ex 
  cout << "Total Number of Levels = " << g_nConLvlTot << endl;
} // PrintConLvl

/////////////////// Transistion Type ///////////////////////////////////////////
int GetTransType(int nSpbI, int nParI, int nSpbF, int nParF) {
  // only integer spins
  int nTransType = 0; 
    // 0 No trans; 1 Pure E1; 2 Mixed M1+E2; 3 Pure M1; 4 Pure E2
  int ndSpb = TMath::Abs(nSpbF - nSpbI);
  int ndPar = TMath::Abs(nParF - nParI);

  if(nSpbF < 0 || nSpbF >= g_nConSpbMax) return 0; // not possible 
  if(nSpbI < 0 ) cerr << "err: Negative input spins" << endl;

  // there are faster ways to compute transistion type than what is below, 
  // but the function is not called that many times 
  // so its better to be fully explicit for clarity - trust me u dont want gotos
  // the only real head-scratchers occur when spin bin 0 is involved

  if(g_bIsEvenA) {
    if(ndSpb == 0) { 
      if(nSpbF > 0 && nSpbI > 0) { 
        if(ndPar == 0) {
          nTransType = 2; // M1+E2
        } else {
          nTransType = 1; // E1
        }
      } else { 
        nTransType = 0; // no 0 -> 0 with gammas, ignore E0 Internal Conversion
      }
    } else if(ndSpb == 1) { 
      if(nSpbF > 0 && nSpbI > 0) { // triangle check
        if(ndPar == 0) { 
          nTransType = 2; // M1+E2
        } else { 
          nTransType = 1; // E1
        }
      } else {
        if(ndPar == 0) { // 0+ -> 1+ or 0- -> 1- or 1- -> 0- or 1+ -> 0+
          nTransType = 3; // M1 Pure
        } else { // 0+ -> 1- or 0- -> 1+ or 1+ -> 0- or 1- -> 0+
          nTransType = 1; // E1
        }
      }
    } else if(ndSpb == 2) {
      if(ndPar == 0) {
        nTransType = 4; // E2
      } else {
        nTransType = 0; // no M2,E3
      }
    } else { // ndSpb > 2 
      nTransType = 0; // no Octupole
    }
  } else { //////////////////// odd A //////////////////////////////////////////
    // 0 No trans; 1 Pure E1; 2 Mixed M1+E2; 3 Pure M1; 4 Pure E2
    if(ndSpb == 0) {
      if(nSpbF > 0 && nSpbI > 0) {
        if(ndPar == 0) {
          nTransType = 2; // M1+E2
        } else {
          nTransType = 1; // E1
        }
      } else {
        if(ndPar == 0) {
          nTransType = 3; // 1/2+ -> 1/2+ pure M1, no E2 via triangle condition
        } else {
          nTransType = 1; // 1/2+ -> 1/2- E1
        }
      }
    } else if(ndSpb == 1) { 
      if(nSpbF > 0 && nSpbI > 0) { // triangle check
        if(ndPar == 0) {
          nTransType = 2; // M1+E2
        } else {
          nTransType = 1; // E1
        }
      } else {
        if(ndPar == 0) { // 1/2+ -> 3/2+ can be quadrupole unlike 0+ -> 1-
          nTransType = 2; // M1+E2 
        } else { // 1/2+ -> 3/2- , etc.
          nTransType = 1; // E1
        }
      }
    } else if(ndSpb == 2) {
      if(ndPar == 0) {
        nTransType = 4; // E2
      } else {
        nTransType = 0; // no M2,E3
      }
    } else { // ndSpb > 2
      nTransType = 0; // no Octupole
    }
  } // A determination

  return nTransType;
} // GetTransType 

/////////////////////// Gamma Strength * Eg^(2L+1) /////////////////////////////
// Physical Constants
const double g_d4Pi2 = 39.4784176; // 4*pi^2
const double g_dKX1 = 8.673592583E-08; // mb^-1 MeV^-2;  = 1/(3*(pi*hbar*c)^2) 
const double g_dKX2 = 5.204155555E-08; // mb^-1 MeV^-2;  = 1/(5*(pi*hbar*c)^2)

double GetTemp(double dEx) {
  double dEff = GetEff(dEx);
  double dLDa = GetLDa(dEx);
  return sqrt(dEff / dLDa);
} // GetTemp

double GetStrE1(double dEx, double dEg) { 
  double dStr = 0.0;
  double dTemp = GetTemp(dEx - dEg);

  for(int set=0; set<g_nParE1; set++) { // sum over split dipoles in applicable
    double dGam = g_adGamE1[set] * (dEg*dEg + g_d4Pi2 * dTemp*dTemp) 
      / (g_adEneE1[set]*g_adEneE1[set]); // energy dependent width

    #ifdef bE1_GenLor // Kopecky and Uhl Gen Lorentzian; model 5 dicebox
    double dTerm1 = dEg * dGam 
      / ( pow( (dEg*dEg - g_adEneE1[set]*g_adEneE1[set]),2) + pow(dEg*dGam,2) );
    double dTerm2 = 0.7 * g_adGamE1[set] * g_d4Pi2 * dTemp*dTemp 
      / pow(g_adEneE1[set],5);
      // non zero limit, F=0.7 is the fermi liquid quasiparticle collision fact
    double dTerm = dTerm1 + dTerm2;
    #endif

    #ifdef bE1_KMF // Kadmenskij Markushev Furman
    double dTerm = 0.7 * g_adEneE1[set] * dGam
      / pow( dEg*dEg - g_adEneE1[set]*g_adEneE1[set], 2);
    #endif

    #ifdef bE1_KopChr // Kopecky Chrien
    double dTerm = dEg * dGam
      / ( pow( (dEg*dEg - g_adEneE1[set]*g_adEneE1[set]),2) + pow(dEg*dGam,2) );
    #endif

    #ifdef bE1_StdLor // Standard Lorentzian
    double dTerm = dEg * g_adGamE1[set]
      / ( pow(dEg*dEg - g_adEneE1[set]*g_adEneE1[set], 2)
      + pow(dEg*g_adGamE1[set], 2) );
    #endif    

    dStr +=  g_dKX1 * g_adSigE1[set] * g_adGamE1[set] * dTerm;
  } // E1 parameter set

  if(dStr < 0) {cerr << "err: Negative strength" << endl;}
  return dStr * pow(dEg,3); // Eg^(2L+1) so this in not formally gamma strength
} // GetStrE1

double GetStrM1(double dEg) {
  // Standard Lorentzian
  double dStr = 0.0;
  for(int set=0; set<g_nParM1; set++) {
    dStr += g_dKX1 * g_adSigM1[set] * dEg * g_adGamM1[set]*g_adGamM1[set]
      / ( pow(dEg*dEg - g_adEneM1[set]*g_adEneM1[set], 2) 
      + pow(dEg*g_adGamM1[set], 2) );
  } // M1 parameter set
    
  #ifdef bM1StrUpbend
  double dUpbend = g_dUpbendM1Const * exp(-g_dUpbendM1Exp * dEg);
  dStr += dUpbend;
  #endif
  if(dStr < 0) {cerr << "err: Negative strength" << endl;}
  return dStr * pow(dEg,3);
}

double GetStrE2(double dEg) {
  #ifdef bE2_StdLor
  // Standard Lorentzian
  double dStr = g_dKX2 * g_dSigE2 * g_dGamE2*g_dGamE2 
    / (dEg * (pow(dEg*dEg - g_dEneE2*g_dEneE2,2) + pow(dEg*g_dGamE2,2)));
  // divide by Eg so units work out: TALYS formula units don't work
  #endif

  #ifdef bE2_SingPart
  double dStr = g_dSpSigE2;
  #endif
  if(dStr < 0) {cerr << "err: Negative strength" << endl;}
  return dStr * pow(dEg,5);
}

double g_de = 2.7182818284590452353; // Euler's number
double g_dlam = 0.5; // need to convert Gamma dist to Chi2 dist
double GetRanChi2(TRandom2 &ran, double dnu) {
  // ROOT doesn't supply Chi-Squared Dist tied to a TRandom
  // http://pdg.lbl.gov/2013/reviews/rpp2013-rev-monte-carlo-techniques.pdf
  double dAcceptX = -1;
  bool bAccept = false;
  double dk = dnu / 2.0;
  if(dk == 1) { // exponential dist
    double du = ran.Uniform();
    bAccept = true;
    dAcceptX = -log(du);
  } else if(dk > 0 && dk < 1) { // pole at 0, could return 0 due to underflow
    double dv1 = (g_de + dk) / g_de;
    while(!bAccept) {
      double du1 = ran.Uniform();
      double du2 = ran.Uniform();
      double dv2 = dv1 * du1;
      if(dv2 <= 1) {
        double dx = pow(dv2, 1.0/dk);
        if(du2 <= exp(-dx) ) {
          bAccept = true;
          dAcceptX = dx;
        } // else restart with new u1, u2
      } else { // dv2 > 1
        double dx = -log((dv1 - dv2) / dk);
        if(du2 <= pow(dx, dk-1)) {
          bAccept = true;
          dAcceptX = dx;
        } // else restart with new u1, u2
      } // v2 condition
    } // accept
  } else if(dk > 1) { // closed like gaussian
    double dc = 3 * dk - 0.75;
    while(!bAccept) {
      double du1 = ran.Uniform();
      double dv1 = du1 * (1 - du1);
      double dv2 = (du1 - 0.5) * sqrt(dc / dv1);
      double dx = dk + dv2 - 1;
      if(dx > 0) {
        double du2 = ran.Uniform();
        double dv3 = 64 * dv1*dv1*dv1 * du2*du2;
        if(  (dv3 <= 1 - 2 * dv2*dv2 / dx)
          || (log(dv3) <= 2 * ( (dk -1) * log(dx / (dk - 1)) - dv2) ) ) {
          bAccept = true;
          dAcceptX = dx;
        } // v3 condtion
      } // x condtion
    } // accept
  } else {
    cerr << "err: negative chi2 degree freedom" << endl;
  } // k conditon
  if(dAcceptX < 0) cerr << "err: no ch2 random assigned" << endl;
  return dAcceptX / g_dlam;
} // GetRanChi2

double GetStr(double dEx, double dEg, int nTransType, double &dMixDelta2,
  TRandom2 &ranStr) { 
  // returns sum_XL( Str_XL * Eg^(2L+1) )
  dMixDelta2 = 0.0;
  switch(nTransType) {
    case 0: return 0.0; // save some flops for impossible transistions
    case 1: { // int nX = 0, nL = 1;
              #ifdef bWFD_PTD
              double dGaus = ranStr.Gaus();
              double dFluct = dGaus * dGaus;
              #endif
              #ifdef bWFD_nu
              double dFluct = GetRanChi2(ranStr, g_dNu);
              #endif
              #ifdef bWFD_Off
              double dFluct = 1.0;
              #endif
              return dFluct * GetStrE1(dEx, dEg);
            }
    case 2: { 
              // both decay branches fluctuate independently
              // int nX = 1, nL = 1;
              #ifdef bWFD_PTD
              double dGausM1 = ranStr.Gaus();
              double dFluctM1 = dGausM1 * dGausM1;
              #endif
              #ifdef bWFD_nu
              double dFluctM1 = GetRanChi2(ranStr, g_dNu);
              #endif
              #ifdef bWFD_Off
              double dFluctM1 = 1.0;
              #endif
              double dStrM1 = dFluctM1 * GetStrM1(dEg);

              // int nX = 0, nL = 2;
              #ifdef bWFD_PTD
              double dGausE2 = ranStr.Gaus();
              double dFluctE2 = dGausE2 * dGausE2;
              #endif
              #ifdef bWFD_nu
              double dFluctE2 = GetRanChi2(ranStr, g_dNu);
              #endif
              #ifdef bWFD_Off
              double dFluctE2 = 1.0; // PT fluct off
              #endif
              double dStrE2 = dFluctE2 * GetStrE2(dEg);
              
              // delta = <I1||E2||I0> / <I1||M1||I0>; I0,1 = init,final w.f.
              //       = sqrt(GamE2 / GamM1) 
              //       = sqrt(strE2 * Eg^5 / strM1 * Eg^3)
              dMixDelta2 = dStrE2 / dStrM1; // fact of Eg in GetStrXL
              return dStrM1 + dStrE2; 
            }
    case 3: { // int nX = 1, nL = 1;
              #ifdef bWFD_PTD
              double dGaus = ranStr.Gaus();
              double dFluct = dGaus * dGaus;
              #endif
              #ifdef bWFD_nu
              double dFluct = GetRanChi2(ranStr, g_dNu);
              #endif
              #ifdef bWFD_Off
              double dFluct = 1.0;
              #endif
              return dFluct * GetStrM1(dEg);
            }
    case 4: { // int nX = 0, nL = 2;
              #ifdef bWFD_PTD
              double dGaus = ranStr.Gaus();
              double dFluct = dGaus * dGaus;
              #endif
              #ifdef bWFD_nu
              double dFluct = GetRanChi2(ranStr, g_dNu);
              #endif
              #ifdef bWFD_Off
              double dFluct = 1.0;
              #endif
              return dFluct * GetStrE2(dEg);
            }
    default : cerr << "err: Invaild strength type" << endl; return 0.0;
  } // TransType
} // GetStr

////////////////////// Internal Conversion /////////////////////////////////////
double GetBrICC(double dEg, int nTransMade=1, double dMixDelta2=0.0) {
  int nSuccess = -7;
  int nReadLine = 0;  // BrIcc output changes based on input
  switch(nTransMade) {
    case 0: return 0.0; break;
    case 1: nSuccess = system(
      Form("./%s -Z %d -g %f -L E1 -w %s > oAlpha.briccs", 
      cbriccs, g_nZ, dEg*1000, g_sBrIccModel) ); nReadLine = 8; break; // MeV
    case 2: nSuccess = system(
      Form("./%s -Z %d -g %f -L M1+E2 -d %f -w %s > oAlpha.briccs", 
      cbriccs, g_nZ, dEg*1000, sqrt(dMixDelta2), g_sBrIccModel) ); 
      nReadLine = 11; break;
    case 3: nSuccess = system(
      Form("./%s -Z %d -g %f -L M1 -w %s > oAlpha.briccs", 
      cbriccs, g_nZ, dEg*1000, g_sBrIccModel) ); nReadLine = 8; break;
    case 4: nSuccess = system(
      Form("./%s -Z %d -g %f -L E2 -w %s > oAlpha.briccs", 
      cbriccs, g_nZ, dEg*1000, g_sBrIccModel) ); nReadLine = 8; break;
    default: cerr << "err: impossible transistion" << endl;
  } // transistion
  if(nSuccess) cerr << "err: BrIcc failure" << endl;
  ifstream fileAlpha("oAlpha.briccs");
  string sAlphaLine;
  for(int line=0; line<nReadLine; line++) { getline(fileAlpha,sAlphaLine); }
  getline(fileAlpha,sAlphaLine);
  double dAlpha;
  istringstream issAlpha(sAlphaLine);
  issAlpha >> dAlpha;
  return dAlpha;
} // GetBrICC

// unmixed conv coeff from BrICC, hist interpolation faster than graph
TH1D *g_hE1ICC;
TH1D *g_hM1ICC;
TH1D *g_hE2ICC;
// plot after run: g_hE1ICC->Draw()
void InitICC() {
  cout << "Initializing internal conversion coefficients: alphas" << endl;
  g_hE1ICC = new TH1D("g_hE1ICC","g_hE1ICC", g_nBinICC,g_dICCMin,g_dICCMax);
  g_hM1ICC = new TH1D("g_hM1ICC","g_hM1ICC", g_nBinICC,g_dICCMin,g_dICCMax);
  g_hE2ICC = new TH1D("g_hE2ICC","g_hE2ICC", g_nBinICC,g_dICCMin,g_dICCMax);
  for(int bin=1; bin<=g_nBinICC; bin++) {
    double dBinCenter = g_hE1ICC->GetBinCenter(bin);
    g_hE1ICC->SetBinContent(bin, GetBrICC(dBinCenter,1,0));
    g_hM1ICC->SetBinContent(bin, GetBrICC(dBinCenter,3,0));
    g_hE2ICC->SetBinContent(bin, GetBrICC(dBinCenter,4,0));
    cout << bin << " / " << g_nBinICC << "\r" << flush;
  } // bin
  int nSuccess = system("rm oAlpha.briccs");
  cout << endl;
} // InitICC

double GetICC(double dEg, int nTransMade=1, double dMixDelta2=0.0) {
  #ifdef bUseICC
  switch(nTransMade) {
    case 0: return 0; break;
    case 1: return g_hE1ICC->Interpolate(dEg); break;
    case 2: {
      double dM1 = g_hM1ICC->Interpolate(dEg);
      double dE2 = g_hE2ICC->Interpolate(dEg);
      return (dM1 + dMixDelta2 * dE2) / (1 + dMixDelta2); } break;
    case 3: return g_hM1ICC->Interpolate(dEg); break;
    case 4: return g_hE2ICC->Interpolate(dEg); break;
    default: cerr << "err: impossible transistion" << endl;
  }
  return 0.0;
  #else
  double dEg1 = dEg; // to avoid unused warnings
  int nTransMade1 = nTransMade; 
  double dMixDelta21 = dMixDelta2; 
  return 0.0;
  #endif
} // GetICC

///////////////// Calculate Widths /////////////////////////////////////////////
double GetWidth(int nExI, int nSpbI, int nParI, int nLvlInBinI, int nReal,
  double *adConWid, double *adDisWid, TRandom2 *arConState) {
  
  TRandom2 ranStr(1 + nReal + g_nReal * ( g_anConCumul[EJP(nExI,nSpbI,nParI)] 
    - g_anConLvl[EJP(nExI,nSpbI,nParI)] + nLvlInBinI )); // seed with lvl num

  double dTotWid = 0.0;
  // should not calculate widths out of a discrete state
  // already have BR, discrete width not used in TakeStep
  if(nExI >= 0) { // in constructed scheme
    #ifdef bExSingle
    double dExI;
    if( nExI == (g_nConEBin - 1) )
      dExI = g_adConExCen[nExI] + 0.5 * g_dConESpac; // start top of 1st bin
    else 
      dExI = g_adConExCen[nExI];
    #else
    //double dExI = g_adConExCen[nExI]; // fast, good approximation
    double dExI = GetInBinE(nReal, nExI, nSpbI, nParI, nLvlInBinI);
    #endif

    double dSpI;
    if(g_bIsEvenA){ dSpI = nSpbI; } else { dSpI = nSpbI + 0.5; }
    double dLvlSpac = 1.0 / GetDensity(dExI,dSpI,nParI);
    if(dLvlSpac <= 0){cerr << "err: Level spacing" << endl; dLvlSpac = 0.00001;}

    /////// to constructed scheme ///////
    for(int spb=nSpbI-2; spb<=nSpbI+2; spb++) { // dipole and quadrapole
      for(int par=0; par<2; par++) {
        int nTransType = GetTransType(nSpbI,nParI,spb,par);
        if(nTransType) { // possible
          for(int ex=0; ex<nExI; ex++) {
            // Could look at transisitions within same energy bin, 
            // but so tiny effect. Would have to GetInBinE() of both in and 
            // out state so that we don't accidentally go up in energy, 

            int nLvlTrans = g_anConLvl[EJP(ex,spb,par)]; // #lvls in final bin
            if(nLvlTrans) { // might save some flops by ignoring empties
  
              double dExF = g_adConExCen[ex];
                // could get precise E but takes comp time, hardly effects str
              double dEg = dExI - dExF;

              double dStr = 0.0;
              arConState[EJP(ex,spb,par)] = ranStr; // want same rand #s
                // backing out same in-bin widths in TakeStep
              for(int outlvl=0; outlvl<nLvlTrans; outlvl++) {
                double dMixDelta2; // for ICC
                double dStrTmp = // need to get delta2
                  GetStr(dExI, dEg, nTransType, dMixDelta2, ranStr);
                dStr += dStrTmp * (1.0 + GetICC(dEg,nTransType,dMixDelta2));
              } // outlvl
  
              adConWid[EJP(ex,spb,par)] = dStr * dLvlSpac;
              dTotWid += dStr * dLvlSpac;
            } // final bin has lvls
          } // ex
        } // possible
      } // par
    } // sp bin

    /////// to discrete ///////
    for(int lvl=0; lvl<g_nDisLvlMax; lvl++) {
      double dExF = g_adDisEne[lvl];
      int nParF   = g_anDisPar[lvl];
      double dEg = dExI - dExF;

      double dSpF = g_adDisSp[lvl];
      int nSpbF; // convert half-int to int for binned transistion type
      if(g_bIsEvenA) nSpbF = int(dSpF + 0.001); else nSpbF = int(dSpF - 0.499);
        // protection against double imprecision
      int nTransType = GetTransType(nSpbI,nParI,nSpbF,nParF);
      adDisWid[lvl] = 0.0; // dont want uninitialized
      if(nTransType) { // possible
        double dMixDelta2; // for ICC
        double dStrTmp = // need to get delta2
          GetStr(dExI, dEg, nTransType, dMixDelta2, ranStr);
        double dStr = dStrTmp * (1.0 + GetICC(dEg,nTransType,dMixDelta2));

        adDisWid[lvl] = dStr * dLvlSpac;
        dTotWid += dStr * dLvlSpac;
      } // possible
    } // discrete lvl
  } else cerr << "err: Requesting discrete width" << endl; 
  return dTotWid;
} // GetWidth

/////////////////////// Lifetime ///////////////////////////////////////////////
const double g_dHBar = 6.5821195e-7; // MeV fs
double GetDecayTime(double dTotWid, TRandom2 &ranEv) {
  double dLifeT = g_dHBar / dTotWid; // fs
  return ranEv.Exp(dLifeT);
}

////////////////////////// Take Step ///////////////////////////////////////////
bool TakeStep(int &nConEx, int &nSpb, int &nPar, int &nDisEx, int &nLvlInBin, 
  int &nTransMade, double &dMixDelta2, double dTotWid, int nReal,
  double *adConWid, double *adDisWid, TRandom2 *arConState, TRandom2 &ranEv) {
  // TakeStep will change variables with &
  // Initial variables not marked with "I" because they are also the final state
  // use ConEx instead of nEx to be explicit that this is for constructed region
    // wasn't a problem in GetWidth since constructed region was prerequisite

  int nToConEx;
  int nToDisEx = g_nDisLvlMax;
  int nToSpb = -7, nToPar = -7, nToLvlInBin = -7; // junk initializers
  double dToMixDelta2 = 0.0;

  /////// in constructed scheme ///////
  if(nConEx >= 0) {
    double dSp; 
    if(g_bIsEvenA){ dSp = nSpb; } else { dSp = nSpb + 0.5; }
    #ifdef bExSingle
    double dExI;
    if( nConEx == (g_nConEBin - 1) )
      dExI = g_adConExCen[nConEx] + 0.5 * g_dConESpac; // start top of 1st bin
    else
      dExI = g_adConExCen[nConEx];
    #else
    //double dExI = g_adConExCen[nExI];
    double dExI = GetInBinE(nReal, nConEx, nSpb, nPar, nLvlInBin);
    #endif

    if(dTotWid <= 0.0) return false; // dead state in constructed: no E2 options
    double dRanWid = ranEv.Uniform(dTotWid);
    double dWidCumulative = 0.0;
    bool bFoundLvl = false;

    ///// decay to discrete? /////
    for(int lvl=0; lvl<g_nDisLvlMax; lvl++) {
      dWidCumulative += adDisWid[lvl]; // already has ICC
      if(adDisWid[lvl] > 1e-4) cerr << "err: discrete width uninit" << endl;
      if(dWidCumulative >= dRanWid) {//once adds to more than dRanWid, it decays
        if(!bFoundLvl) { // for safety check
          bFoundLvl = true;
          nToConEx = -1;
          nToDisEx = lvl;
          if(g_bIsEvenA) nToSpb = int(g_adDisSp[lvl] + 0.001);
          else nToSpb = int(g_adDisSp[lvl] - 0.499);
          nToPar = g_anDisPar[lvl];
          nToLvlInBin = 0;
          nTransMade = GetTransType(nSpb,nPar,nToSpb,nToPar);
          lvl = g_nDisLvlMax; // break out of loop, for speed
        } // found lvl
      } // Cumulative >= Rand
    } // discrete lvl

    ///// decay to constructed scheme? /////
    for(int spb=nSpb-2; spb<=nSpb+2; spb++) { // dipole and quadrapole
      for(int par=0; par<2; par++) {
        int nTransType = GetTransType(nSpb,nPar,spb,par);
        if(nTransType != 0) {
          for(int ex=0; ex<nConEx; ex++) {
            double dConWid = adConWid[EJP(ex,spb,par)];
            if(dConWid > 1e-4) cerr << "err: con width uninit" << endl;
            dWidCumulative += dConWid;

            if(dWidCumulative >= dRanWid) {// once adds up to dRanWid, it decays
              if(!bFoundLvl) { // possibly already decayed to discrete
                bFoundLvl = true;
                nToConEx = ex;
                nToDisEx = g_nDisLvlMax;
                nToSpb = spb;
                nToPar = par;
                nTransMade = nTransType;
                // need to backtrack width sum and find out which individual 
                // level in EJP bin it decayed to since many levels in a bin 
                // each with random width according to PT distribution
                dWidCumulative -= dConWid;
                bool bFoundLvlInBin = false;

                double dLvlSpac = 1.0 / GetDensity(dExI,dSp,nPar);
                int nLvlTrans = g_anConLvl[EJP(ex,spb,par)];
                if(nLvlTrans) {
                  double dExF = g_adConExCen[ex];
                  double dEg = dExI - dExF;
                  TRandom2 ranStr = arConState[EJP(ex,spb,par)];

                  for(int outlvl=0; outlvl<nLvlTrans; outlvl++) {
                    double dMixDelta2Tmp;
                    double dStrTmp = // need to get delta2
                      GetStr(dExI, dEg, nTransType, dMixDelta2Tmp, ranStr);
                    double dStr = dStrTmp * (1.0 + 
                      GetICC(dEg,nTransType,dMixDelta2Tmp));
                    dWidCumulative += dStr * dLvlSpac;
                    if(dWidCumulative >= dRanWid) {
                      if(!bFoundLvlInBin) { // possibly decayed prev inbin lvl
                        bFoundLvlInBin = true;
                        nToLvlInBin = outlvl;
                        dToMixDelta2 = dMixDelta2Tmp;
                        outlvl = nLvlTrans; // break loop for speed
                      } // found lvl in bin
                    } // Cumulative >= Rand
                  } // outlvl

                  ex = nConEx; // break out of loops, for speed
                  spb = g_nConSpbMax;
                  par = 2;
                } // final bin has lvls
              } // found lvl
            } // Cumulative >= Rand
          } // ex
        } // possible
      } // par
    } // sp bin
  } ///// in constructed scheme /////
  else { /////// in discrete ///////
    nToConEx = -1;
    double dRanBR = ranEv.Uniform(1.0); // should all add up to 1
    double dBRCumulative = 0.0;
    bool bFoundLvl = false;
    for(int gam=0; gam<g_anDisGam[nDisEx]; gam++) {
      dBRCumulative += g_adDisGamBR[nDisEx][gam]; // already has ICC in it
      if(dBRCumulative >= dRanBR) { // once adds to more than dRanBR, it decays
        if(!bFoundLvl) { // for safety check
          bFoundLvl = true;
          nToConEx = -1;
          nToDisEx = g_anDisGamToLvl[nDisEx][gam];
          if(g_bIsEvenA) nToSpb = int(g_adDisSp[nToDisEx] + 0.001);
          else nToSpb = int(g_adDisSp[nToDisEx] - 0.499);
          nToPar = g_anDisPar[nToDisEx];
          nToLvlInBin = 0;
          nTransMade = GetTransType(nSpb,nPar,nToSpb,nToPar);
          gam = g_anDisGam[nDisEx]; // break loop for speed
        } // found lvl
      } // Cumulative >= Rand
    } // gam choices
  } ///// in discrete /////

  // this is where user errors usually turn up the most
  if(nToSpb < 0 || nToPar < 0 || nToLvlInBin < 0 ) { // error check
    cerr << endl << "err: JP: from " << endl 
         << nSpb << (nPar?"+":"-") << ", Lvl: " << nDisEx << ", ConBin: " 
         << nConEx << ";" << nLvlInBin << endl << "To " << endl
         << nToSpb << (nToPar?"+":"-") << ", Lvl: " << nToDisEx << ", ConBin: "
         << nToConEx << ";" << nToLvlInBin << endl 
         << "Likely branching ratios from file don't add to 1.000000" << endl 
         << "Check level " << nDisEx << " in " << "zFile" << " manually.";
         // else write code to normalize
    return false;
  } // err
  nConEx = nToConEx;
  nDisEx = nToDisEx;
  nSpb   = nToSpb;
  nPar   = nToPar;
  nLvlInBin = nToLvlInBin;
  dMixDelta2 = dToMixDelta2;
  return true;
} // TakeStep

TH2D *g_h2PopI;
///////////////////////// Initial Excitation ///////////////////////////////////
void GetExI(int &nExI, int &nSpbI, int &nParI, int &nDisEx, int &nLvlInBinI,
  TRandom2 &ranEv, double dExIMean, double dExISpread) {
  // GetExI will change above variables marked with &

  ///// DICEBOX-like initial state /////
  #ifdef bExSingle
  // one starting state
  #ifdef bForceBinNum
  nExI = g_nConEBin - 1;
  #else 
  double dExI = g_dExIMax;
  nExI = round( (dExI - g_dECrit) / g_dConESpac );
  #endif
  nSpbI = int(g_dSpI);
  nParI = g_dParI;
  nDisEx = g_nDisLvlMax;
  nLvlInBinI = 0;
  #endif

  ///// Beta-decay like selection of states /////
  #ifdef bExSelect
  double dRanState = ranEv.Uniform(1.0);
  double dBRSum = 0.0;
  for(int state=0; state<g_nStateI; state++) {
    dBRSum += g_adBRI[state];
    if(dBRSum > dRanState) {
      double dExI = g_adExI[state];
      nExI = round( (dExI - g_dECrit) / g_dConESpac ); 
      nSpbI = int(g_adSpI[state]);
      nParI = g_anParI[state];
      nDisEx = g_nDisLvlMax;
      nLvlInBinI = 0;
    } // BR > RanState
  } // state
  #endif

  #ifdef bExSpread
  ///// Energy spread /////
  nDisEx = g_nDisLvlMax; // start with all discrete levels as possibilities
  double dSpIMean = g_dSpIMean;
  bool bFoundSpin = false;
  nParI = ranEv.Integer(2);

  while(!bFoundSpin) {
    nSpbI = ranEv.Poisson(dSpIMean);

    int nAttempt = 0;
    int nMaxAttempt = 1000;
    bool bFoundLvl = false;
    while(!bFoundLvl && nAttempt < nMaxAttempt){ // dont pop not existent lvls
      nAttempt++;
      nExI = round( (dExIMean - g_dECrit + ranEv.Gaus(0.0, dExISpread) ) 
        / g_dConESpac);
      if(nExI > g_nConEBin) cerr << "err: ExI above constructed max" << endl;
      if(nExI < 0) cerr << "err: ExI below constructed scheme" << endl;

      int nLvlAvail = g_anConLvl[EJP(nExI,nSpbI,nParI)];
      if(nLvlAvail > 0) {
        nLvlInBinI = ranEv.Integer(nLvlAvail); 
        bFoundLvl = true;
        bFoundSpin = true;
      } // lvl avail
    } // attempts at finding a level at given spin

    // if is no lvl with the given spin in the level scheme, there is not
    // much more you can do than to repick a spin, throws off given spin dist

  } // found spin
  #endif // specific initial EJP range

  #ifdef bExFullRxn
  // Randomly selects EJP bin from input file population distribution
  // if no level in corresponding bin, searches nearby E bins
  // - not much more you can do when matching continuum and discrete physics
  double dSp = 0.0, dEx = 0.0;
  double dPopIntegral = g_h2PopDist->Integral(); // could calc outside 
  bool bLvlMatch = false;

  // find a lvl from the histogram:
  bool bEJLvlSuggested = false;
  // Dont g_h2PopDist->GetRandom2(dSp, dEx)! ROOT's GetRandom2() is real crap!
  double dRanPop = ranEv.Uniform(dPopIntegral);
  double dPopSum = 0.0;
  for(int binE=1; binE<=g_nExPopI-1; binE++) {
    for(int binJ=1; binJ<=g_nSpPopIBin; binJ++) {
      // Does not do any interpolation of TALYS hist - room for improvement
      double dPop = g_h2PopDist->GetBinContent(binJ,binE);
      dPopSum += dPop; // populate if PopSum > RanPop
      if(!bEJLvlSuggested && (dPopSum > dRanPop) ) {
        bEJLvlSuggested = true;
        dSp = binJ-1; // 1st bin is spin 0 or 0.5
        double dLowEBdy = g_h2PopDist->GetYaxis()->GetBinLowEdge(binE);
        double dUpEBdy  = g_h2PopDist->GetYaxis()->GetBinUpEdge(binE);
        // matching the constructed and discrete regions aint pretty:
        if(dLowEBdy < g_dECrit + 0.001) { // is discrete
          // need to be very careful of doublets and precision
          nExI = -1;
          bool bDisBinFound = false;
          for(int lvl=0; lvl<g_nDisLvlMax; lvl++) { // find discrete
            double dLvlE = g_adDisEne[lvl];
            if(TMath::Abs(dLowEBdy - dLvlE) < 0.001) {
              nDisEx = lvl;
              bDisBinFound = true;
              bLvlMatch = true;
            } // match E
          } // lvl
          if( !bDisBinFound ) cerr << "err: no discrete match" << endl;
          if(g_bIsEvenA) nSpbI = int(g_adDisSp[nDisEx] + 0.001);
          else nSpbI = int(g_adDisSp[nDisEx] - 0.499);
          nParI = g_anDisPar[nDisEx];
          nLvlInBinI = 0;
        } else { // is constructed
          dEx = dLowEBdy + ranEv.Uniform(dUpEBdy - dLowEBdy);
          // assume equipartition
          nParI = ranEv.Integer(2); // equal positive:1 and negative:0 
          nSpbI = int(dSp); // should work for both even and odd A
          nExI = round( (dEx - g_dECrit) / g_dConESpac);
          if(nExI > g_nConEBin) cerr << "err: ExI above constructed max"
            << endl;
          if(nExI < 0) cerr << "err: ExI below constructed scheme" << endl;

          // search in bin then nearby
          bool bPlacedLvl = false;
          while(!bPlacedLvl) {
            int nLvlAvail = g_anConLvl[EJP(nExI,nSpbI,nParI)];
            if(nLvlAvail > 0) { // dont pop lvls that dont exist
              bLvlMatch = true;
              bPlacedLvl = true;
              nLvlInBinI = ranEv.Integer(nLvlAvail); // any lvls in bin fine
              nDisEx = g_nDisLvlMax; // start with all discrete lvls
            } else { // random walk in energy space, avoid ECrit line
              int nEBinStep = ranEv.Integer(2);
              if(nEBinStep) { // 0=down; 1=up
                nExI++;
                if(nExI >= g_nConEBin) { // stay below maximum energy
                  nExI -= 2;
                } // ex max
              } else { // step down
                nExI--;
                if(nExI < 0 ) { // stay above ECrit
                  nExI += 2;
                } // ex min
              } // increase or decrease E
            } // lvl avail
          } // placed lvl
        } // dis or con

      } // suggestion found
    } // bin J
  } // bin E
  if(!bEJLvlSuggested) cerr << "err: lvl not suggested" << endl;
  if(!bLvlMatch) cerr << "err: lvl not found" << endl;  
  #endif // large swath of EJP
} // GetExI

TH2D *g_ah2PopLvl  [g_nReal][g_nExIMean]; 
TH1D *g_ahDisPop   [g_nReal][g_nExIMean];
TH1D *g_ahJPop     [g_nReal][g_nExIMean];
TH1D *g_ahTSC      [g_nReal][g_nExIMean][g_nDisLvlMax];
TH1D *g_ahPrim2    [g_nReal][g_nExIMean][g_nPrim2];
TH2D *g_ah2FeedTime[g_nReal][g_nExIMean];
TH2D *g_ah2ExEg    [g_nReal][g_nExIMean];
TH2D *g_ah21Gen    [g_nReal][g_nExIMean];
TH1D *g_ahGSpec    [g_nReal][g_nExIMean];
TGraph *g_grTotWidAvg       [g_nExIMean];

/******************************************************************************/
/**************************** MAIN LOOP ***************************************/
/******************************************************************************/
void RAINIER(TString sWriteFileName = "SaveRAINIER") {
  cout << "Starting RAINIER" << endl;
  TTimeStamp tBegin;
  TFile *fRAINIER = new TFile(sWriteFileName+".root","recreate",sWriteFileName);
  #ifdef bLvlTALYS
  ReadTALYS_DisInputFile();
  #endif
  #ifdef bLvlDICEBOX
  ReadDICEBOX_DisInputFile();
  #endif
  #ifdef bPrintLvl
  PrintDisLvl();
  #endif
  #ifdef bExFullRxn
  ReadPopFile();
  #endif
  #ifdef bUseICC
  InitICC();
  #endif
  g_h2PopI = new TH2D("g_h2PopI","Events Populated in Simulation",
    g_dPlotSpMax,0,g_dPlotSpMax, 900,0,g_dExIMax);

  for(int exim=0; exim<g_nExIMean; exim++) {
    g_grTotWidAvg[exim] = new TGraph(g_nReal); // bench
  }

  ///////// Realization Loop /////////
  for(int real=0; real<g_nReal; real++) {
    BuildConstructed(real); 
    #ifdef bPrintLvl
    PrintConLvl();
    #endif
    cout << "Realization " << real << endl;

    ///////// Initial Excitation loop /////////
    for(int exim=0; exim<g_nExIMean; exim++) {
      double dExIMean = g_adExIMean[exim];
      double dExISpread = g_dExISpread; // could make resolution dep on ExIMean
      #ifdef bExSpread
      cout << "  Initial Excitation Mean: " << dExIMean << " +- "
        << dExISpread << " MeV" << endl;
      #endif

      ///// Initialize Histograms /////
      g_ah2PopLvl[real][exim] = new TH2D(
        Form("h2ExI%dPopLvl_%d",exim,real), 
        Form("Population of Levels: %2.1f MeV, Real%d",dExIMean,real),
        2*g_dPlotSpMax, -g_dPlotSpMax, g_dPlotSpMax,
        g_dExIMax / g_dConESpac, 0, g_dExIMax);

      double dFeedTimeMax = (270-20) / (5.5-11.0) * dExIMean + 520; // fs
        // harder to pick out multistep decay at short times and low ExI
      int nFeedTimeBin = 300;
      g_ah2FeedTime[real][exim] = new TH2D(
        Form("h2ExI%dFeedTime_%d",exim,real), 
        Form("Feeding Levels: %2.1f MeV, Real%d",dExIMean,real), 
        g_nDisLvlMax,0,g_nDisLvlMax, nFeedTimeBin,0.0,dFeedTimeMax);

      int nBinEx = 300, nBinEg = 300; // mama, rhosigchi, etc. purposes
      g_ah2ExEg[real][exim] = new TH2D(
        Form("h2ExI%dEg_%d",exim,real),
        Form("E_{x,i} = %2.1f MeV vs. E_{#gamma}, Real%d",dExIMean,real),
        nBinEg,0,g_dExIMax*1000, nBinEx,0,g_dExIMax*1000);

      g_ah21Gen[real][exim] = new TH2D(
        Form("h2ExI%d1Gen_%d",exim,real),
        Form("E_{x,i} = %2.1f MeV 1st Generation, Real%d",dExIMean,real),
        nBinEg,0,g_dExIMax*1000, nBinEx,0,g_dExIMax*1000);

      double dEgMax = g_dExIMax;
      for(int dis=0; dis<g_nDisLvlMax; dis++) {
        double dDisEne = g_adDisEne[dis];
        double dDisSp  = g_adDisSp [dis];
        int nDisPar    = g_anDisPar[dis];
        if(nDisPar == 0) // negative parity
          g_ahTSC[real][exim][dis] = new TH1D(
            Form("hExI%dto%dTSC_%d",exim,dis,real), 
            Form("TSC to lvl %2.1f- %2.3f MeV, Real%d",
            dDisSp,dDisEne,real),
            g_nEgBinTSC,0.0,dEgMax);
        else // positive parity
          g_ahTSC[real][exim][dis] = new TH1D(
            Form("hExI%dto%dTSC_%d",exim,dis,real), 
            Form("TSC to lvl %2.1f+ %2.3f MeV, Real%d",
            dDisSp,dDisEne,real),
            g_nEgBinTSC,0.0,dEgMax);
      } // TSC to discrete lvl

      for(int prim2=0; prim2<g_nPrim2; prim2++) {
        // dont want to make into th2d and do projections later 
        // like I did with feedAnalysis, was too time coding time costly
        g_ahPrim2[real][exim][prim2] = new TH1D(
          Form("hExI%dPrim2_%d_%d",exim,prim2,real), 
          Form("Primary 2^{+}: %2.1f MeV, Real%d",dExIMean,real),
          g_nEgBinTSC,0.0,dEgMax);
      } // prim2

      g_ahGSpec[real][exim] = new TH1D(
        Form("hExI%dGSpec_%d",exim,real),
        Form("Gamma Spectrum: %2.1f MeV, Real%d",dExIMean,real),
        g_nEgBinTSC,0.00,dEgMax);

      g_ahDisPop[real][exim] = new TH1D(
        Form("hExI%dDisPop_%d",exim,real),
        Form("Discrete Populations: %2.1f MeV, Real%d",dExIMean,real),
        g_nDisLvlMax,0,g_nDisLvlMax);

      g_ahJPop[real][exim] = new TH1D(
        Form("hExI%dJPop_%d",exim,real),
        Form("Spin Initial Pop: %2.1f MeV, Real%d",dExIMean,real),
        int(g_dPlotSpMax),0,g_dPlotSpMax); // wont have this plot for half int J

      #ifdef bExSingle
      // save initial widths and rands so dont need to recompute
      double *adDisWid1     = new double  [g_nDisLvlMax]; 
      double *adConWid1     = new double  [g_nConEBin * g_nConSpbMax * 2](); 
      TRandom2 *arConState1 = new TRandom2[g_nConEBin * g_nConSpbMax * 2];

      TRandom2 ranEv1(1); // unused
      int nExI1,nSpbI1,nParI1,nDisEx1,nLvlInBinI1;
      GetExI(nExI1, nSpbI1, nParI1, nDisEx1, nLvlInBinI1, ranEv1, 
        1.0, 1.0); // unused
      int nConEx1 = nExI1, nSpb1 = nSpbI1, nPar1 = nParI1,
        nLvlInBin1 = nLvlInBinI1;            
      cout << "Getting 1st step widths" << endl;
      double dTotWid1 = GetWidth(
        nConEx1, nSpb1, nPar1, nLvlInBin1, real,
        adConWid1, adDisWid1, arConState1);
      cout << "Starting decay events" << endl;
      #endif

      #ifdef bParallel
      #pragma omp parallel // if resize cmd window while running - will stall
      #endif
      { // each active processor gets an allocated array
        // dont have to renew with each event
        //double adDisWid[g_nDisLvlMax]; // width to each discrete level
        double *adDisWid;
        adDisWid = new double[g_nDisLvlMax]; // wid to each discrete lvl
        double *adConWid; // width to each EJP bin (summed over in-bin lvls)
        adConWid   = new double  [g_nConEBin * g_nConSpbMax * 2](); // 0 init
        TRandom2 *arConState; // TRandom2 state for randoms
        arConState = new TRandom2[g_nConEBin * g_nConSpbMax * 2];

        int nEle = 0;
        #ifdef bParallel
        int nCount = 0;
        int nThreads = omp_get_num_threads();
        #pragma omp for 
        #endif
        ////////////////////////////////////////////////////////////////////////
        /////// EVENT LOOP /////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        for(int ev=0; ev<g_nEvent; ev++) {
          #ifdef bParallel
          nCount++;
          if( !( (nCount * nThreads) % g_nEvUpdate) )
            fprintf(stderr,"    approx: %d / %d \r", 
              nCount * nThreads, g_nEvent);
          #else
          if( !(ev % g_nEvUpdate) ) 
            cout << "    " << ev << " / " << g_nEvent << "\r" << flush;
          #endif
          TRandom2 ranEv(1 + real + ev * g_nReal);

          /////// initial state formation ///////
          int nExI,nSpbI,nParI,nDisEx,nLvlInBinI; // Initial state variables
          GetExI(nExI, nSpbI, nParI, nDisEx, nLvlInBinI, ranEv, // init vars
            dExIMean, dExISpread); // with these inputs from experiment
          int nConEx = nExI, nSpb = nSpbI, nPar = nParI,
            nLvlInBin = nLvlInBinI; // set initial to working state variables
          g_ahJPop[real][exim]->Fill(nSpbI);

          double dExI;
          if(nConEx < 0) { // in discrete
            dExI = g_adDisEne[nDisEx];
          } else { // in constructed scheme
            dExI = GetInBinE(real,nExI,nSpbI,nParI,nLvlInBinI);
          }
          g_h2PopI->Fill(nSpbI,dExI);

          double dTimeToLvl = 0.0;
          int nStep = 0;
          bool bIsAlive = true;
          double dEg1 = 0.0, dEg2 = 0.0; // steps for TSC
          bool bHadEle = false;
          /////// Decay till ground state ////////
          while(nDisEx > 0 && bIsAlive) { // nDisEx > g.s. (i.e. excited)
            // bIsAlive: could get stuck in constructed high spin state and have
              // no dipole or quadrupole decay option: effectively an isomer

            ///// pre decay /////
            double dExPre, dMixDelta2;
            if(nConEx < 0) { // in discrete
              dExPre  = g_adDisEne[nDisEx];
            } else { // in constructed scheme
              dExPre  = GetInBinE(real,nConEx,nSpb,nPar,nLvlInBin); 
            } // pre decay

            // Fill populations
            if(nPar == 1) g_ah2PopLvl[real][exim]->Fill(nSpb + 0.5, dExPre); 
            else g_ah2PopLvl[real][exim]->Fill(-nSpb - 0.5, dExPre); 
              // +-0.5 so can plot 0+ and 0-

            ///// during decay /////
            int nTransMade = 0; // want to know multipole and character for ICC
            if(nConEx < 0) { ///// in discrete /////
              double dTotWid = 0.0;
              double dLifeT = g_adDisT12[nDisEx] / log(2); // lifetime from file
              double dDecayTime = ranEv.Exp(dLifeT);
              dTimeToLvl += dDecayTime;
              bIsAlive = TakeStep( // no variable change if !bIsAlive
                nConEx, nSpb, nPar, nDisEx, nLvlInBin,
                nTransMade, dMixDelta2, dTotWid, real,
                adConWid, adDisWid, arConState, ranEv);
            } else { ///// in constructed scheme //////
              #ifdef bExSingle
              double dTotWid;
              if(nConEx == g_nConEBin - 1) { // use saved init width
                dTotWid = dTotWid1;
              } else { // not at intial state
                dTotWid = GetWidth(
                  nConEx, nSpb, nPar, nLvlInBin, real,
                  adConWid, adDisWid, arConState);
              } // Ex single
              #else
              double dTotWid = GetWidth(
                nConEx, nSpb, nPar, nLvlInBin, real,
                adConWid, adDisWid, arConState);
              #endif
              if(ev == 0 && nStep == 0) { // bench
                if(real != 0) {
                  double dOldAvg, dReal;
                  g_grTotWidAvg[exim]->GetPoint(real-1, dReal, dOldAvg);
                  double dNewAvg = (dOldAvg * real + dTotWid) / double(real+1);
                  g_grTotWidAvg[exim]->SetPoint(real, real, dNewAvg);
                } else { // 1st real
                  g_grTotWidAvg[exim]->SetPoint(real, real, dTotWid);
                } 
              } // bench
              #ifdef bExSingle
              if(nConEx == g_nConEBin - 1) { // use saved init widths and rands
                dTimeToLvl += GetDecayTime(dTotWid1, ranEv); 
                bIsAlive = TakeStep(
                  nConEx, nSpb, nPar, nDisEx, nLvlInBin,
                  nTransMade, dMixDelta2, dTotWid1, real,
                  adConWid1, adDisWid1, arConState1, ranEv);
              } else { // not at initial state
                dTimeToLvl += GetDecayTime(dTotWid, ranEv); 
                bIsAlive = TakeStep(
                  nConEx, nSpb, nPar, nDisEx, nLvlInBin,
                  nTransMade, dMixDelta2, dTotWid, real,
                  adConWid, adDisWid, arConState, ranEv);
              } // Ex Single
              #else
              dTimeToLvl += GetDecayTime(dTotWid, ranEv); 
              bIsAlive = TakeStep(
                nConEx, nSpb, nPar, nDisEx, nLvlInBin,
                nTransMade, dMixDelta2, dTotWid, real,
                adConWid, adDisWid, arConState, ranEv);
              #endif
            } // end of decay
            nStep++;

            ///// post decay /////
            if(bIsAlive) { // not stuck in isomeric state: emission ignored
              double dExPost;
              if(nConEx < 0) { // in discrete
                dExPost = g_adDisEne[nDisEx];
                // Fill low lying feeding time
                g_ah2FeedTime[real][exim]->Fill(nDisEx, dTimeToLvl);
                g_ahDisPop[real][exim]->Fill(nDisEx);
              } else { // in constructed scheme
                dExPost = GetInBinE(real,nConEx,nSpb,nPar,nLvlInBin); 
              }

              double dEg = dExPre - dExPost;

              ///// Internal Conversion /////
              double dICC = GetICC(dEg,nTransMade,dMixDelta2);
              double dProbEle = dICC / (1.0 + dICC);
              double dRanICC = ranEv.Uniform(1.0);
              bool bIsElectron = false;
              if(dProbEle > dRanICC) bIsElectron = true;

              if( !bIsElectron ) { // emitted gamma
                double dExRes = g_dExRes; // ~ particle resolution
                double dExDet = dExI + ranEv.Gaus(0.0, dExRes);
                g_ahGSpec[real][exim]->Fill(dEg);
                g_ah2ExEg[real][exim]->Fill(dEg*1000, dExDet*1000);
                if(nStep == 1)
                  g_ah21Gen[real][exim]->Fill(dEg*1000, dExDet*1000);

                // TSC spectra
                if(nStep == 1) dEg1 = dEg;
                if(nStep == 2) dEg2 = dEg;
                // Prim2 spectra
                if(nStep == 1 && nDisEx < g_nDisLvlMax) { // 1st step discrete
                  for(int prim2=0; prim2<g_nPrim2; prim2++) {
                    if(g_anPrim2[prim2] == nDisEx) { // primary is known 2+
                      g_ahPrim2[real][exim][prim2]->Fill(dEg); // neglect ICC
                    } // primary is known 2+
                  } // check if primary 
                } // 1st step to discrete
              } else { // was electron
                bHadEle = true; // at least one electron ruins TSC
                nEle++;
              } // ICC check

              // TSC spectra to specific states
              if(nStep == 2 && !bHadEle) {
                for(int dis=0; dis<g_nDisLvlMax; dis++) {
                  if(nDisEx == dis) {
                    g_ahTSC[real][exim][dis]->Fill(dEg1); 
                    g_ahTSC[real][exim][dis]->Fill(dEg2);
                  } // discrete match
                } // end on discrete
              } // 2 steps

            } // IsAlive

          } // no longer excited
        } //////////////////////////////////////////////////////////////////////
        ////////////////////////// EVENTS //////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        // deallocate memory
        delete[] adConWid;
        delete[] arConState;
        adConWid   = 0;
        arConState = 0;
        #ifdef bExSingle
        delete[] adConWid;
        delete[] arConState;
        adConWid   = 0;
        arConState = 0;
        #endif
        //cout << "    " << nEle << " internal conversions" << endl;
      } // parallel
      cout << endl << "    " << g_nEvent << " / " << g_nEvent << " processed" 
        << endl << endl;

      ///// plotting preferences /////
      TH2D *g_ah2Temp[] = {
        g_ah2PopLvl  [real][exim], 
        g_ah2FeedTime[real][exim], 
        g_ah2ExEg    [real][exim],  
        g_ah21Gen    [real][exim] };
      int nh2Temp = 4;
      for(int h=0; h<nh2Temp; h++) {
        g_ah2Temp[h]->GetXaxis()->SetTitleSize(0.055);
        g_ah2Temp[h]->GetXaxis()->SetTitleFont(132);
        g_ah2Temp[h]->GetXaxis()->SetTitleOffset(0.8);
        g_ah2Temp[h]->GetXaxis()->CenterTitle();

        g_ah2Temp[h]->GetYaxis()->SetTitleSize(0.055);
        g_ah2Temp[h]->GetYaxis()->SetTitleFont(132);
        g_ah2Temp[h]->GetYaxis()->SetTitleOffset(0.65);    
        g_ah2Temp[h]->GetYaxis()->CenterTitle();

        g_ah2Temp[h]->GetZaxis()->SetTitleSize(0.045);
        g_ah2Temp[h]->GetZaxis()->SetTitleFont(132);
        g_ah2Temp[h]->GetZaxis()->SetTitleOffset(-0.4);    
        g_ah2Temp[h]->GetZaxis()->CenterTitle();
        g_ah2Temp[h]->GetZaxis()->SetTitle("Counts");    
      } // hists
      g_ah2PopLvl  [real][exim]->GetXaxis()->SetTitle("J#Pi");
      g_ah2PopLvl  [real][exim]->GetXaxis()->SetNdivisions(20,0,0,kFALSE);
      g_ah2PopLvl  [real][exim]->GetYaxis()->SetTitle("E_{x} (MeV)");
      g_ah2FeedTime[real][exim]->GetXaxis()->SetTitle("Discrete Level Number");
      g_ah2FeedTime[real][exim]->GetYaxis()->SetTitle("Feeding Time (fs)");
      g_ah2ExEg    [real][exim]->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
      g_ah2ExEg    [real][exim]->GetYaxis()->SetTitle("E_{x,I} (MeV)");

      // TH1D:
      TH1D *g_ah1Temp[] = {
        g_ahGSpec [real][exim], 
        //g_ahTSC   [real][exim],
        g_ahDisPop[real][exim], 
        g_ahJPop  [real][exim],
        g_ahPrim2 [real][exim][0] // 1st determines axes in plot
      };
      int nh1Temp = 4;
      for(int h=0; h<nh1Temp; h++) {
        g_ah1Temp[h]->GetXaxis()->SetTitleSize(0.055);
        g_ah1Temp[h]->GetXaxis()->SetTitleFont(132);
        g_ah1Temp[h]->GetXaxis()->SetTitleOffset(0.8);
        g_ah1Temp[h]->GetXaxis()->CenterTitle();
        g_ah1Temp[h]->GetXaxis()->SetTitle("E_{#gamma} (MeV)");    

        g_ah1Temp[h]->GetYaxis()->SetTitleSize(0.055);
        g_ah1Temp[h]->GetYaxis()->SetTitleFont(132);
        g_ah1Temp[h]->GetYaxis()->SetTitleOffset(0.85);    
        g_ah1Temp[h]->GetYaxis()->CenterTitle();
        g_ah1Temp[h]->GetYaxis()->SetTitle("Counts");    
      } // hists

    } // Excitation mean
  } // realization

  fRAINIER->Write();
  gROOT->ProcessLine(".L Analyze.h");
  TTimeStamp tEnd;
  double dElapsedSec = double(tEnd.GetSec() - tBegin.GetSec());
  cout << "Elapsed RAINIER time: " << dElapsedSec << " sec" << endl;
} // main
