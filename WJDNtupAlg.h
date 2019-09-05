#ifndef WJDNTUP_WJDNTUPALG_H
#define WJDNTUP_WJDNTUPALG_H 1

//#include "WJDTriggers.h"

#include "AthAnalysisBaseComps/AthAnalysisAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "DataModel/DataVector.h"

#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
#include "InDetConversionFinderTools/VertexPointEstimator.h"
//Example ROOT Includes
//#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TVector3.h"

//triggers
//#include "AsgTools/AnaToolHandle.h"
#include "TrigDecisionTool/TrigDecisionTool.h"

//check these
#include "AthenaBaseComps/AthAlgorithm.h"
#include "TrkVKalVrtFitter/TrkVKalVrtFitter.h"
#include "TrkVKalVrtFitter/IVertexCascadeFitter.h"
#include "TrkVKalVrtFitter/VxCascadeInfo.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>


// GRL include(s):
#include "AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h"

class StoreGateSvc;

namespace Trk {class IVertexFitter;
    class IVertexFitter;
    class TrkVKalVrtFitter;
    class IVertexCascadeFitter;
    class VxCascadeInfo; //new and improved
/*class RungeKuttaPropagator;*/}
namespace InDet {class VertexPointEstimator;
                 class IInDetTrackSelectionTool;
               class IVertexCascadeFitter;}

namespace Trig {class TrigDecisionTool;}

class WJDNtupAlg: public ::AthAnalysisAlgorithm { 
 public: 
  WJDNtupAlg( const std::string& name, ISvcLocator* pSvcLocator );
  virtual ~WJDNtupAlg(); 

  ///uncomment and implement methods as required

                                        //IS EXECUTED:
  virtual StatusCode  initialize();     //once, before any input is loaded
  virtual StatusCode  beginInputFile(); //start of each input file, only metadata loaded
  //virtual StatusCode  firstExecute();   //once, after first eventdata is loaded (not per file)
  virtual StatusCode  execute();        //per event
  //virtual StatusCode  endInputFile();   //end of each input file
  //virtual StatusCode  metaDataStop();   //when outputMetaStore is populated by MetaDataTools
  virtual StatusCode  finalize();       //once, after all events processed
  

  ///Other useful methods provided by base class are:
  ///evtStore()        : ServiceHandle to main event data storegate
  ///inputMetaStore()  : ServiceHandle to input metadata storegate
  ///outputMetaStore() : ServiceHandle to output metadata storegate
  ///histSvc()         : ServiceHandle to output ROOT service (writing TObjects)
  ///currentFile()     : TFile* to the currently open input file
  ///retrieveMetadata(...): See twiki.cern.ch/twiki/bin/view/AtlasProtected/AthAnalysisBase#ReadingMetaDataInCpp


 private: 

  //Example algorithm property, see constructor for declaration:
  //int m_nProperty = 0;

  //Example histogram, see initialize method for registration to output histSvc
  TH1D* m_myHist = 0, *h_JpsiMass = 0, *h_JpsiZFitted = 0, *h_DsMass = 0;
  TH1D *H0[20];
  TH1D *HT[20];
  TH1D *HT_V[28];
  TH2D *H2T[10];
  TH1D *H1Ds[10];
  TH1D *H1J[10];
  TH2D *H2[5], *H2J[2], *H2Ds[5];
  
  //TVector3 v3PrimaryDs,v3PrimaryTrkPerigee,v3PrimaryDsP,v3PrimaryTrkP;
  int cm_count = 0, evt_count = 0, indet_count = 0, indet_primcut = 0,Ds_Candidates = 0, Ds_Candidates1 = 0,Ds_Candidates2 = 0, Jpsi_Candidates = 0, k = 0, d1 = 0, d2 = 0,
  muon_count = 0, Ds_flag = 1, T_nJpsi = 0, T_nDs = 1, T_nKaons = 0, i = 0, j = 0, T_nPi = 0, T_nPhi = 0;
  const double pi = std::acos(-1);
  double  Ds_mass1, Ds_mass2, DsX1, DsY1, DsZ1, DsX2, DsY2, DsZ2, PrimX, PrimY, PrimZ, PrimXE, PrimYE, PrimZE, dPrim_Jpsi_temp=65535, dPrim_Jpsi=65535,
  mPi = 139, mK = 493, mMu = 105, casChi2 = 0, JpsiPtCut=10, DsPtCut=10, DsdR = 0, T_pTmax = 0, T_pTmin = 65535, Iso_rel = 1;
  
  float Jpsi_z0;
  bool T_DsPass = true, T_Ds = false, T_Jpsi = false, T_WJD = false, Pass_ColCut = false, Pass_IsoCut = false;
  TLorentzVector vmu1,vmu2,vj,vk1,vk2,vd11,vd12,vd21,vd22,vds1,vds2,vw1,vw2,pre_vds,vseed,vseed_sub,vphi1,vphi2;
  TLorentzVector Tvmu1,Tvmu2,Tvj,Tvds, Tvds1, Tvds2, Tvds3, Tvseed, TvW, TvPhi, TvPhiK1, TvPhiK2;
  //bool T_Pi_exists = false;
  bool SeedTrackCut = false;
  float Jpsi_range = 10;
  float SeedpTMin = 1.5;
  float SeedpTMax = 2;
  float SeedJpsidRMin = 1.5;
  double dRSeedCone = 1;
  bool isMC = false;
  bool dofit = false;
  bool doTruth = false;
  bool doData = true;
  bool DsIsoCheck = true;
  bool Verbose = true;
  bool Iso_cut = true;
  bool doNtup = true;
  double Collinear_cut = 0.5;
  float chi2 = 11;
  TTree* tree=0;

  std::vector<std::string> *matchTrigNames=0, *TrigNames=0;
  std:: vector<int>     *TrigRes=0, *L1TrigRes=0;
  
  UInt_t          evtNum=0, runNum=0, lumiNum=0, TrkNum=0,trkSize=0;
  UInt_t itr=0;

  std::vector<float>   *TrackPx=0,*TrackPy=0,*TrackPz=0,*TrackEnergy=0,*TrackD0=0,*TrackCharge=0,*TrackPt=0,*TrackPhi=0,*TrackEta=0, *trfHits=0, *trDzVtx=0;
  std::vector<int>     *TrackNDF=0,*TrackChi2=0;


  std::vector<float>  *JpsiTruthPx=0, *JpsiTruthPy=0, *JpsiTruthPz=0, *JpsiTruthE=0;
  std::vector<float> *DsTruthPx=0, *DsTruthPy=0, *DsTruthPz=0, *DsTruthE=0;
  
  std::vector<float> *allpriVtxX=0, *allpriVtxY=0, *allpriVtxZ=0, *allpriVtxChiNorm=0,*allpriVtxChi = 0, *allpriVtxXE = 0,*allpriVtxYE = 0,*allpriVtxZE = 0;
  Float_t         priVtxX = 0, priVtxY = 0, priVtxZ = 0, priVtxChiNorm = 0, priVtxChi = 0, priVtxXE = 0, priVtxYE = 0, priVtxZE = 0;
 
  std::vector<float>  *JpsiTrackPx=0, *JpsiTrackPy=0, *JpsiTrackPz=0, *JpsiTrackE=0, *JMass=0;
  std::vector<float>  *DsTrackPx1=0, *DsTrackPy1=0, *DsTrackPz1=0, *DsTrackE1=0,*DsTrackPx2=0, *DsTrackPy2=0, *DsTrackPz2=0, *DsTrackE2=0;
  std::vector<double>  *DMass1=0,*DMass2=0;
  std::vector<int>     *DCharge1=0,*DCharge2=0,*KCHARGE=0;
  std::vector<float>   *DDecayVtxX=0,*DDecayVtxY=0,*DDecayVtxZ=0,*DDecayVtxXE=0,*DDecayVtxYE=0,*DDecayVtxZE=0;
  Float_t  trackZ=0;
  std::vector<double>  *DVtxC2=0;
  UInt_t          nJ=0, nDs1 = 0, nDs2 = 0;
  std::vector<float>   *JDecayVtxX=0,*JDecayVtxY = 0,*JDecayVtxZ = 0,*JVtxC2 =0,*JDecayVtxXE = 0,*JDecayVtxYE = 0,*JDecayVtxZE = 0;
  std::vector<float>   *coneJpsiKplus1=0,*coneJpsiKminus1=0,*coneJpsiPi1=0,*coneJpsiKplus2=0,*coneJpsiKminus2=0,*coneJpsiPi2=0;
  UInt_t          nMC=0;
  std::vector<int>     *JpsiIndex=0,*DsIndex=0;

  std::vector<int>     *MCPdgIdAll=0;
  std::vector<float>   *MCEtaAll=0, *MCPhiAll=0; 
  std::vector<int>     *MCCharge=0;

  std::vector<int>     *MCPdgId=0, *MCnDecay=0,*MCNDaughters=0,*MCNParent=0,*MCStatus=0;
  std::vector<float>   *MCPx=0,*MCPy=0,*MCPz=0,*MCPhi=0,*MCeta=0,*MCE=0,*MCPt=0,*MCP=0,*MCVx=0,*MCVy=0,*MCVz=0,*MCtheta=0;
  
  UInt_t          nMu=0;
  std::vector<float>   *allmuPx=0,*allmuPy=0,*allmuPz=0;
  std::vector<float>   *allmuCharge=0;
  std::vector<int>   *allmuType=0,*allmuNDF = 0, *allmuQual=0;
  std::vector<float>   *allmuD0 = 0, *allmuDz = 0, *allmuChi2 = 0, *allmufHits=0,*muDzVtx=0;
  TLorentzVector    mupP,  mumP;
  std::vector<float>   *mumPx = 0, *mumPy = 0, *mumPz = 0, *mupPx = 0,*mupPy=0,*mupPz=0;
  std::vector<float>   *muDz=0;
  std::vector<float>   *mumfChi2 = 0, *mupfChi2 = 0;
  std::vector<int>      *mupfNDF = 0,*mumfNDF = 0,*mupIdx=0,*mumIdx=0;
  std::vector<bool> *muFirstSource =0;
  std::vector<double>  *trToJpsiTagDZ =0;

  std::vector<float>   *KplusDsPx1=0,*KplusDsPy1=0,*KplusDsPz1=0,*KplusDsE1=0,*KminusDsPx1=0,*KminusDsPy1=0,*KminusDsPz1=0,*KminusDsE1=0,*PiDsPx1=0,*PiDsPy1=0,*PiDsPz1=0,*PiDsE1=0;  
  std::vector<float>   *KplusDsD01=0,*KminusDsD01=0,*PiDsD01=0;
  std::vector<float>   *KplusDsPx2=0,*KplusDsPy2=0,*KplusDsPz2=0,*KplusDsE2=0,*KminusDsPx2=0,*KminusDsPy2=0,*KminusDsPz2=0,*KminusDsE2=0,*PiDsPx2=0,*PiDsPy2=0,*PiDsPz2=0,*PiDsE2=0; 
  std::vector<float>   *KplusDsD02=0,*KminusDsD02=0,*PiDsD02=0;

  UInt_t          nW1=0,nW2=0;
  std::vector<double>  *WPxLor1=0,*WPyLor1=0,*WPzLor1=0,*WMassLor1=0,*WPxLor2=0,*WPyLor2=0,*WPzLor2=0,*WMassLor2=0,*WVtxC2=0;
  //std::vector<float>   *WDecayVtxX=0,*WDecayVtxY=0,*WDecayVtxZ=0;



  Int_t ev;
  //std::vector<TLorentzVector> *JpsiTrackV,*DsTruthV,*DsTrackV,*JpsiTruthV;
  std::vector<const xAOD::TrackParticle*> DsTracks;
  std::vector<double> JpsiMassVector;
  std::vector<double> DsMassVector;

  //Trk::VxCascadeInfo *CascadeInfo;
  ToolHandle < Trk::IVertexFitter > m_iVertexFitter;
  ToolHandle < InDet::VertexPointEstimator > m_vertexEstimator;
  ToolHandle < Trk::IVertexFitter > m_JpsiFitter;
  ToolHandle < Trk::IVertexFitter > m_DsFitter;
  ToolHandle < Trk::IVertexFitter > m_WFitter;
  //std::string m_DsContainerName; //!< Name of Ds container
  //std::string m_JspiContainerName;
  ToolHandle < InDet::IInDetTrackSelectionTool > m_selTool;
//  ToolHandle < Trk::RungeKuttaPropagator > m_Propagator;
//  ToolHandle < Trk::IVertexCascadeFitter > m_CascadeFitter;
  ToolHandle < Trk::TrkVKalVrtFitter > m_CascadeFitter; //new and improved
  ToolHandle <Trig::TrigDecisionTool> m_tdt;
  ToolHandle <IGoodRunsListSelectionTool> m_grlTool;

   //TTree* m_myTree = 0;
   //ToolHandle < InDet::VertexPointEstimator > m_vertexEstimator;

}; 

#endif //> !WJD_WJDNTUPALG_H
