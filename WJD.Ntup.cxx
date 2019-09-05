// WJD includes

#include "WJDNtupAlg.h"

#include "xAODEventInfo/EventInfo.h"

#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/VertexAuxContainer.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/Muon.h"

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthEventContainer.h"

//triggers
#include "TrigDecisionTool/TrigDecisionTool.h"

// Amg include
#include "EventPrimitives/EventPrimitives.h"
//check these


#include "TrkVertexFitterInterfaces/IVertexFitter.h"
#include "TrkVKalVrtFitter/TrkVKalVrtFitter.h"
#include "TrkVKalVrtFitter/IVertexCascadeFitter.h"
#include "TrkVKalVrtFitter/VxCascadeInfo.h"

#include "InDetConversionFinderTools/VertexPointEstimator.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IPartPropSvc.h"

#include "DataModel/ElementLink.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"



// GRL include(s):
#include "AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h"

#include "EventKernel/PdtPdg.h"

WJDNtupAlg::WJDNtupAlg( const std::string& name, ISvcLocator* pSvcLocator ) : AthAnalysisAlgorithm( name, pSvcLocator ),
  JpsiPtCut(10),
  DsPtCut(10),
  SeedTrackCut(false),
  Jpsi_range(10),
  SeedpTMin(1.5),
  SeedpTMax(2),
  SeedJpsidRMin(1.5),
  dRSeedCone(1),
  dofit(false),
  doTruth(true),
  doData(true),
  DsIsoCheck(false),
  Verbose(true),
  Iso_cut(true),
  doNtup(true),
  Collinear_cut(0.8),
  m_iVertexFitter("Trk::TrkVKalVrtFitter"),
  m_vertexEstimator("InDet::VertexPointEstimator"),
  m_JpsiFitter("Trk::TrkVKalVrtFitter"),
  m_DsFitter("Trk::TrkVKalVrtFitter"),
  m_WFitter("Trk::TrkVKalVrtFitter"),
  //m_DsContainer(0), m_JpsiContainer(0)
  m_selTool( "InDet::InDetTrackSelectionTool/TrackSelectionTool"), 
  //m_CascadeFitter("Trk::TrkVKalVrtFitter/IVertexCascadeFitter"),
  m_CascadeFitter("Trk::TrkVKalVrtFitter"), // new and improved
  m_tdt("Trig::TrigDecisionTool/TrigDecisionTool")
  //m_grlTool("GoodRunsLists/IGoodRunsListSelectionTool")
//  m_Propagator("Trk::RungeKuttaPropagator")
{
  //declareInterface<WJDNtupAlg>(this);
  //declareProperty( "Property", m_nProperty = 0, "My Example Integer Property" ); //example property declaration
  declareProperty("JpsiPtCut",JpsiPtCut);
  declareProperty("DsPtCut",DsPtCut);
  declareProperty("SeedTrackCut",SeedTrackCut);
  declareProperty("Jpsi_range",Jpsi_range);
  declareProperty("SeedpTMin",SeedpTMin);
  declareProperty("SeedpTMax",SeedpTMax);
  declareProperty("SeedJpsidRMin",SeedJpsidRMin);
  declareProperty("dRSeedCone",dRSeedCone);
  declareProperty("CascadeFit",dofit);
  declareProperty("doTruth",doTruth);
  declareProperty("doData",doData);
  declareProperty("DsIsoCheck",DsIsoCheck);
  declareProperty("Verbose",Verbose);
  declareProperty("Iso_cut",Iso_cut);
  declareProperty("doNtup",doNtup);
  declareProperty("Collinear_cut",Collinear_cut);
  declareProperty("VertexFitter",m_iVertexFitter);
  declareProperty("VertexPointEstimator",m_vertexEstimator);
  declareProperty("JpsiVertexFitter",m_JpsiFitter);
  declareProperty("DsVertexFitter",m_DsFitter);
  declareProperty("TrackSelectionTool", m_selTool );
  declareProperty("TrackCascadeFitter", m_CascadeFitter );
  //declareProperty("GRLTool",  m_grlTool );
  //declareProperty("JpsiCandidates"  ,m_JpsiContainerName = "JpsiCandidates");
  //declareProperty("DsCandidates"  ,m_DsContainerName = "DsCandidates");
}
//
WJDNtupAlg::~WJDNtupAlg() {}

StatusCode WJDNtupAlg::initialize() {
  ATH_MSG_INFO ("Initializing " << name() << "...");
  //
  //This is called once, before the start of the event loop
  //Retrieves of tools you have configured in the joboptions go here
  //
  // get the tool service
  IToolSvc* toolSvc;
  StatusCode sc = service("ToolSvc",toolSvc);
  if (StatusCode::SUCCESS != sc) {
    ATH_MSG_ERROR("Unable to retrieve ToolSvc");
    return StatusCode::FAILURE;
  }
 
  // get the VertexFitter tool
  if ( m_iVertexFitter.retrieve().isFailure() ) {
    ATH_MSG_FATAL("Failed to retrieve tool " << m_iVertexFitter);
    return StatusCode::FAILURE;
  } else {
    ATH_MSG_INFO("Retrieved tool " << m_iVertexFitter);
  }
 if (!dofit) {casChi2 = 0;}
  if (dofit) {
    JpsiMassVector.push_back(mMu);
    JpsiMassVector.push_back(mMu);
    DsMassVector.push_back(mPi);
    DsMassVector.push_back(mK);
    DsMassVector.push_back(mK);
  }




 // create a ROOT Tree

	TFile *f = new TFile("new.root","recreate");
    tree= new TTree("WJD_TTRee","WJD");
    CHECK( histSvc()->regTree("/MYSTREAM/TTRee",tree) );

 // branches for trigger
    tree->Branch("TrigNames", &TrigNames,16000);
    tree->Branch("MatchTrigNames", &matchTrigNames,16000);
    tree->Branch("L1TrigRes", &L1TrigRes,16000);
    tree->Branch("TrigRes", &TrigRes,16000);
 
 // branches for eventInfo
    tree->Branch("evtNum", &evtNum,16000);
    tree->Branch("runNum", &runNum,16000);
    tree->Branch("lumiNum", &lumiNum,16000);
    tree->Branch("TrkNum", &TrkNum,16000);
    
 
 // branches for Jpsi  
    tree->Branch("JpsiTrackPx", &JpsiTrackPx,16000);
    tree->Branch("JpsiTrackPy", &JpsiTrackPy,16000);
    tree->Branch("JpsiTrackPz", &JpsiTrackPz,16000);
    tree->Branch("JpsiTrackEnergy", &JpsiTrackE,16000);
    tree->Branch("JpsiMass", &JMass,16000);
    tree->Branch("nJ", &nJ,16000);
    tree->Branch("JDecayVtxX",&JDecayVtxX,16000);
    tree->Branch("JDecayVtxY",&JDecayVtxY,16000);
    tree->Branch("JDecayVtxZ",&JDecayVtxZ,16000);
    tree->Branch("JDecayVtxXE",&JDecayVtxXE,16000);
    tree->Branch("JDecayVtxYE",&JDecayVtxYE,16000);
    tree->Branch("JDecayVtxZE",&JDecayVtxZE,16000);
    tree->Branch("JVtxC2",&JVtxC2,16000);
    

    tree->Branch("priVtxX", &priVtxX,16000);
    tree->Branch("priVtxY", &priVtxY,16000);
    tree->Branch("priVtxZ", &priVtxZ,16000);
    tree->Branch("priVtxChiNorm",&priVtxChiNorm,16000);
    tree->Branch("priVtxChi",&priVtxChi,16000);
    tree->Branch("priVtxXE", &priVtxXE,16000);
    tree->Branch("priVtxYE", &priVtxYE,16000);
    tree->Branch("priVtxZE", &priVtxZE,16000);

    
    tree->Branch("allpriVtxX", &allpriVtxX,16000);
    tree->Branch("allpriVtxY", &allpriVtxY,16000);
    tree->Branch("allpriVtxZ", &allpriVtxZ,16000);
    tree->Branch("allpriVtxXE", &allpriVtxXE,16000);
    tree->Branch("allpriVtxYE", &allpriVtxYE,16000);
    tree->Branch("allpriVtxZE", &allpriVtxZE,16000);
    tree->Branch("allpriVtxChiNorm",&allpriVtxChiNorm,16000);
    tree->Branch("allpriVtxChi",&allpriVtxChi,16000);
    

    tree->Branch("TrackPx",&TrackPx,16000);
    tree->Branch("TrackPy",&TrackPy,16000);
    tree->Branch("TrackPz",&TrackPz,16000);  
    tree->Branch("TrackEnergy",&TrackEnergy,16000); 
    tree->Branch("TrackNDF",&TrackNDF,16000);   
    tree->Branch("TrackChi2",&TrackChi2,16000); 
    tree->Branch("TrackPt",&TrackPt,16000);   
    tree->Branch("TrackCharge",&TrackCharge,16000);
    tree->Branch("TrackPhi",&TrackPhi,16000);    
    tree->Branch("TrackEta",&TrackEta,16000); 
    //tree->Branch("trDzVtx",&trDzVtx);
    tree->Branch("TrackD0",&TrackD0,16000); 
    tree->Branch("trfHits",&trfHits,16000);
    //tree->Branch("trToJpsiTagDZ",&trToJpsiTagDZ);
    
    
    
     
 // branches for muons
    tree->Branch("allmuPx", &allmuPx,16000);
    tree->Branch("allmuPy", &allmuPy,16000);
    tree->Branch("allmuPz", &allmuPz,16000);
    tree->Branch("allmuCharge", &allmuCharge,16000);
    tree->Branch("allmuType", &allmuType,16000);
    tree->Branch("allmuD0", &allmuD0,16000);
    tree->Branch("allmuDz", &allmuDz,16000);
    tree->Branch("allmuChi2", &allmuChi2,16000);
    tree->Branch("allmuNDF",&allmuNDF,16000);
    tree->Branch("allmuQual",&allmuQual,16000);
    tree->Branch("allmufHits",&allmufHits,16000);
    //tree->Branch("muDzVtx",&muDzVtx);
    
    
    
    
    tree->Branch("mumPx", &mumPx,16000);
    tree->Branch("mumPy", &mumPy,16000);
    tree->Branch("mumPz", &mumPz,16000);
    tree->Branch("mumfChi2", &mumfChi2,16000);
    tree->Branch("mumfNDF", &mumfNDF,16000);
    tree->Branch("mupPx", &mupPx,16000);
    tree->Branch("mupPy", &mupPy,16000);
    tree->Branch("mupPz", &mupPz,16000);
    tree->Branch("mupfChi2", &mupfChi2,16000);
    tree->Branch("mupfNDF", &mupfNDF,16000);
    //tree->Branch("mupIdx", &mupIdx);
    //tree->Branch("mumIdx", &mumIdx);
    
    
  
  //branches for Ds
    tree->Branch("DsTrackPxtheory1", &DsTrackPx1,16000);
    tree->Branch("DsTrackPytheory1", &DsTrackPy1,16000);
    tree->Branch("DsTrackPztheory1", &DsTrackPz1,16000);
    tree->Branch("DsTrackEnergytheory1", &DsTrackE1,16000);
    tree->Branch("DsMasstheory1", &DMass1,16000);   
    tree->Branch("nDstheory1", &nDs1,16000);   
    
    tree->Branch("DDecayVtxX",&DDecayVtxX,16000);
    tree->Branch("DDecayVtxY",&DDecayVtxY,16000);
    tree->Branch("DDecayVtxZ",&DDecayVtxZ,16000);
    tree->Branch("DDecayVtxXE",&DDecayVtxXE,16000);
    tree->Branch("DDecayVtxYE",&DDecayVtxYE,16000);
    tree->Branch("DDecayVtxZE",&DDecayVtxZE,16000);
    tree->Branch("DVtxC2",&DVtxC2,16000);
    
  

    tree->Branch("DsTrackPxtheory2", &DsTrackPx2,16000);
    tree->Branch("DsTrackPytheory2", &DsTrackPy2,16000);
    tree->Branch("DsTrackPztheory2", &DsTrackPz2,16000);
    tree->Branch("DsTrackEnergytheory2", &DsTrackE2,16000);
    tree->Branch("DsMasstheory2", &DMass2,16000);
    tree->Branch("nDstheory2", &nDs2,16000);
    

    tree->Branch("KplusDsPx1",&KplusDsPx1,16000);
    tree->Branch("KplusDsPy1",&KplusDsPy1,16000);
    tree->Branch("KplusDsPz1",&KplusDsPz1,16000);
    tree->Branch("KplusDsD01",&KplusDsD01,16000);
    tree->Branch("KplusDsD01",&KplusDsE1,16000);
    tree->Branch("KminusDsPx1",&KminusDsPx1,16000);
    tree->Branch("KminusDsPy1",&KminusDsPy1,60001);
    tree->Branch("KminusDsPz1",&KminusDsPz1,16000);
    tree->Branch("KminusDsD01",&KminusDsD01,16000);
    tree->Branch("KminusDsD01",&KminusDsE1,16000);
    tree->Branch("PiDsPx1",&PiDsPx1,16000);
    tree->Branch("PiDsPy1",&PiDsPy1,16000);
    tree->Branch("PiDsPz1",&PiDsPz1,16000);
    tree->Branch("PiDsD01",&PiDsD01,16000);
    tree->Branch("PiDsE1",&PiDsE1,16000);
    tree->Branch("coneJpsiKplus1",&coneJpsiKplus1,16000);
    tree->Branch("coneJpsiKminus1",&coneJpsiKminus1,16000);
    tree->Branch("coneJpsiPi1",&coneJpsiPi1,16000);
    

    tree->Branch("KplusDsPx2",&KplusDsPx2,16000);
    tree->Branch("KplusDsPy2",&KplusDsPy2,16000);
    tree->Branch("KplusDsPz2",&KplusDsPz2,16000);
    tree->Branch("KplusDsD02",&KplusDsD02,16000);
    tree->Branch("KplusDsD02",&KplusDsE2,16000);
    tree->Branch("KminusDsPx2",&KminusDsPx2,16000);
    tree->Branch("KminusDsPy2",&KminusDsPy2,16000);
    tree->Branch("KminusDsPz2",&KminusDsPz2,16000);
    tree->Branch("KminusDsD02",&KminusDsD02,16000);
    tree->Branch("KminusDsD02",&KminusDsE2,16000);
    tree->Branch("PiDsPx2",&PiDsPx2,16000);
    tree->Branch("PiDsPy2",&PiDsPy2,16000);
    tree->Branch("PiDsPz2",&PiDsPz2,16000);
    tree->Branch("PiDsD02",&PiDsD02,16000);
    tree->Branch("PiDsE2",&PiDsE2,16000);
    tree->Branch("coneJpsiKplus2",&coneJpsiKplus2,16000);
    tree->Branch("coneJpsiKminus2",&coneJpsiKminus2,16000);
    tree->Branch("coneJpsiPi2",&coneJpsiPi2,16000);

    
    tree->Branch("nW1",&nW1,16000);
    tree->Branch("WMassLor1",&WMassLor1,16000);
    tree->Branch("WPxLor1",&WPxLor1,16000);
    tree->Branch("WPyLor1",&WPyLor1,16000);
    tree->Branch("WPzLor1",&WPzLor1,16000);
    tree->Branch("nW2",&nW2,16000);
    tree->Branch("WMassLor2",&WMassLor2,16000);
    tree->Branch("WPxLor2",&WPxLor2,16000);
    tree->Branch("WPyLor2",&WPyLor2,16000);
    tree->Branch("WPzLor2",&WPzLor2,16000);

    //tree->Branch("WDecayVtxX",&WDecayVtxX);
    //tree->Branch("WDecayVtxY",&WDecayVtxY);
    //tree->Branch("WDecayVtxZ",&WDecayVtxZ);
    //tree->Branch("WVtxC2",&WVtxC2);
    


  // branches for monte carlo  
    tree->Branch("MCNumber", &nMC,16000);
    tree->Branch("MCPdgIdAll", &MCPdgIdAll,16000);
    tree->Branch("MCEtaAll", &MCEtaAll,16000);
    tree->Branch("MCPhiAll", &MCPhiAll,16000);
    tree->Branch("MCCharge", &MCCharge,16000);
    tree->Branch("MCPdgId", &MCPdgId,16000);
    //tree->Branch("MCP", &MCP);     
    tree->Branch("MCPx", &MCPx,16000);
    tree->Branch("MCPy", &MCPy,16000);
    tree->Branch("MCPz", &MCPz,16000);
    tree->Branch("MCPt", &MCPt,16000);   
    tree->Branch("MCE", &MCE,16000); 
    tree->Branch("MCPhi", &MCPhi,16000);
    tree->Branch("MCeta", &MCeta,16000); 
    tree->Branch("MCtheta", &MCtheta,16000); 
    tree->Branch("MCStatus", &MCStatus,16000); 
    tree->Branch("MCVx", &MCVx,16000);
    tree->Branch("MCVy", &MCVy,16000);
    tree->Branch("MCVz", &MCVz,16000);
    tree->Branch("MCNDaughters", &MCNDaughters,16000);
    tree->Branch("MCNParent", &MCNParent,16000);

    



    tree->Branch("JpsiTruthPx", &JpsiTruthPx,16000);
    tree->Branch("JpsiTruthPy", &JpsiTruthPy,16000);
    tree->Branch("JpsiTruthPz", &JpsiTruthPz,16000);
    tree->Branch("JpsiTruthEnergy", &JpsiTruthE,16000);
    tree->Branch("DsTruthPx", &DsTruthPx,16000);
    tree->Branch("DsTruthPy", &DsTruthPy,16000);
    tree->Branch("DsTruthPz", &DsTruthPz,16000);
    tree->Branch("DsTruthEnergy", &DsTruthE,16000);
    
   



  return StatusCode::SUCCESS;
}



StatusCode WJDNtupAlg::finalize() {
  ATH_MSG_INFO ("Finalizing " << name() << "...");
  //
  //Things that happen once at the end of the event loop go here
  //

  return StatusCode::SUCCESS;
}



StatusCode WJDNtupAlg::execute() {  
  ATH_MSG_DEBUG ("Executing " << name() << "...");
  setFilterPassed(false); //optional: start with algorithm not passed

//-------------------------------------------------------Event Counter-----------------------------------------------------------

  evt_count++;
  std::cout << "===== " << evt_count << " event =====" << std::endl;
//  if (evt_count % 100 == 0) {std::cout << evt_count << " events processed" <<std::endl;}

  const xAOD::EventInfo* eventInfo = 0;
  CHECK(evtStore()->retrieve( eventInfo, "EventInfo"));  
  runNum = eventInfo->runNumber();
  evtNum = eventInfo->eventNumber();
  lumiNum  = eventInfo->lumiBlock();

 



//-------------------------------------------------------Trigger List-------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
   auto allChains = m_tdt->getChainGroup(".*");
   //std::cout << "TriggerNames is passed : ";
   for(auto& trig : allChains->getListOfTriggers()) 
   {
   auto cg = m_tdt->getChainGroup(trig);
   std::string thisTrig = trig;
   TrigNames->push_back(thisTrig);
   if (cg->isPassed()) {matchTrigNames->push_back((thisTrig));}
   TrigRes->push_back(cg->getPrescale());
      //ANA_MSG_INFO ("execute(): " << thisTrig << ", chain passed(1)/failed(0) = " << cg->isPassed() << ", total chain prescale (L1*HLT) = " << cg->getPrescale());
    }

  
   auto L1Chains = m_tdt->getChainGroup("L1_.*");
   for(auto& trig1 : L1Chains->getListOfTriggers()) 
   { 
   auto cg1 = m_tdt->getChainGroup(trig1);
   L1TrigRes->push_back(cg1->getPrescale());
   }
  
   
   



//-------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------Truth Particles----------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
if(doTruth) {
  // check if the event is MC
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
  isMC = true; // can do something with this later
  }
}



if(isMC && doTruth)
{  
  std::cout<<"==================DO TRUTH=========================="<<std::endl; //start truth
  nMC = eventInfo->mcEventNumber(); 

   const xAOD::TruthVertexContainer * truevertices = 0;
   CHECK(evtStore()->retrieve( truevertices, "TruthVertices" ));
  
   xAOD::TruthVertexContainer::const_iterator  truevertices_itr =  truevertices->begin();
   xAOD::TruthVertexContainer::const_iterator  truevertices_end =  truevertices->end();
   
    for( ; truevertices_itr != truevertices_end; ++truevertices_itr ) 
    {
     MCVx->push_back((*truevertices_itr)->x());
     MCVy->push_back((*truevertices_itr)->y());
     MCVz->push_back((*truevertices_itr)->z());
     }


   const xAOD::TruthParticleContainer* particles = 0;
   CHECK(evtStore()->retrieve( particles, "TruthParticles" ));
   
   for (const auto& truth: *particles) 
     {
      MCPdgIdAll->push_back(truth->pdgId());
      MCEtaAll->push_back(truth->eta());
      MCPhiAll->push_back(truth->phi());
      }
   
  
   xAOD::TruthParticleContainer::const_iterator particle_itr = particles->begin();
   xAOD::TruthParticleContainer::const_iterator particle_end = particles->end();
   
    for( ; particle_itr != particle_end; ++particle_itr ) 
    { //start loop
    	 itr++;
       if ((*particle_itr)->pt()*0.001 > T_pTmax) T_pTmax = (*particle_itr)->pt()*0.001;
       if ((*particle_itr)->pt()*0.001 < T_pTmin && (*particle_itr)->pt() != 0) T_pTmin = (*particle_itr)->pt()*0.001;
       if ((std::abs((*particle_itr)->pdgId()) == 24) && ((*particle_itr)->pt() != 0))  
       //if (std::abs((*particle_itr)->pdgId()) == 333)  {std::cout << " Phi found " << std::endl;}

       if (std::abs((*particle_itr)->pdgId()) == 24) {
          for(i = 0; i < (*particle_itr)->nChildren(); i++) {
         
          if ((*particle_itr)->child(i)->pdgId() == 443) 
          	  {T_Jpsi = true; T_nJpsi = i;}
          if (std::abs((*particle_itr)->child(i)->pdgId()) == 431) 
          	 { T_Ds = true; T_nDs = i;}    
          }
          if (T_Jpsi && T_Ds) T_WJD = true;
          T_Jpsi = false; T_Ds = false;
         }

  
        if (T_WJD)
        {
          T_WJD = false;
          //std::cout << " W found " << std::endl;

          
          if (Verbose) std::cout << " TRUTH: W -> Jpsi + Ds found" << std::endl;


          if (std::abs((*particle_itr)->child(T_nJpsi)->child(0)->pdgId()) != 13 
          && std::abs((*particle_itr)->child(T_nJpsi)->child(1)->pdgId()) != 13) continue;

          if (Verbose) std::cout << " TRUTH: Jpsi -> mu + mu found" << std::endl;

          if (((*particle_itr)->child(T_nDs)->nChildren() != 3) && ((*particle_itr)->child(T_nDs)->nChildren() != 2)) continue;
          if (4) std::cout << " TRUTH: Ds has " << (*particle_itr)->child(T_nDs)->nChildren() << " daughters " << std::endl;

           //case Ds -> KKPi
           T_nKaons = 0; j = 0;
             if ((*particle_itr)->child(T_nDs)->nChildren() == 3) {
               for (i = 0; i < 3; i++) {
                 if (std::abs((*particle_itr)->child(T_nDs)->child(i)->pdgId()) == 321) {T_nKaons++;} 
                 else  {j = i;}
                }
              if((T_nKaons == 2) && (std::abs((*particle_itr)->child(T_nDs)->child(j)->pdgId()) == 211)) {T_nKaons = 0; j = 0;}
              if (Verbose) std::cout << " TRUTH: Ds -> Pi + K + K found " << std::endl;
              else T_DsPass = false;
              }
            
            //case Ds -> PhiPi
            //here T_nKaons actually counts pions
            if ((*particle_itr)->child(0)->nChildren() == 2) {
              for (i = 0; i < 2; i++) {
                if (std::abs((*particle_itr)->child(T_nDs)->child(i)->pdgId()) == 211) {T_nKaons++; T_nPi = i;} 
                else  {j = i;}
               }
            if((T_nKaons == 1) && (std::abs((*particle_itr)->child(T_nDs)->child(j)->pdgId()) == 333)) {
               T_nKaons = 0; j = 0; T_nPhi = j;
              if (Verbose) std::cout << " TRUTH: Ds -> Pi + Phi found " << std::endl;
             }
            else T_DsPass = false;
             }


          if (T_DsPass) continue;

          TvW = (*particle_itr)->p4(); 

           if ((*particle_itr)->child(T_nDs)->nChildren() == 3) {
            // Ds part
            Tvds = (*particle_itr)->child(T_nDs)->p4();
            Tvds1 = (*particle_itr)->child(T_nDs)->child(0)->p4();
            Tvds2 = (*particle_itr)->child(T_nDs)->child(1)->p4();
            Tvds3 = (*particle_itr)->child(T_nDs)->child(2)->p4();           
            
      //K1Vec.DeltaR(K2Vec)
      
            DsdR = Tvds1.DeltaR(Tvds2);
      
            if (DsdR < Tvds1.DeltaR(Tvds3)) DsdR = Tvds1.DeltaR(Tvds3);
            if (DsdR < Tvds2.DeltaR(Tvds3)) DsdR = Tvds2.DeltaR(Tvds3);
            
            Tvseed = Tvds1;
            if (Tvseed.Pt() < Tvds2.Pt()) Tvseed = Tvds2;
            if (Tvseed.Pt() < Tvds3.Pt()) Tvseed = Tvds3;
          
            }



// for pi phi
      if ((*particle_itr)->child(T_nDs)->nChildren() == 2) {
            // Ds part
            Tvds = (*particle_itr)->child(T_nDs)->p4();
            Tvds1 = (*particle_itr)->child(T_nDs)->child(0)->p4();
            Tvds2 = (*particle_itr)->child(T_nDs)->child(1)->p4();

            TvPhiK1 = (*particle_itr)->child(T_nDs)->child(T_nPhi)->child(0)->p4();
            TvPhiK2 = (*particle_itr)->child(T_nDs)->child(T_nPhi)->child(1)->p4();

      TvPhi = TvPhiK1 + TvPhiK2;

        DsTruthPx->push_back(Tvds.Px());
        DsTruthPy->push_back(Tvds.Py());
        DsTruthPz->push_back(Tvds.Pz());
        DsTruthE->push_back(Tvds.E());
        
        
      }

      //start of Jpsi part
      Tvj = (*particle_itr)->child(T_nJpsi)->p4();
      Tvmu1 = (*particle_itr)->child(T_nJpsi)->child(0)->p4();
      Tvmu2 = (*particle_itr)->child(T_nJpsi)->child(1)->p4();
      
      
        JpsiTruthPx->push_back(Tvj.Px());
        JpsiTruthPy->push_back(Tvj.Py());
        JpsiTruthPz->push_back(Tvj.Pz());
        JpsiTruthE->push_back(Tvj.E());
          

          MCPdgId->push_back((*particle_itr)->pdgId());
          //MCP->push_back((*particle_itr)->p4());
          MCPx->push_back((*particle_itr)->p4().Px());
          MCPy->push_back((*particle_itr)->p4().Py());
          MCPz->push_back((*particle_itr)->p4().Pz());
          MCPt->push_back((*particle_itr)->pt());
          MCE->push_back((*particle_itr)->p4().E());
          MCPhi->push_back((*particle_itr)->phi());
          MCeta->push_back((*particle_itr)->eta());
          float temp = (*particle_itr)->eta();
          MCtheta->push_back(2*atan(1/exp(temp)));
          
          MCCharge->push_back((*particle_itr)->charge());
          MCStatus->push_back((*particle_itr)->status());
          MCNParent->push_back((*particle_itr)->nParents());
          MCNDaughters->push_back((*particle_itr)->nChildren()); 

      
          TvW = Tvj+Tvds;
          //end Jpsi
      
   }     
        


}//stop loop over truth particles



  T_pTmax = 0;
  T_pTmin = 65535;
std::cout<<"==================END TRUTH=========================="<<std::endl; 
}//stop truth




//----------------------------------------------Primaries And General Muons-----------------------------------------------------
if (doData) {

std::cout<<"==================DO TRACK=========================="<<std::endl; 
  
const xAOD::VertexContainer* pvtx;
CHECK (evtStore()->retrieve (pvtx, "PrimaryVertices"));

xAOD::VertexContainer::const_iterator pvtx_itr = pvtx->begin();
xAOD::VertexContainer::const_iterator pvtx_end = pvtx->end();

PrimX = (*pvtx_itr)->x();
PrimY = (*pvtx_itr)->y();
PrimZ = (*pvtx_itr)->z();

Trk::Vertex* vtx = new Trk::Vertex((*pvtx_itr)->position());
xAOD::Vertex* jvtx = new xAOD::Vertex();


for( ; pvtx_itr != pvtx_end; ++pvtx_itr ) 
{
     allpriVtxX->push_back((*pvtx_itr)->position().x());
     allpriVtxY->push_back((*pvtx_itr)->position().y());
     allpriVtxZ->push_back((*pvtx_itr)->position().z());
     allpriVtxChi->push_back((*pvtx_itr)->chiSquared());
     allpriVtxChiNorm->push_back(((*pvtx_itr)->chiSquared())/((*pvtx_itr)->numberDoF()));
     allpriVtxXE->push_back((*pvtx_itr)->covariance()[0]);
     allpriVtxYE->push_back((*pvtx_itr)->covariance()[2]);
     allpriVtxZE->push_back((*pvtx_itr)->covariance()[5]);
    
} // end for loop over primaries

const xAOD::MuonContainer* muons = 0;
CHECK (evtStore()->retrieve (muons, "Muons"));

  // loop over the muons in the container
for (auto muon : *muons) {
//    ANA_MSG_INFO ("execute(): original muon pt = " << ((muon)->pt() * 0.001) << " GeV"); // just to print out something
muon_count++;
//hist("h_muonpt")->Fill ((muon)->pt() * 0.001); // GeV
 xAOD::MuonContainer::const_iterator mu_tracks_itr1 = muons->begin();
 xAOD::MuonContainer::const_iterator mu_tracks_end = muons->end();
 for( ; mu_tracks_itr1 != mu_tracks_end; ++mu_tracks_itr1 ) { allmuQual->push_back((*mu_tracks_itr1)->quality());}
  } // end for loop over muons
//ANA_MSG_INFO ("execute(): muons this event = " << muon_count); // just to print out something

nMu = muon_count;
muon_count=0;


//-------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------SEARCH FOR J/psi----------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------

  std::vector<const xAOD::TrackParticle*> CMuonTracks; CMuonTracks.clear(); //Muon tracks storage
  

  const xAOD::TrackParticleContainer* cm_tracks;
  CHECK (evtStore()->retrieve (cm_tracks, "CombinedMuonTrackParticles"));

 
  xAOD::TrackParticleContainer::const_iterator cm_tracks_itr1 = cm_tracks->begin();
  xAOD::TrackParticleContainer::const_iterator cm_tracks_end = cm_tracks->end();
  cm_count = 0;
  std::cout << "combined muon tracks found: " << cm_tracks->size() << std::endl;
  Jpsi_Candidates = 0;
  //loop over the combined muon tracks
  for( ; cm_tracks_itr1 != cm_tracks_end; ++cm_tracks_itr1 ) 
  {
    cm_count++;
          
          allmuPx->push_back(((*cm_tracks_itr1)->p4()).Px());
          allmuPy->push_back(((*cm_tracks_itr1)->p4()).Py());
          allmuPz->push_back(((*cm_tracks_itr1)->p4()).Pz());
          allmuCharge->push_back((*cm_tracks_itr1)->charge());
          allmuType->push_back((*cm_tracks_itr1)->type());
          allmuD0->push_back((*cm_tracks_itr1)->d0());
          allmuDz->push_back((*cm_tracks_itr1)->z0());
          allmuChi2->push_back((*cm_tracks_itr1)->chiSquared());
          allmuNDF->push_back((*cm_tracks_itr1)->numberDoF());
          allmufHits->push_back((*cm_tracks_itr1)->identifierOfFirstHit());
          //muDzVtx->push_back(std::abs((*cm_tracks_itr1)->vertex()->z() - PrimZ));

         
    if (cm_tracks->size() < 2) continue; //less than two tracks - skipping
    
       CMuonTracks.push_back((*cm_tracks_itr1)); //Filling the Muon tracks vector
   
    if (cm_tracks->size() >= 2)
    { //size check
      for( xAOD::TrackParticleContainer::const_iterator cm_tracks_itr2 = cm_tracks_itr1+1 ; cm_tracks_itr2 != cm_tracks_end; ++cm_tracks_itr2 ) {
        if (cm_tracks_itr2 == cm_tracks_itr1) {std::cout << " something doesn't work" <<std::endl;}
        if ((*cm_tracks_itr1)->charge() == (*cm_tracks_itr2)->charge()) {/*std::cout << evt_count << " - event, same charge, skipping" <<std::endl;*/continue;}
        //std::cout << evt_count << " - event, we got a Jpsi candidate" <<std::endl;
        //here's the part where we fit the vertices
        CMuonTracks.push_back((*cm_tracks_itr2)); //2nd track for the muon track vector
        

        const Trk::Perigee& aPerigee2 = (*cm_tracks_itr1)->perigeeParameters();
        const Trk::Perigee& aPerigee1 = (*cm_tracks_itr2)->perigeeParameters();
        int sflag = 0;
        int errorcode = 0;
        Amg::Vector3D startingPoint = m_vertexEstimator->getCirclesIntersectionPoint(&aPerigee1,&aPerigee2,sflag,errorcode);
        if (errorcode != 0) {startingPoint(0) = 0.0; startingPoint(1) = 0.0; startingPoint(2) = 0.0;}
        
        xAOD::Vertex* myVxCandidate = m_iVertexFitter->fit(CMuonTracks, startingPoint); //fitting the two muon tracks into a Jpsi
        //if(!(myVxCandidate != 0)) {std::cout << "vertex cannot be fit" <<std::endl;}
        if(myVxCandidate != 0){
          std::vector<ElementLink<DataVector<xAOD::TrackParticle> > > newLinkVector;
          for(unsigned int i=0; i< myVxCandidate->trackParticleLinks().size(); i++)
            { ElementLink<DataVector<xAOD::TrackParticle> > mylink=myVxCandidate->trackParticleLinks()[i]; //makes a copy (non-const) 
            mylink.setStorableObject(*cm_tracks, true); 
            newLinkVector.push_back( mylink ); }
          myVxCandidate->clearTracks();
          myVxCandidate->setTrackParticleLinks( newLinkVector );
          jvtx = myVxCandidate;


          //need to look if there's a differance between Lz vec sum and summing the linked tracks
          vmu1 = (*cm_tracks_itr1)->p4();
          vmu2 = (*cm_tracks_itr2)->p4();

          vj = vmu1 + vmu2;
          
          
          if (((*cm_tracks_itr1)->charge() > 0) && ((*cm_tracks_itr2)->charge() < 0)) {
             mupPx->push_back((*cm_tracks_itr1)->p4().Px());
             mupPy->push_back((*cm_tracks_itr1)->p4().Py());
             mupPz->push_back((*cm_tracks_itr1)->p4().Pz());
             mupfChi2->push_back((*cm_tracks_itr1)->chiSquared());
             mupfNDF->push_back((*cm_tracks_itr1)->numberDoF());
             mumPx->push_back((*cm_tracks_itr2)->p4().Px());
             mumPy->push_back((*cm_tracks_itr2)->p4().Py());
             mumPz->push_back((*cm_tracks_itr2)->p4().Pz());
             mumfChi2->push_back((*cm_tracks_itr2)->chiSquared());
             mumfNDF->push_back((*cm_tracks_itr2)->numberDoF());
               
             }
               else {
                    mupPx->push_back((*cm_tracks_itr2)->p4().Px());
                    mupPy->push_back((*cm_tracks_itr2)->p4().Py());
                    mupPz->push_back((*cm_tracks_itr2)->p4().Pz());
                    mupfChi2->push_back((*cm_tracks_itr2)->chiSquared());
                    mupfNDF->push_back((*cm_tracks_itr2)->numberDoF());
                    mumPx->push_back((*cm_tracks_itr1)->p4().Px());
                    mumPy->push_back((*cm_tracks_itr1)->p4().Py());
                    mumPz->push_back((*cm_tracks_itr1)->p4().Pz());
                    mumfChi2->push_back((*cm_tracks_itr1)->chiSquared());
                    mumfNDF->push_back((*cm_tracks_itr1)->numberDoF());
                    
                    }

   
          


          //for leading mu > 20 GeV pT
          //if ( (vmu1.Pt()*0.001 < 20) && (vmu2.Pt()*0.001 < 20)) continue;
          
          //data mumu pt compare plot
          


          if (dofit) {/*m_CascadeFitter->cleanCascade();*/ m_CascadeFitter->nextVertex(CMuonTracks,JpsiMassVector,3096.916 );}


          if (vj.Pt()*0.001 > JpsiPtCut) {
            //Which primary is the closest?
            for( ; pvtx_itr != pvtx_end; ++pvtx_itr ) {
              dPrim_Jpsi_temp = std::sqrt( std::pow((*pvtx_itr)->x()-myVxCandidate->x(),2) + std::pow((*pvtx_itr)->y()-myVxCandidate->y(),2) + std::pow((*pvtx_itr)->z()-myVxCandidate->z(),2));
              if (dPrim_Jpsi > dPrim_Jpsi_temp) {
                PrimX = (*pvtx_itr)->x(); PrimY = (*pvtx_itr)->y(); PrimZ = (*pvtx_itr)->z(); priVtxChi = (*pvtx_itr)->chiSquared();  priVtxChiNorm = ((*pvtx_itr)->chiSquared()/(*pvtx_itr)->numberDoF());
                //PrimXE=((*pvtx_itr)->covariance()[0]);
      	        //PrimYE=((*pvtx_itr)->covariance()[2]);
      	        //PrimZE=((*pvtx_itr)->covariance()[5]);
                dPrim_Jpsi = dPrim_Jpsi_temp;
                *vtx = Trk::Vertex((*pvtx_itr)->position());
              }

               priVtxX = PrimX;
               priVtxY = PrimY;
               priVtxZ = PrimZ;
               priVtxXE = PrimXE;
               priVtxYE = PrimYE;
               priVtxZE = PrimZE;
               
             
            } // end for loop over primaries


            Jpsi_z0 = myVxCandidate->z();
          
                JpsiTrackPx->push_back(vj.Px());
                JpsiTrackPy->push_back(vj.Py());
                JpsiTrackPz->push_back(vj.Pz());
                JpsiTrackE->push_back(vj.Energy());
                JMass->push_back(vj.M());
                JDecayVtxX->push_back(myVxCandidate->x());
                JDecayVtxY->push_back(myVxCandidate->y());
                JDecayVtxZ->push_back(myVxCandidate->z());

                JVtxC2->push_back(myVxCandidate->chiSquared());
                
                JDecayVtxXE->push_back(myVxCandidate->covariance()[0]);
      	        JDecayVtxYE->push_back(myVxCandidate->covariance()[2]);
      	        JDecayVtxZE->push_back(myVxCandidate->covariance()[5]);
                
      
          }
          //xAOD::CompositeParticle* Reco_Jpsi;
          //Reco_Jpsi.setP4(vj);
          //Reco_Jpsi.addPart
          Jpsi_Candidates++;
          nJ++;
         
        
        }

         CMuonTracks.pop_back(); //remove the last track (it should be the 2nd muon)
       
      } //end loop over 2nd muons 
          
    }// end of size check    
   // reset the vector for a new 1st muon, if everything is right this should be == CMuonTracks.pop_back();
           
    } //end loop over 1st muons

    nMu=cm_count++;

  std::cout << Jpsi_Candidates << " Jpsi candidates found" << std::endl;

       
       
//-------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------SEARCH FOR Ds-----------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------

Ds_Candidates = 0;
  std::vector<const xAOD::TrackParticle*> InDetTracks; InDetTracks.clear(); //Muon tracks storage
  //std::vector<const xAOD::TrackParticle*> CMuonTrackPair; CMuonTrackPair.clear(); //Muon tracks storage

  const xAOD::TrackParticleContainer* indet_tracks;
  CHECK (evtStore()->retrieve (indet_tracks, "InDetTrackParticles"));

  if(indet_tracks->size() != 0) {trkSize=(indet_tracks->size());}
  xAOD::TrackParticleContainer::const_iterator indet_tracks_itr1 = indet_tracks->begin();
  xAOD::TrackParticleContainer::const_iterator indet_tracks_end = indet_tracks->end();
  indet_count = 0;
  indet_primcut = 0;
  if (Jpsi_Candidates != 0) std::cout << "inner detector tracks found:" << indet_tracks->size() << std::endl;

  //loop over the combined indet tracks
  for( ; indet_tracks_itr1 != indet_tracks_end; ++indet_tracks_itr1 ) 
  { // -1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1--1-
    
    
   // trToJpsiTagDZ->push_back((std::abs(((*indet_tracks_itr1)->vertex()->z()) - PrimZ)));
        
    if ( !m_selTool->accept((*indet_tracks_itr1)) ) continue;
    if (Jpsi_Candidates == 0) continue;

    //range cut out of seed track cut ---!!!!!!!!!!!!!!!!!!!!!!!!---!!!!!!!!!!!!!!!!!!!!!!!---!!!!!!!!!!!!!!!!!!!!!!---
    if (std::abs((*indet_tracks_itr1)->z0() - Jpsi_z0) <= Jpsi_range) indet_primcut++;// std::cout << indet_primcut;
    
    if (indet_tracks->size() < 3) continue; //less than three tracks - skipping
    if (SeedTrackCut) { // ===================================SEED TRACK CUT=========================================
      //std::cout << " Seed track test ..." << std::endl;
      vseed = (*indet_tracks_itr1)->p4();
      if (std::abs((*indet_tracks_itr1)->z0() - Jpsi_z0) > Jpsi_range) continue;


    if ( !((*indet_tracks_itr1)->pt() > SeedpTMax || (*indet_tracks_itr1)->pt() < SeedpTMin)) continue;
    if (vseed.DeltaR(vj) < SeedJpsidRMin /*|| vseed.Angle(vj.Vect()) > 1.7*/) continue;

    }
                 
          
          TrackPx->push_back((*indet_tracks_itr1)->p4().Px());
          TrackPy->push_back((*indet_tracks_itr1)->p4().Py());
          TrackPz->push_back((*indet_tracks_itr1)->p4().Pz());
          TrackEnergy->push_back((*indet_tracks_itr1)->p4().E());
          TrackNDF->push_back((*indet_tracks_itr1)->numberDoF());
          TrackChi2->push_back((*indet_tracks_itr1)->chiSquared());
          TrackCharge->push_back((*indet_tracks_itr1)->charge());
          TrackPt->push_back((*indet_tracks_itr1)->pt());
          TrackPhi->push_back((*indet_tracks_itr1)->phi());
          TrackEta->push_back((*indet_tracks_itr1)->eta());
          TrackD0->push_back((*indet_tracks_itr1)->d0());
          trfHits->push_back((*indet_tracks_itr1)->identifierOfFirstHit()); 
          //trDzVtx->push_back(std::abs(((*indet_tracks_itr1)->vertex()->z())-PrimZ));
   

    //std::cout << " Seed track test passed" << std::endl;
    indet_count++;

    InDetTracks.push_back((*indet_tracks_itr1)); //Filling the indet tracks vector with the 1st particle

    xAOD::TrackParticleContainer::const_iterator indet_tracks_itr2 = indet_tracks_itr1+1;
    
    if(SeedTrackCut) {xAOD::TrackParticleContainer::const_iterator indet_tracks_itr2 = indet_tracks->begin();}

    for( ; indet_tracks_itr2 != indet_tracks_end; ++indet_tracks_itr2 ) 
    { //2nd particle -2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2--2-
      if ( !m_selTool->accept((*indet_tracks_itr2))) continue;

      if(SeedTrackCut) {
        if (std::abs((*indet_tracks_itr2)->z0() - Jpsi_z0) > Jpsi_range) continue;
        if (indet_tracks_itr2 == indet_tracks_itr1) continue;
        vseed_sub = (*indet_tracks_itr2)->p4();
        if (vseed.DeltaR(vseed_sub) > dRSeedCone) continue;
          }

      //!!!!! The "continues" must be BEFORE push_back's, otherwise, a pop_back will be missed and there will be more than 3 trackpartilces in the InDetTracks vector !!!!!
      InDetTracks.push_back((*indet_tracks_itr2)); //Filling the indet tracks vector with the 2nd particle


      const Trk::Perigee& aPerigeePrefit1 = (*indet_tracks_itr1)->perigeeParameters();
      const Trk::Perigee& aPerigeePrefit2 = (*indet_tracks_itr2)->perigeeParameters();
      int sflag_Ds = 0;
      int errorcode_Ds = 0;
      Amg::Vector3D sP0 = m_vertexEstimator->getCirclesIntersectionPoint(&aPerigeePrefit1,&aPerigeePrefit2,sflag_Ds,errorcode_Ds);
      if (errorcode_Ds != 0) {sP0(0) = 0.0; sP0(1) = 0.0; sP0(2) = 0.0;}
      xAOD::Vertex* PrefitVx = m_JpsiFitter->fit(InDetTracks, sP0);
      if (PrefitVx != 0) {chi2 = PrefitVx->chiSquared(); delete PrefitVx;}
      if (chi2 > 10) {InDetTracks.pop_back(); continue;}


      //if (PrefitVx->chiSquared() > 10) {InDetTracks.pop_back(); delete PrefitVx; continue;}

      //std::cout << "not deleted" << std::endl;
      xAOD::TrackParticleContainer::const_iterator indet_tracks_itr3 = indet_tracks_itr2+1;
      for( ; indet_tracks_itr3 != indet_tracks_end; ++indet_tracks_itr3 ) 
      { //3d particle -3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3--3-
        //if ((*indet_tracks_itr3)->pt() < 200) continue; //0.2 GeV pT cut
        if ( !m_selTool->accept((*indet_tracks_itr3))) continue;
        if (((*indet_tracks_itr1)->charge() == (*indet_tracks_itr2)->charge()) && ((*indet_tracks_itr1)->charge() == (*indet_tracks_itr3)->charge())) continue; //if +++ or --- - skip
        //std::cout << (*indet_tracks_itr1)->charge() << (*indet_tracks_itr2)->charge() << (*indet_tracks_itr3)->charge() << std::endl;
        //assuming it's not a all equal charge situation
        pre_vds = (*indet_tracks_itr1)->p4() + (*indet_tracks_itr2)->p4() + (*indet_tracks_itr3)->p4();
        
        if(SeedTrackCut) {
          vseed_sub = (*indet_tracks_itr3)->p4();
          if (indet_tracks_itr3 == indet_tracks_itr1) continue;
          if (indet_tracks_itr3 == indet_tracks_itr2) continue;
          if (std::abs((*indet_tracks_itr3)->z0() - Jpsi_z0) > Jpsi_range) continue;
          if (vseed.DeltaR(vseed_sub) > dRSeedCone) continue;

           }

        

        if (pre_vds.Pt()*0.001 < DsPtCut) continue;
        //!!!!! The "continues" must be BEFORE push_back's, otherwiser there will be more than 3 trackpartilces in the InDetTracks vector !!!!! 
        InDetTracks.push_back((*indet_tracks_itr3)); //Filling the indet tracks vector with the 3d particle
        
        
        if (InDetTracks[0]->charge() == InDetTracks[1]->charge()) { k = 2; d1 = 0; d2 = 1;
        } else if (InDetTracks[0]->charge() == InDetTracks[2]->charge()) { k = 1; d1 = 0; d2 = 2;
        } else { k = 0; d1 = 1; d2 = 2;}



      
        //it seems by default all particle masses are Pi masses
        //also there has to be a smarter way to do the mass theory;
        

        DsTracks.push_back(InDetTracks[d1]);
        DsTracks.push_back(InDetTracks[k]);
        DsTracks.push_back(InDetTracks[d2]);

  

          const Trk::Perigee& aPerigeeDs1 = (*indet_tracks_itr1)->perigeeParameters();
          const Trk::Perigee& aPerigeeDs2 = (*indet_tracks_itr2)->perigeeParameters();
          const Trk::Perigee& aPerigeeDs3 = (*indet_tracks_itr3)->perigeeParameters();

          sflag_Ds = 0;
          errorcode_Ds = 0;

          Amg::Vector3D sP1 = m_vertexEstimator->getCirclesIntersectionPoint(&aPerigeeDs1,&aPerigeeDs2,sflag_Ds,errorcode_Ds);
          Amg::Vector3D sP2 = m_vertexEstimator->getCirclesIntersectionPoint(&aPerigeeDs2,&aPerigeeDs3,sflag_Ds,errorcode_Ds);
          Amg::Vector3D sP3 = m_vertexEstimator->getCirclesIntersectionPoint(&aPerigeeDs1,&aPerigeeDs3,sflag_Ds,errorcode_Ds);
          //Average???
          Amg::Vector3D startingPointDs( (sP1(0)+sP2(0)+sP3(0))/3,(sP1(1)+sP2(1)+sP3(1))/3,(sP1(2)+sP2(2)+sP3(2))/3 );
          if (errorcode_Ds != 0) {startingPointDs(0) = 0.0; startingPointDs(1) = 0.0; startingPointDs(2) = 0.0;}

          xAOD::Vertex* myDsVxCandidate = m_DsFitter->fit(DsTracks, startingPointDs);

          if(myDsVxCandidate != 0)
           {         
             DVtxC2->push_back(myDsVxCandidate->chiSquared());
             DDecayVtxX->push_back(myDsVxCandidate->x());
             DDecayVtxY->push_back(myDsVxCandidate->y());
             DDecayVtxZ->push_back(myDsVxCandidate->z());
   
                DDecayVtxXE->push_back(myDsVxCandidate->covariance()[0]);
      	        DDecayVtxYE->push_back(myDsVxCandidate->covariance()[2]);
      	        DDecayVtxZE->push_back(myDsVxCandidate->covariance()[5]);
                
             
             delete myDsVxCandidate;
             }
              DsTracks.clear();


         //==================FIRST==============================================================================================       
        vk1.SetVectM(InDetTracks[k]->p4().Vect(), 493);
        vd11.SetVectM(InDetTracks[d1]->p4().Vect(), 139);
        vd21.SetVectM(InDetTracks[d2]->p4().Vect(), 493);
        
        vds1 = vk1 + vd11 + vd21;
        vphi1 = vk1 + vd21;

             nDs1++;
             Ds_Candidates++;
             DsTrackPx1->push_back(vds1.Px());
             DsTrackPy1->push_back(vds1.Py());
             DsTrackPz1->push_back(vds1.Pz());
             DsTrackE1->push_back(vds1.Energy());
             DMass1->push_back(vds1.M());  


       if (InDetTracks[k]->charge() > 0)
         {
         KplusDsPx1->push_back(InDetTracks[k]->p4().Px());
         KplusDsPy1->push_back(InDetTracks[k]->p4().Py());
         KplusDsPz1->push_back(InDetTracks[k]->p4().Pz());
         KplusDsD01->push_back(InDetTracks[k]->d0());
         KplusDsE1->push_back(InDetTracks[k]->p4().E());
         coneJpsiKplus1->push_back(vk1.Angle(vj.Vect()));
         }   
          else {
          	KminusDsPx1->push_back(InDetTracks[k]->p4().Px());
         	KminusDsPy1->push_back(InDetTracks[k]->p4().Py());
         	KminusDsPz1->push_back(InDetTracks[k]->p4().Pz());
         	KminusDsD01->push_back(InDetTracks[k]->d0());
         	KminusDsE1->push_back(InDetTracks[k]->p4().E());
         	coneJpsiKminus1->push_back(vk1.Angle(vj.Vect()));
            }
          if (InDetTracks[d2]->charge()> 0)
          	{
             KplusDsPx1->push_back(InDetTracks[d2]->p4().Px());
             KplusDsPy1->push_back(InDetTracks[d2]->p4().Py());
             KplusDsPz1->push_back(InDetTracks[d2]->p4().Pz());
             KplusDsD01->push_back(InDetTracks[d2]->d0());
             KplusDsE1->push_back(InDetTracks[d2]->p4().E());
             coneJpsiKplus1->push_back(vd21.Angle(vj.Vect()));
              }   
              else {
          	    KminusDsPx1->push_back(InDetTracks[d2]->p4().Px());
         	    KminusDsPy1->push_back(InDetTracks[d2]->p4().Py());
         	    KminusDsPz1->push_back(InDetTracks[d2]->p4().Pz());
         	    KminusDsD01->push_back(InDetTracks[d2]->d0());
         	    KminusDsE1->push_back(InDetTracks[d2]->p4().E());
         	    coneJpsiKminus1->push_back(vd21.Angle(vj.Vect()));
                 }

                PiDsPx1->push_back(InDetTracks[d1]->p4().Px());
                PiDsPy1->push_back(InDetTracks[d1]->p4().Py());
                PiDsPz1->push_back(InDetTracks[d1]->p4().Pz());
                PiDsD01->push_back(InDetTracks[d1]->d0());
                PiDsE1->push_back(InDetTracks[d1]->p4().E());
                coneJpsiPi1->push_back(vd11.Angle(vj.Vect()));

        Ds_mass1 = vds1.M();

         if ((std::abs(Ds_mass1-1968.47)<100))
               {
                 vw1 = vds1+vj;
                 nW1++;           
                 WMassLor1->push_back(vw1.M());
                 WPxLor1->push_back(vw1.Px());
                 WPyLor1->push_back(vw1.Py());
                 WPzLor1->push_back(vw1.Pz());
                
                 }
          


//==================SECOND==============================================================================================        
        vk2.SetVectM(InDetTracks[k]->p4().Vect(), 493);
        vd12.SetVectM(InDetTracks[d1]->p4().Vect(), 493);
        vd22.SetVectM(InDetTracks[d2]->p4().Vect(), 139);
        
        vds2 = vk2 + vd12 + vd22;
        vphi2 = vk2 + vd12;
        Ds_mass2 = vds2.M();

             nDs2++;

             DsTrackPx2->push_back(vds2.Px());
             DsTrackPy2->push_back(vds2.Py());
             DsTrackPz2->push_back(vds2.Pz());
             DsTrackE2->push_back(vds2.Energy());
             DMass2->push_back(vds2.M());  

       if (InDetTracks[k]->charge() > 0)
         {
         KplusDsPx2->push_back(InDetTracks[k]->p4().Px());
         KplusDsPy2->push_back(InDetTracks[k]->p4().Py());
         KplusDsPz2->push_back(InDetTracks[k]->p4().Pz());
         KplusDsD02->push_back(InDetTracks[k]->d0());
         KplusDsE2->push_back(InDetTracks[k]->p4().E());
         coneJpsiKplus2->push_back(vk2.Angle(vj.Vect()));
         }   
          else {
          	KminusDsPx2->push_back(InDetTracks[k]->p4().Px());
         	KminusDsPy2->push_back(InDetTracks[k]->p4().Py());
         	KminusDsPz2->push_back(InDetTracks[k]->p4().Pz());
         	KminusDsD02->push_back(InDetTracks[k]->d0());
         	KminusDsE2->push_back(InDetTracks[k]->p4().E());
         	coneJpsiKminus2->push_back(vk2.Angle(vj.Vect()));
            }
          if (InDetTracks[d1]->charge()> 0)
          	{
             KplusDsPx2->push_back(InDetTracks[d1]->p4().Px());
             KplusDsPy2->push_back(InDetTracks[d1]->p4().Py());
             KplusDsPz2->push_back(InDetTracks[d1]->p4().Pz());
             KplusDsD02->push_back(InDetTracks[d1]->d0());
             KplusDsE2->push_back(InDetTracks[d1]->p4().E());
             coneJpsiKplus2->push_back(vd12.Angle(vj.Vect()));
              }   
              else {
          	    KminusDsPx2->push_back(InDetTracks[d1]->p4().Px());
         	    KminusDsPy2->push_back(InDetTracks[d1]->p4().Py());
         	    KminusDsPz2->push_back(InDetTracks[d1]->p4().Pz());
         	    KminusDsD02->push_back(InDetTracks[d1]->d0());
         	    KminusDsE2->push_back(InDetTracks[d1]->p4().E());
         	    coneJpsiKminus2->push_back(vd12.Angle(vj.Vect()));
                 }

                PiDsPx2->push_back(InDetTracks[d2]->p4().Px());
                PiDsPy2->push_back(InDetTracks[d2]->p4().Py());
                PiDsPz2->push_back(InDetTracks[d2]->p4().Pz());
                PiDsD02->push_back(InDetTracks[d2]->d0());
                PiDsE2->push_back(InDetTracks[d2]->p4().E());
                coneJpsiPi2->push_back(vd22.Angle(vj.Vect()));
       
             
             if ((std::abs(Ds_mass2-1968.47)<100)){
                 vw2 = vds2+vj;
                 nW2++;
                 WMassLor2->push_back(vw2.M());
                 WPxLor2->push_back(vw2.Px());
                 WPyLor2->push_back(vw2.Py());
                 WPzLor2->push_back(vw2.Pz());   
                 }
            
      
        InDetTracks.pop_back();
      } // end of loop over 3d track
     
      InDetTracks.pop_back();
    } //end loop over 2nd track      
 
  InDetTracks.clear(); 
  CMuonTracks.clear();
 
   // reset the vector for a new 1st track, if everything is right this should be == InDetTracks.pop_back();
 } //end loop over 1st track
  
  if(indet_count != 0) {TrkNum=indet_count;}

  if (Jpsi_Candidates != 0) std::cout << Ds_Candidates << " Ds candidates found" << std::endl;
  setFilterPassed(true); //if got here, assume that means algorithm passed
  delete vtx;
  delete jvtx;
 

if ((Jpsi_Candidates != 0)&&(Ds_Candidates != 0))  {tree->Fill();}

std::cout<<"==================END TRACK=========================="<<std::endl;
} //doData END

      

          TrigNames->clear();
          TrigRes->clear();
          matchTrigNames->clear();
          L1TrigRes->clear();

          JpsiTruthPx->clear();
          JpsiTruthPy->clear();
          JpsiTruthPz->clear();
          JpsiTruthE->clear();
          DsTruthPx->clear();
          DsTruthPy->clear();
          DsTruthPz->clear();
          DsTruthE->clear();

          TrackPx->clear();
          TrackPy->clear();
          TrackPz->clear();
          TrackEnergy->clear();
          TrackNDF->clear();
          TrackChi2->clear();
          TrackCharge->clear();
          TrackPt->clear();
          TrackPhi->clear();
          TrackEta->clear();
          trfHits->clear();
          //trDzVtx->clear();
          TrackD0->clear();
          //trToJpsiTagDZ->clear();
          
          JpsiTrackPx->clear();
          JpsiTrackPy->clear();
          JpsiTrackPz->clear();
          JpsiTrackE->clear();
          JMass->clear();   
          JDecayVtxX->clear();
          JDecayVtxY->clear();
          JDecayVtxZ->clear();
          JDecayVtxXE->clear();
          JDecayVtxYE->clear();
          JDecayVtxZE->clear();
          JVtxC2->clear();
          
          
          DDecayVtxX->clear();
          DDecayVtxY->clear();
          DDecayVtxZ->clear();
          DDecayVtxXE->clear();
          DDecayVtxYE->clear();
          DDecayVtxZE->clear();
          DVtxC2->clear();
        
          DsTrackPx1->clear();
          DsTrackPy1->clear();
          DsTrackPz1->clear();
          DsTrackE1->clear();
          DMass1->clear();       

          DsTrackPx2->clear();
          DsTrackPy2->clear();
          DsTrackPz2->clear();
          DsTrackE2->clear();
          DMass2->clear();
          


          KplusDsPx1->clear();
          KplusDsPy1->clear();
          KplusDsPz1->clear();
          KplusDsD01->clear();
          KplusDsE1->clear();
          KminusDsPx1->clear();
          KminusDsPy1->clear();
          KminusDsPz1->clear();
          KminusDsD01->clear();
          KminusDsE1->clear();
          PiDsPx1->clear();
          PiDsPy1->clear();
          PiDsPz1->clear();
          PiDsD01->clear();
          PiDsE1->clear();
          coneJpsiKplus1->clear();
          coneJpsiKminus1->clear();
          coneJpsiPi1->clear();

          KplusDsPx2->clear();
          KplusDsPy2->clear();
          KplusDsPz2->clear();
          KplusDsD02->clear();
          KplusDsE2->clear();
          KminusDsPx2->clear();
          KminusDsPy2->clear();
          KminusDsPz2->clear();
          KminusDsD02->clear();
          KminusDsE2->clear();
          PiDsPx2->clear();
          PiDsPy2->clear();
          PiDsPz2->clear();
          PiDsD02->clear();
          PiDsE2->clear();
          coneJpsiKplus2->clear();
          coneJpsiKminus2->clear();
          coneJpsiPi2->clear();

          WPxLor1->clear();
          WPyLor1->clear();
          WPzLor1->clear();
          WMassLor1->clear();
          WPxLor2->clear();
          WPyLor2->clear();
          WPzLor2->clear();
          WMassLor2->clear(); 
          
          //WVtxC2->clear();
          //WDecayVtxX->clear();
          //WDecayVtxY->clear();
          //WDecayVtxZ->clear();

          
          MCPdgIdAll->clear();
          MCEtaAll->clear();
          MCPhiAll->clear();
          MCCharge->clear();
          MCPdgId->clear();
          MCPx->clear();
          MCPy->clear();
          MCPz->clear();
          MCPt->clear();
          MCE->clear();
          MCPhi->clear();
          MCeta->clear();
          MCtheta->clear();
          MCStatus->clear();
          MCVx->clear();
          MCVy->clear();
          MCVz->clear();
          MCNDaughters->clear();
          MCNParent->clear();

          
                
          allmuPx->clear();
          allmuPy->clear();
          allmuPz->clear();
          allmuCharge->clear();
          allmuType->clear();
          allmuD0->clear();
          allmuDz->clear();
          allmuChi2->clear();
          allmuNDF->clear();
          allmuQual->clear();
          allmufHits->clear();
          //muDzVtx->clear();

          mupPx->clear();
          mupPy->clear();
          mupPz->clear();
          mumPx->clear();
          mumPy->clear();
          mumPz->clear();
          mumfChi2->clear();
          mupfChi2->clear();
          mumfNDF->clear();
          mupfNDF->clear();

          allpriVtxX->clear();
          allpriVtxY->clear();
          allpriVtxZ->clear();
          allpriVtxXE->clear();
          allpriVtxYE->clear();
          allpriVtxZE->clear();
          allpriVtxChi->clear();
          allpriVtxChiNorm->clear();
          
          TrkNum=0; priVtxZE=0; priVtxYE =0; priVtxXE = 0;PrimZE=0,PrimYE=0,PrimXE=0;
                   
               

  return StatusCode::SUCCESS;
}



StatusCode WJDNtupAlg::beginInputFile() { 
  //
  //This method is called at the start of each input file, even if
  //the input file contains no events. Accumulate metadata information here
  //

  //example of retrieval of CutBookkeepers: (remember you will need to include the necessary header files and use statements in requirements file)
  // const xAOD::CutBookkeeperContainer* bks = 0;
  // CHECK( inputMetaStore()->retrieve(bks, "CutBookkeepers") );

  //example of IOVMetaData retrieval (see https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/AthAnalysisBase#How_to_access_file_metadata_in_C)
  //float beamEnergy(0); CHECK( retrieveMetadata("/TagInfo","beam_energy",beamEnergy) );
  //std::vector<float> bunchPattern; CHECK( retrieveMetadata("/Digitiation/Parameters","BeamIntensityPattern",bunchPattern) );



  return StatusCode::SUCCESS;
}
