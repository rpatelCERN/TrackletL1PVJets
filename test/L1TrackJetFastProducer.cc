//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "DataFormats/L1TVertex/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/CaloClusterer.h"
#include "L1Trigger/TrackFindingTracklet/interface/StubPtConsistency.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TLorentzVector.h>
//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

/////////////
//// Fast Jet
#include "RecoJets/JetProducers/plugins/VirtualJetProducer.h"
//////////////
// NAMESPACES
//#include "tracklet_em_disp.h"
using namespace std;
using namespace edm;
using namespace l1t;

//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TrackJetFastProducer : public EDAnalyzer
{
public:

typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
  // Constructor/destructor
  explicit L1TrackJetFastProducer(const ParameterSet& iConfig);
  virtual ~L1TrackJetFastProducer();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const Event& iEvent, const EventSetup& iSetup);
  virtual bool TrackQualityCuts(float trk_pt,int trk_nstub, double trk_chi2);
  virtual void FillFastJets(std::vector<float>pt, std::vector<float>eta, std::vector<float>phi, std::vector<float>p, std::vector<float>z0, std::vector<int> TruthID, float conesize, std::vector<fastjet::PseudoJet> &JetOutputs_, std::vector<int>&JetNtracks, std::vector<float>&JetVz, std::vector<float>&TrueSumPt);
  virtual void FillCaloJets(Handle<reco::PFJetCollection> PFCaloJetHandle,Handle< L1TkJetParticleCollection> CaloTkJetHandle);
 // virtual etaphibin * L1_cluster(etaphibin * phislice);
 // virtual void  L2_cluster(vector< Ptr< L1TTTrackType > > L1TrackPtrs, vector<int>ttrk, vector<int>tdtrk,vector<int>ttdtrk, maxzbin &mzb);
protected:
  
private:
  int Zbins=60;
  float zstep=15.*2./Zbins;  
  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  ParameterSet config; 
  
  int MyProcess;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;       // lots of debug printout statements
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  int L1Tk_nPar;        // use 4 or 5 parameter track fit? 
  int NPS_minStubs;
  int TP_minNStub;      // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer; // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;      // save TPs with pt > minPt 
  double TP_maxEta;     // save TPs with |eta| < maxEta 
  double TP_maxZ0;      // save TPs with |z0| < maxZ0 
  double CHI2MAX;       // save Chi2 per ndof <5
  double DeltaZ0Cut;    // save with |L1z-z0| < maxZ0
  double CONESize;      // Use anti-kt with this cone size
  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)
  double PTMAX; 
  InputTag L1TrackInputTag;        // L1 track collection
  InputTag MCTruthTrackInputTag;   // MC truth collection
  InputTag L1StubInputTag;
  InputTag MCTruthStubInputTag;
  InputTag TrackingParticleInputTag;
  InputTag TwoLayerTkJetInputTag; 
  InputTag TkJetInputTag; 
  InputTag CaloTkJetInputTag; 
  InputTag CaloJetInputTag; 
  InputTag PFJetInputTag; 
  InputTag RecoVertexInputTag;
  InputTag  GenParticleInputTag;
  InputTag GenJetAK4;
/* 
      vector<int> zbincount;
       vector<int> ttrk;
       vector<int> tdtrk;
       vector<int> ttdtrk;
vector< Ptr< L1TTTrackType > > L1TwoLayerInputPtrs;
*/
  EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;

  EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;

  EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  EDGetTokenT< vector<reco::GenParticle> > HEPMCVertexToken_;
  EDGetTokenT<std::vector<reco::GenJet> >GenJetCollectionToken_;
  EDGetTokenT< VertexCollection>L1VertexToken_;
  //EDGetTokenT< L1TkJetParticleCollection >L1TwoLayerTkJetsToken_;
  EDGetTokenT< L1TkJetParticleCollection >L1TkJetsToken_;
  EDGetTokenT< L1TkJetParticleCollection >L1CaloTkJetsToken_;
  EDGetTokenT<reco::PFJetCollection>L1PFCaloClusterToken_;
  EDGetTokenT<reco::PFJetCollection>L1PFJetToken_;
  std::vector<fastjet::PseudoJet>  GenJetInputs_;
  std::vector<fastjet::PseudoJet>  JetInputs_;
  std::vector<float>JetVz_;
  std::vector<int>JetNtracks_;
  std::vector<float>JetEventID_;
  std::vector<fastjet::PseudoJet> JetOutputs_;
  std::vector<fastjet::PseudoJet> fjConstituents_;
  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;
  // primary vertex
  //std::vector<float>* m_pv_L1recofakesumpt;
  std::vector<float>* m_pv_L1recotruesumpt;
  std::vector<float>* m_pv_L1recosumpt;
  std::vector<float>* m_pv_L1reco;
  std::vector<float>* m_pv_L1TP;
  std::vector<float>* m_pv_L1TPsumpt;
  std::vector<float>* m_pv_MC;
  std::vector<float>* m_pv_MCChgSumpT;
  std::vector<int>* m_MC_lep;
  // all L1 tracks

  std::vector<float>* m_trk_p;
  std::vector<float>* m_trk_pt;
  std::vector<float>* m_trk_eta;
  std::vector<float>* m_trk_phi;
  std::vector<float>* m_trk_d0;   // (filled if L1Tk_nPar==5, else 999)
  std::vector<float>* m_trk_z0;
  std::vector<float>* m_trk_chi2; 
  std::vector<int>*   m_trk_psnstub;
  std::vector<int>*   m_trk_nstub;
  std::vector<int>*   m_trk_genuine;
  std::vector<int>*   m_trk_loose;
  std::vector<int>*   m_trk_unknown;
  std::vector<int>*   m_trk_combinatoric;
  std::vector<int>*   m_trk_fake; //0 fake, 1 track from primary interaction, 2 secondary track
  std::vector<int>*   m_trk_matchtp_pdgid;
  std::vector<float>* m_trk_matchtp_pt;
  std::vector<float>* m_trk_matchtp_eta;
  std::vector<float>* m_trk_matchtp_phi;
  std::vector<float>* m_trk_matchtp_z0;
  std::vector<float>* m_trk_matchtp_dxy;
  std::vector<float>* m_trk_bconsist;
  std::vector<float>* m_trk_sconsist;
  // all tracking particles
  std::vector<float>* m_tp_p;
  std::vector<float>* m_tp_pt;
  std::vector<float>* m_tp_eta;
  std::vector<float>* m_tp_phi;
  std::vector<float>* m_tp_dxy;
  std::vector<float>* m_tp_d0;
  std::vector<float>* m_tp_z0;
  std::vector<float>* m_tp_d0_prod;
  std::vector<float>* m_tp_z0_prod;
  std::vector<int>*   m_tp_pdgid;
  std::vector<int>*   m_tp_nmatch;
  std::vector<int>*   m_tp_nstub;
  std::vector<int>*   m_tp_nstublayers;
  std::vector<int>*   m_tp_eventid;
  std::vector<int>*   m_tp_charge;
  // *L1 track* properties if m_tp_nmatch > 0
/*
  std::vector<float>* m_matchtrk_p;
  std::vector<float>* m_matchtrk_pt;
  std::vector<float>* m_matchtrk_eta;
  std::vector<float>* m_matchtrk_phi;
  std::vector<float>* m_matchtrk_d0; //this variable is only filled if L1Tk_nPar==5
  std::vector<float>* m_matchtrk_z0;
  std::vector<float>* m_matchtrk_chi2; 
  std::vector<int>*   m_matchtrk_nstub;
*/
  // ALL stubs
  std::vector<float>* m_genak4jet_phi;
  std::vector<float>* m_genak4jet_neufrac;
  std::vector<float>* m_genak4jet_chgfrac;
  std::vector<float>* m_genak4jet_metfrac;
  std::vector<float>* m_genak4jet_eta;
  std::vector<float>* m_genak4jet_pt;
  std::vector<float>* m_genak4jet_p;
  std::vector<float>* m_genak4chgjet_phi;
  std::vector<float>* m_genak4chgjet_eta;
  std::vector<float>* m_genak4chgjet_pt;
  std::vector<float>* m_genak4chgjet_p;
  std::vector<float>* m_genak4chgjet_z;
  std::vector<float>* m_PFjet_vz;
  std::vector<float>* m_PFjet_p;
  std::vector<float>* m_PFjet_phi;
  std::vector<float>* m_PFjet_eta;
  std::vector<float>* m_PFjet_pt;

  std::vector<float>* m_2ltrkjet_vz;
  std::vector<float>* m_2ltrkjet_p;
  std::vector<float>* m_2ltrkjet_phi;
  std::vector<float>* m_2ltrkjet_eta;
  std::vector<float>* m_2ltrkjet_pt;
  std::vector<int>* m_2ltrkjet_ntracks;
  std::vector<int>* m_2ltrkjet_ndtrk;
  std::vector<int>* m_2ltrkjet_nttrk;
  std::vector<int>* m_2ltrkjet_ntdtrk;

  std::vector<float>* m_trkjet_vz;
  std::vector<float>* m_trkjet_p;
  std::vector<float>* m_trkjet_phi;
  std::vector<float>* m_trkjet_eta;
  std::vector<float>* m_trkjet_pt;
  std::vector<int>* m_trkjet_ntracks;
  std::vector<float>* m_trkjet_tp_sumpt;
  std::vector<float>* m_trkjet_truetp_sumpt;
  std::vector<float>* m_calotrkjet_vz;
  std::vector<float>* m_calotrkjet_p;
  std::vector<float>* m_calotrkjet_phi;
  std::vector<float>* m_calotrkjet_eta;
  std::vector<float>* m_calotrkjet_pt;

  std::vector<float>* m_calojet_p;
  std::vector<float>* m_calojet_phi;
  std::vector<float>* m_calojet_eta;
  std::vector<float>* m_calojet_pt;
  std::vector<float>* m_tpjet_vz;
  std::vector<float>* m_tpjet_p;
  std::vector<float>* m_tpjet_phi;
  std::vector<float>* m_tpjet_eta;
  std::vector<int>* m_tpjet_ntracks;
  std::vector<float>* m_tpjet_tp_sumpt;
  std::vector<float>* m_tpjet_truetp_sumpt;
  std::vector<float>* m_tpjet_pt;
};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1TrackJetFastProducer::L1TrackJetFastProducer(ParameterSet const& iConfig) : 
  config(iConfig)
{

  MyProcess        = iConfig.getParameter< int >("MyProcess");
  DebugMode        = iConfig.getParameter< bool >("DebugMode");
  SaveAllTracks    = iConfig.getParameter< bool >("SaveAllTracks");
  L1Tk_nPar        = iConfig.getParameter< int >("L1Tk_nPar");
  //Cuts on Tracks
  NPS_minStubs= iConfig.getParameter<int>("nStubsPSmin");
  TP_minNStub      = iConfig.getParameter< int >("TP_minNStub");
  TP_minNStubLayer = iConfig.getParameter< int >("TP_minNStubLayer");
  TP_minPt         = iConfig.getParameter< double >("PTMINTRA");
  TP_maxEta        = 2.4;
  TP_maxZ0         = iConfig.getParameter< double >("ZMAX");
  DeltaZ0Cut        = iConfig.getParameter< double >("DeltaZ0Cut");
  CHI2MAX 	    = iConfig.getParameter< double >("CHI2MAX");
  PTMAX 	   =iConfig.getParameter< double >("PTMAX");
  CONESize	    =iConfig.getParameter<double>("CONESize");
  L1TrackInputTag      = iConfig.getParameter<InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<InputTag>("MCTruthTrackInputTag");
  RecoVertexInputTag    = iConfig.getParameter<InputTag>("RecoVertexInputTag");
  GenParticleInputTag   = iConfig.getParameter<InputTag >("GenParticleInputTag");
  MCTruthStubInputTag = iConfig.getParameter<InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = iConfig.getParameter<InputTag>("TrackingParticleInputTag");
  CaloJetInputTag = iConfig.getParameter<InputTag>("CaloJetInputTag");
  CaloTkJetInputTag = iConfig.getParameter<InputTag>("CaloTkJetInputTag");
  TkJetInputTag = iConfig.getParameter<InputTag>("TkJetInputTag");
  //TwoLayerTkJetInputTag = iConfig.getParameter<InputTag>("TwoLayerJetInputTag");
  GenJetAK4=iConfig.getParameter<InputTag>("GenJetAK4");
  PFJetInputTag=iConfig.getParameter<InputTag>("PFJetInputTag");
  HEPMCVertexToken_=consumes< std::vector< reco::GenParticle> >(GenParticleInputTag);
  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthTrackInputTag);
  ttStubMCTruthToken_ = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);
  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);
  GenJetCollectionToken_=consumes< std::vector<reco::GenJet > >(GenJetAK4);
  L1VertexToken_=consumes<VertexCollection>(RecoVertexInputTag);
  L1TkJetsToken_=consumes<L1TkJetParticleCollection>(TkJetInputTag);
  //L1TwoLayerTkJetsToken_=consumes<L1TkJetParticleCollection>(TwoLayerTkJetInputTag);
  L1CaloTkJetsToken_=consumes<L1TkJetParticleCollection>(CaloTkJetInputTag);
  L1PFCaloClusterToken_=consumes<reco::PFJetCollection>(CaloJetInputTag);
  L1PFJetToken_=consumes<reco::PFJetCollection>(PFJetInputTag);
}

/////////////
// DESTRUCTOR
L1TrackJetFastProducer::~L1TrackJetFastProducer()
{
}  

//////////
// END JOB
void L1TrackJetFastProducer::endJob()
{
  // things to be done at the exit of the event Loop
  cerr << "L1TrackJetFastProducer::endJob" << endl;

}

////////////
// BEGIN JOB
void L1TrackJetFastProducer::beginJob()
{

  // things to be done before entering the event Loop
  cerr << "L1TrackJetFastProducer::beginJob" << endl;

  //-----------------------------------------------------------------------------------------------
  // book histograms / make ntuple
  Service<TFileService> fs;


  // initilize
  m_trk_p    = new std::vector<float>;
  m_trk_pt    = new std::vector<float>;
  m_trk_eta   = new std::vector<float>;
  m_trk_phi   = new std::vector<float>;
  m_trk_z0    = new std::vector<float>;
  m_trk_d0    = new std::vector<float>;
  m_trk_chi2  = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
  m_trk_psnstub = new std::vector<int>;
  m_trk_genuine      = new std::vector<int>;
  m_trk_loose        = new std::vector<int>;
  m_trk_unknown      = new std::vector<int>;
  m_trk_combinatoric = new std::vector<int>;
  m_trk_fake = new std::vector<int>;
  m_trk_matchtp_pdgid = new std::vector<int>;
  m_trk_matchtp_pt = new std::vector<float>;
  m_trk_matchtp_eta = new std::vector<float>;
  m_trk_matchtp_phi = new std::vector<float>;
  m_trk_matchtp_z0 = new std::vector<float>;
  m_trk_matchtp_dxy = new std::vector<float>;
  m_trk_bconsist=new std::vector<float>;
  m_trk_sconsist=new std::vector<float>;
  m_tp_p     = new std::vector<float>;
  m_tp_pt     = new std::vector<float>;
  m_tp_eta    = new std::vector<float>;
  m_tp_phi    = new std::vector<float>;
  m_tp_dxy    = new std::vector<float>;
  m_tp_d0     = new std::vector<float>;
  m_tp_z0     = new std::vector<float>;
  m_tp_d0_prod = new std::vector<float>;
  m_tp_z0_prod = new std::vector<float>;
  m_tp_pdgid  = new std::vector<int>;
  m_tp_nmatch = new std::vector<int>;
  m_tp_nstub  = new std::vector<int>;
  m_tp_nstublayers  = new std::vector<int>;
  m_tp_eventid = new std::vector<int>;
  m_tp_charge = new std::vector<int>;
/*
  m_matchtrk_p    = new std::vector<float>;
  m_matchtrk_pt    = new std::vector<float>;
  m_matchtrk_eta   = new std::vector<float>;
  m_matchtrk_phi   = new std::vector<float>;
  m_matchtrk_z0    = new std::vector<float>;
  m_matchtrk_d0    = new std::vector<float>;
  m_matchtrk_chi2  = new std::vector<float>;
  m_matchtrk_nstub = new std::vector<int>;
  */
  m_pv_L1recotruesumpt = new std::vector<float>;
  m_pv_L1recosumpt = new std::vector<float>;
  m_pv_L1reco = new std::vector<float>;
  m_pv_L1TP = new std::vector<float>;
  m_pv_L1TPsumpt = new std::vector<float>;
  m_pv_MC = new std::vector<float>;
  m_pv_MCChgSumpT = new std::vector<float>;
  m_MC_lep=new std::vector<int>;  
  m_genak4jet_phi = new std::vector<float>;
  m_genak4jet_neufrac = new std::vector<float>;
  m_genak4jet_chgfrac = new std::vector<float>;
  m_genak4jet_metfrac = new std::vector<float>;
  m_genak4jet_eta = new std::vector<float>;
  m_genak4jet_pt = new std::vector<float>;
  m_genak4jet_p = new std::vector<float>;

  m_genak4chgjet_phi = new std::vector<float>;
  m_genak4chgjet_eta = new std::vector<float>;
  m_genak4chgjet_pt = new std::vector<float>;
  m_genak4chgjet_p = new std::vector<float>;
  m_genak4chgjet_z = new std::vector<float>;

  m_PFjet_eta = new std::vector<float>;
  m_PFjet_vz = new std::vector<float>;
  m_PFjet_phi = new std::vector<float>;
  m_PFjet_p = new std::vector<float>;
  m_PFjet_pt = new std::vector<float>;
  m_2ltrkjet_eta = new std::vector<float>;
  m_2ltrkjet_vz = new std::vector<float>;
  m_2ltrkjet_phi = new std::vector<float>;
  m_2ltrkjet_p = new std::vector<float>;
  m_2ltrkjet_pt = new std::vector<float>;
  m_2ltrkjet_ntracks=new std::vector<int>;
  m_2ltrkjet_ndtrk=new std::vector<int>;
  m_2ltrkjet_nttrk=new std::vector<int>;
   m_2ltrkjet_ntdtrk=new std::vector<int>;

  m_trkjet_eta = new std::vector<float>;
  m_trkjet_vz = new std::vector<float>;
  m_trkjet_phi = new std::vector<float>;
  m_trkjet_p = new std::vector<float>;
  m_trkjet_pt = new std::vector<float>;
  m_trkjet_ntracks = new std::vector<int>;
  m_trkjet_tp_sumpt = new std::vector<float>;
  m_trkjet_truetp_sumpt = new std::vector<float>;

  m_calotrkjet_vz=new std::vector<float>;
  m_calotrkjet_p=new std::vector<float>;
  m_calotrkjet_phi=new std::vector<float>;
  m_calotrkjet_eta=new std::vector<float>;
  m_calotrkjet_pt=new std::vector<float>;

  m_calojet_p=new std::vector<float>;
  m_calojet_phi=new std::vector<float>;
  m_calojet_eta=new std::vector<float>;
  m_calojet_pt=new std::vector<float>;
  
  m_tpjet_eta = new std::vector<float>;
  m_tpjet_vz = new std::vector<float>;
  m_tpjet_phi = new std::vector<float>;
  m_tpjet_p = new std::vector<float>;
  m_tpjet_pt = new std::vector<float>;
  m_tpjet_ntracks = new std::vector<int>;
  m_tpjet_tp_sumpt = new std::vector<float>;
  m_tpjet_truetp_sumpt = new std::vector<float>;
  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");
  
  if (SaveAllTracks) {
    eventTree->Branch("trk_p",    &m_trk_p);
    eventTree->Branch("trk_pt",    &m_trk_pt);
    eventTree->Branch("trk_eta",   &m_trk_eta);
    eventTree->Branch("trk_phi",   &m_trk_phi);
    eventTree->Branch("trk_d0",    &m_trk_d0);
    eventTree->Branch("trk_z0",    &m_trk_z0);
    eventTree->Branch("trk_chi2",  &m_trk_chi2);
    eventTree->Branch("trk_nstub", &m_trk_nstub);
    eventTree->Branch("trk_psnstub", &m_trk_psnstub);
    eventTree->Branch("trk_sconsist", &m_trk_sconsist);
    eventTree->Branch("trk_bconsist", &m_trk_bconsist);

    eventTree->Branch("trk_genuine",      &m_trk_genuine);
    eventTree->Branch("trk_loose",        &m_trk_loose);
    eventTree->Branch("trk_unknown",      &m_trk_unknown);
    eventTree->Branch("trk_combinatoric", &m_trk_combinatoric);
    eventTree->Branch("trk_fake", &m_trk_fake);
    eventTree->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
    eventTree->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
    eventTree->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
    eventTree->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
    eventTree->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
    eventTree->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);

  eventTree->Branch("tp_p",     &m_tp_p);
  eventTree->Branch("tp_pt",     &m_tp_pt);
  eventTree->Branch("tp_eta",    &m_tp_eta);
  eventTree->Branch("tp_phi",    &m_tp_phi);
  eventTree->Branch("tp_dxy",    &m_tp_dxy);
  eventTree->Branch("tp_d0",     &m_tp_d0);
  eventTree->Branch("tp_z0",     &m_tp_z0);
  eventTree->Branch("tp_d0_prod",&m_tp_d0_prod);
  eventTree->Branch("tp_z0_prod",&m_tp_z0_prod);
  eventTree->Branch("tp_pdgid",  &m_tp_pdgid);
  eventTree->Branch("tp_nmatch", &m_tp_nmatch);
  eventTree->Branch("tp_nstub", &m_tp_nstub);
  eventTree->Branch("tp_nstublayers", &m_tp_nstublayers);
  eventTree->Branch("tp_eventid",&m_tp_eventid);
  eventTree->Branch("tp_charge",&m_tp_charge);
/*
  eventTree->Branch("matchtrk_p",      &m_matchtrk_p);
  eventTree->Branch("matchtrk_pt",      &m_matchtrk_pt);
  eventTree->Branch("matchtrk_eta",     &m_matchtrk_eta);
  eventTree->Branch("matchtrk_phi",     &m_matchtrk_phi);
  eventTree->Branch("matchtrk_z0",      &m_matchtrk_z0);
  eventTree->Branch("matchtrk_d0",      &m_matchtrk_d0);
  eventTree->Branch("matchtrk_chi2",    &m_matchtrk_chi2);
  eventTree->Branch("matchtrk_nstub",   &m_matchtrk_nstub);
 */
  }

    //eventTree->Branch("pv_L1recofakesumpt", &m_pv_L1recofakesumpt);
    eventTree->Branch("pv_L1recotruesumpt", &m_pv_L1recotruesumpt);
    eventTree->Branch("pv_L1recosumpt", &m_pv_L1recosumpt);
    eventTree->Branch("pv_L1reco", &m_pv_L1reco);
    eventTree->Branch("pv_L1TP", &m_pv_L1TP);
    eventTree->Branch("pv_L1TPsumpt", &m_pv_L1TPsumpt);
    eventTree->Branch("MC_lep", &m_MC_lep);
    eventTree->Branch("pv_MCChgSumpT", &m_pv_MCChgSumpT);
    eventTree->Branch("pv_MC", &m_pv_MC);
    eventTree->Branch("tpjet_eta", &m_tpjet_eta);
    eventTree->Branch("tpjet_vz", &m_tpjet_vz);
    eventTree->Branch("tpjet_p", &m_tpjet_p);
    eventTree->Branch("tpjet_pt", &m_tpjet_pt);
    eventTree->Branch("tpjet_phi", &m_tpjet_phi);
    eventTree->Branch("tpjet_ntracks", &m_tpjet_ntracks);
    eventTree->Branch("tpjet_tp_sumpt", &m_tpjet_tp_sumpt);
    eventTree->Branch("tpjet_truetp_sumpt", &m_tpjet_truetp_sumpt);
    eventTree->Branch("calojet_eta", &m_calojet_eta);
    eventTree->Branch("calojet_p", &m_calojet_p);
    eventTree->Branch("calojet_pt", &m_calojet_pt);
    eventTree->Branch("calojet_phi", &m_calojet_phi);
    eventTree->Branch("calotrkjet_eta", &m_calotrkjet_eta);
    eventTree->Branch("calotrkjet_p", &m_calotrkjet_p);
    eventTree->Branch("calotrkjet_pt", &m_calotrkjet_pt);
    eventTree->Branch("calotrkjet_phi", &m_calotrkjet_phi);
    eventTree->Branch("calotrkjet_vz", &m_calotrkjet_vz);
    eventTree->Branch("PFjet_eta", &m_PFjet_eta);
    eventTree->Branch("PFjet_vz", &m_PFjet_vz);
    eventTree->Branch("PFjet_p", &m_PFjet_p);
    eventTree->Branch("PFjet_pt", &m_PFjet_pt);
    eventTree->Branch("PFjet_phi", &m_PFjet_phi);
    eventTree->Branch("2ltrkjet_eta", &m_2ltrkjet_eta);
    eventTree->Branch("2ltrkjet_vz", &m_2ltrkjet_vz);
    eventTree->Branch("2ltrkjet_p", &m_2ltrkjet_p);
    eventTree->Branch("2ltrkjet_pt", &m_2ltrkjet_pt);
    eventTree->Branch("2ltrkjet_phi", &m_2ltrkjet_phi);
    eventTree->Branch("2ltrkjet_ntracks", &m_2ltrkjet_ntracks);
    eventTree->Branch("trkjet_eta", &m_trkjet_eta);
    eventTree->Branch("trkjet_vz", &m_trkjet_vz);
    eventTree->Branch("trkjet_p", &m_trkjet_p);
    eventTree->Branch("trkjet_pt", &m_trkjet_pt);
    eventTree->Branch("trkjet_phi", &m_trkjet_phi);
    eventTree->Branch("trkjet_ntracks", &m_trkjet_ntracks);
    eventTree->Branch("trkjet_truetp_sumpt", m_trkjet_truetp_sumpt);
    eventTree->Branch("genjetak4_neufrac", &m_genak4jet_neufrac);
    eventTree->Branch("genjetak4_chgfrac", &m_genak4jet_chgfrac);
    eventTree->Branch("genjetak4_metfrac", &m_genak4jet_metfrac);
    eventTree->Branch("genjetak4_eta", &m_genak4jet_eta);
    eventTree->Branch("genjetak4_phi", &m_genak4jet_phi);
    eventTree->Branch("genjetak4_p", &m_genak4jet_p);
    eventTree->Branch("genjetak4_pt", &m_genak4jet_pt);
    eventTree->Branch("genjetchgak4_eta", &m_genak4chgjet_eta);
    eventTree->Branch("genjetchgak4_phi", &m_genak4chgjet_phi);
    eventTree->Branch("genjetchgak4_p", &m_genak4chgjet_p);
    eventTree->Branch("genjetchgak4_z", &m_genak4chgjet_z);
    eventTree->Branch("genjetchgak4_pt", &m_genak4chgjet_pt);
}


//////////
// ANALYZE
void L1TrackJetFastProducer::analyze(const Event& iEvent, const EventSetup& iSetup)
{


  if ( !(L1Tk_nPar==4 || L1Tk_nPar==5) ) {
    cout << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar << " but only 4/5 are valid options! Exiting..." << endl;
    return;
  }
  // clear variables
  if (SaveAllTracks) {
    m_trk_p->clear();
    m_trk_pt->clear();
    m_trk_eta->clear();
    m_trk_phi->clear();
    m_trk_d0->clear();
    m_trk_z0->clear();
    m_trk_chi2->clear();
    m_trk_nstub->clear();
    m_trk_psnstub->clear();
    m_trk_sconsist->clear();
    m_trk_bconsist->clear();
    m_trk_genuine->clear();
    m_trk_loose->clear();
    m_trk_unknown->clear();
    m_trk_combinatoric->clear();
    m_trk_fake->clear();
   
    m_trk_matchtp_pdgid->clear();
    m_trk_matchtp_pt->clear();
    m_trk_matchtp_eta->clear();
    m_trk_matchtp_phi->clear();
    m_trk_matchtp_z0->clear();
    m_trk_matchtp_dxy->clear();
   }
  
  m_tp_p->clear();
  m_tp_pt->clear();
  m_tp_eta->clear();
  m_tp_phi->clear();
  m_tp_dxy->clear();
  m_tp_d0->clear();
  m_tp_z0->clear();
  m_tp_d0_prod->clear();
  m_tp_z0_prod->clear();
  m_tp_pdgid->clear();
  m_tp_nmatch->clear();
  m_tp_nstub->clear();
  m_tp_nstublayers->clear();
  m_tp_eventid->clear();
  m_tp_charge->clear();
  m_genak4chgjet_phi->clear();
  m_genak4chgjet_eta->clear() ;
  m_genak4chgjet_pt->clear() ;
  m_genak4chgjet_p->clear() ;
  m_genak4chgjet_z->clear() ;
  m_genak4jet_phi->clear();
  m_genak4jet_neufrac->clear() ;
  m_genak4jet_chgfrac->clear() ;
  m_genak4jet_metfrac->clear() ;
  m_genak4jet_eta->clear() ;
  m_genak4jet_pt->clear() ;
  m_genak4jet_p->clear() ;
  m_PFjet_eta->clear();
  m_PFjet_pt->clear();
  m_PFjet_vz->clear();
  m_PFjet_phi->clear();
  m_PFjet_p->clear();
  m_2ltrkjet_eta->clear();
  m_2ltrkjet_pt->clear();
  m_2ltrkjet_vz->clear();
  m_2ltrkjet_phi->clear();
  m_2ltrkjet_p->clear();
  m_2ltrkjet_ntracks->clear();
  m_2ltrkjet_ndtrk->clear();
  m_2ltrkjet_nttrk->clear();
  m_2ltrkjet_ntdtrk->clear();

  m_trkjet_eta->clear();
  m_trkjet_pt->clear();
  m_trkjet_vz->clear();
  m_trkjet_phi->clear();
  m_trkjet_p->clear();
  m_trkjet_ntracks->clear();
  m_trkjet_truetp_sumpt->clear();
  m_trkjet_tp_sumpt->clear();
  m_tpjet_eta->clear();
  m_tpjet_pt->clear();
  m_tpjet_vz->clear();
  m_tpjet_phi->clear();
  m_tpjet_p->clear();
  m_tpjet_ntracks->clear();
  m_tpjet_tp_sumpt->clear();
  m_tpjet_truetp_sumpt->clear();
  //m_pv_L1recofakesumpt->clear();
  m_pv_L1recotruesumpt->clear();
  m_pv_L1recosumpt->clear();
  m_pv_L1reco->clear();
  m_pv_L1TPsumpt->clear();
  m_pv_L1TP->clear();
  m_pv_MC->clear();
  m_MC_lep->clear();
  m_pv_MCChgSumpT->clear();
  // -----------------------------------------------------------------------------------------------
  /*
  zbincount.clear();
  ttrk.clear();
  tdtrk.clear();
  ttdtrk.clear();
  L1TwoLayerInputPtrs.clear();
*/ 
 // retrieve various containers
  // -----------------------------------------------------------------------------------------------

  // L1 tracks
  Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);
  Handle<std::vector<reco::GenJet> >GenJetsAK4Handle;
  iEvent.getByToken(GenJetCollectionToken_,GenJetsAK4Handle); 

   Handle< std::vector< reco::GenParticle> > GenParticleHandle;
   iEvent.getByToken(HEPMCVertexToken_,GenParticleHandle);
    	int leptonicCount=0;
if(GenParticleHandle.isValid()){  
      vector<reco::GenParticle>::const_iterator genpartIter ;
        for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
		if (abs(genpartIter ->pdgId())!=24)continue;
		//std::cout<<"W-boson mother "<<genpartIter ->mother(0)->pdgId()<<std::endl;

		if(abs(genpartIter ->mother(0)->pdgId())!=6 && abs(genpartIter ->mother(0)->pdgId())!=24 )continue;
		if( abs(genpartIter ->daughter(0)->pdgId())==11  ||  abs(genpartIter ->daughter(0)->pdgId())==13 ||  abs(genpartIter ->daughter(0)->pdgId())==15){++leptonicCount;}
	}
	float zvtx_gen = -999;
	for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
           int status = genpartIter -> status() ;
       	   if (status< 21 || status>29) continue;
       //if (status!=1) continue;
              if ( genpartIter -> numberOfMothers() == 0) continue;
                     zvtx_gen = genpartIter -> vz() ;
       //                     
                                break;
        }
	m_pv_MC->push_back(zvtx_gen);
 }
  m_MC_lep->push_back(leptonicCount);
  // MC truth association maps
   Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
   iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);
 
  Handle<reco::PFJetCollection> PFCaloJetHandle;
  iEvent.getByToken(L1PFCaloClusterToken_,PFCaloJetHandle); //L1 Calo Jets made from PF Clusters;

 
  Handle<reco::PFJetCollection> PFJetHandle;
  iEvent.getByToken(L1PFJetToken_,PFJetHandle); //L1 Calo Jets made from PF Clusters;
  
 // Handle<L1TkJetParticleCollection> TwoLayerTkJetHandle;
  //iEvent.getByToken(L1TwoLayerTkJetsToken_,TwoLayerTkJetHandle); 
 
  Handle<L1TkJetParticleCollection> TkJetHandle;
  iEvent.getByToken(L1TkJetsToken_,TkJetHandle); 
 
  Handle<L1TkJetParticleCollection> CaloTkJetHandle;
  iEvent.getByToken(L1CaloTkJetsToken_,CaloTkJetHandle); //L1 Calo +tk Jets made from TP algo;

  Handle< VertexCollection >L1TkPrimaryVertexHandle;
  iEvent.getByToken(L1VertexToken_, L1TkPrimaryVertexHandle);
 if(L1TkPrimaryVertexHandle.isValid()){
	m_pv_L1reco->push_back(L1TkPrimaryVertexHandle->begin()->z0());
	//add Vertex True quality, Vertex Fake Content, total sumpT, and number of tracks 
	std::vector<Ptr< TTTrack<Ref_Phase2TrackerDigi_> > >Vtxtracks=L1TkPrimaryVertexHandle->begin()->tracks();
	//std::cout<<"Ntracks in Vertex "<<Vtxtracks.size()<<std::endl;
	float sumpt=0;
	float trueContent=0;
	float fakeContent=0;
	for(unsigned int t=0; 	t<Vtxtracks.size(); ++t){
		sumpt=Vtxtracks[t]->getMomentum(L1Tk_nPar).perp()+sumpt;
		//if(Vtxtracks[t].isAvailable())std::cout<<"Has a collection in memory "<<Vtxtracks[t].id()<<std::endl;
		//Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(Vtxtracks[t]);	
		//if(MCTruthTTTrackHandle->isGenuine(l1track_ptr))std::cout<<"MC Handle is valid "<<std::endl;//trueContent=trueContent+Vtxtracks[t]->getMomentum(L1Tk_nPar).perp();
		if(MCTruthTTTrackHandle->isGenuine(Vtxtracks[t]))trueContent=trueContent+Vtxtracks[t]->getMomentum(L1Tk_nPar).perp();
		else fakeContent=fakeContent+Vtxtracks[t]->getMomentum(L1Tk_nPar).perp();
	}
	m_pv_L1recotruesumpt->push_back(trueContent);
	m_pv_L1recosumpt->push_back(sumpt);	
}
  // -----------------------------------------------------------------------------------------------
  // more for TTStubs
  ESHandle<TrackerGeometry> geometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);

  ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

  ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
//  if(!magneticFieldHandle.isValid())std::cout<<" Mag field not present "<<std::endl;
//  //  else {
    const MagneticField* theMagneticField = magneticFieldHandle.product();
      double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
//        //std::cout<<" Mag field "<<mMagneticFieldStrength<<std::endl;

//Gen Jet branches
std::vector<float>genConstpt;
std::vector<float>genConsteta;
std::vector<float>genConstphi;
std::vector<float>genConstp;
std::vector<float>genConstz;
std::vector<int>genID;
if(GenJetsAK4Handle.isValid()){
        for (unsigned iGenJet = 0; iGenJet < GenJetsAK4Handle->size(); ++iGenJet) {

         const reco::GenJet& genJet = (*GenJetsAK4Handle) [iGenJet];
   if(fabs(genJet.eta())<TP_maxEta){
          m_genak4jet_p->push_back(genJet.energy());
          m_genak4jet_pt->push_back(genJet.pt());
	  m_genak4jet_metfrac->push_back(genJet.invisibleEnergy()/genJet.energy());
          m_genak4jet_eta->push_back(genJet.eta());
          m_genak4jet_phi->push_back(genJet.phi());
	float NeuEnergy=0;
	float ChgEnergy=0;

	for(unsigned g=0; g<genJet.getGenConstituents().size(); ++g){
		if(genJet.getGenConstituent(g)->charge()!=0 && genJet.getGenConstituent(g)->pt()>2){
		genConstpt.push_back(genJet.getGenConstituent(g)->pt());
		genConsteta.push_back(genJet.getGenConstituent(g)->eta());
		genConstphi.push_back(genJet.getGenConstituent(g)->phi());
		genConstp.push_back(genJet.getGenConstituent(g)->p());
		genConstz.push_back(genJet.getGenConstituent(g)->vz());	
		ChgEnergy=ChgEnergy+genJet.getGenConstituent(g)->energy();
		 genID.push_back(1);
	}
		else NeuEnergy=NeuEnergy+genJet.getGenConstituent(g)->energy();
        	}
		m_genak4jet_chgfrac->push_back(ChgEnergy);
	        m_genak4jet_neufrac->push_back(NeuEnergy);
    }
  }
//fill chg particle fastjets
JetOutputs_.clear();
JetVz_.clear();
JetNtracks_.clear();
JetEventID_.clear();
FillFastJets(genConstpt, genConsteta, genConstphi, genConstp, genConstz, genID, CONESize, JetOutputs_, JetNtracks_, JetVz_,JetEventID_);
for (unsigned int ijet=0;ijet<JetOutputs_.size();++ijet) {
  m_genak4chgjet_phi->push_back(JetOutputs_[ijet].phi_std());
  m_genak4chgjet_eta->push_back(JetOutputs_[ijet].eta());
  m_genak4chgjet_pt->push_back(JetOutputs_[ijet].pt());
  m_genak4chgjet_p->push_back(JetOutputs_[ijet].modp());
  m_genak4chgjet_z->push_back(JetVz_[ijet]);
  }
}
  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------

  if (SaveAllTracks) {
    
    if (DebugMode) {
      cout << endl << "Loop over L1 tracks!" << endl;
      cout << endl << "Looking at " << L1Tk_nPar << "-parameter tracks!" << endl;
    }
    
    int this_l1track = 0;
    std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;

    for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
      Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);

      this_l1track++;
      float tmp_trk_p   = iterL1Track->getMomentum(L1Tk_nPar).mag();
      float tmp_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
      float tmp_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
      float tmp_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
      float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm
          float tmp_trk_stubPtConsistency = StubPtConsistency::getConsistency(TTTrackHandle->at(this_l1track-1), theTrackerGeom, tTopo,mMagneticFieldStrength,L1Tk_nPar); 
      float tmp_trk_d0 = -999;
      if (L1Tk_nPar == 5) {
	float tmp_trk_x0   = iterL1Track->getPOCA(L1Tk_nPar).x();
	float tmp_trk_y0   = iterL1Track->getPOCA(L1Tk_nPar).y();	
	tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
      }

      float tmp_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);
      int tmp_trk_nstub  = (int) iterL1Track->getStubRefs().size();


      if (tmp_trk_pt < TP_minPt) continue;
      if (fabs(tmp_trk_eta) > TP_maxEta) continue;
      if (fabs(tmp_trk_z0) > TP_maxZ0) continue;
        int nPS=0;
        float tmp_trk_signedPt = 0.3*3.811202/100.0/(iterL1Track->getRInv());
	        std::vector< Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > >  theStubs =  iterL1Track-> getStubRefs() ;
           float tmp_trk_bendchi2 = 0;    
  float sigma_bend1 = 0.463;
  float sigma_bend2 = 0.463;
    for (unsigned int istub=0; istub<(unsigned int)theStubs.size(); istub++) {
	  DetId detIdStub = theTrackerGeom->idToDet( ( theStubs.at(istub)->getClusterRef(0))->getDetId() )->geographicalId();	
          MeasurementPoint coords = theStubs.at(istub)->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	 const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detIdStub );
        const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
        Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );
        bool isBarrel = false;
        int layer=-999999;
        if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
          isBarrel = true;
          layer  = static_cast<int>(tTopo->layer(detIdStub));
        }
        else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
          isBarrel = false;
          layer  = static_cast<int>(tTopo->layer(detIdStub));
          layer+=5;
        }

        float pitch = 0.089;

        if (theTrackerGeom->getDetectorType(detIdStub)==TrackerGeometry::ModuleType::Ph2PSP){
          pitch = 0.099;
	  ++nPS;
        }
        double tmp_stub_z=posStub.z();
        double tmp_stub_r=posStub.perp();
        //double tmp_stub_phi=posStub.phi();
        //TMTT's qOverPt variables
        const GeomDetUnit* det1 = theTrackerGeom->idToDetUnit( tTopo->partnerDetId( detIdStub ) );
        const PixelGeomDetUnit* unit = reinterpret_cast<const PixelGeomDetUnit*>( det0 );
        const PixelTopology& topo = unit->specificTopology();
        bool tiltedBarrel = (isBarrel && tTopo->tobSide(detIdStub)!=3);
        float stripPitch = topo.pitch().first;
        float R0 = det0->position().perp();
        float R1 = det1->position().perp();
        float Z0 = det0->position().z();
        float Z1 = det1->position().z();
        float modMinR = std::min(R0,R1);
        float modMaxR = std::max(R0,R1);
        float modMinZ = std::min(Z0,Z1);
        float modMaxZ = std::max(Z0,Z1);


        float sensorSpacing = sqrt((modMaxR-modMinR)*(modMaxR-modMinR) + (modMaxZ-modMinZ)*(modMaxZ-modMinZ));
        float pitchOverSep = stripPitch/sensorSpacing;//stripPitch = 0.01
        float dphiOverBendCorrection;
        if (tiltedBarrel) dphiOverBendCorrection= 0.886454*fabs(tmp_stub_z)/tmp_stub_r+0.504148;
        else if (isBarrel) dphiOverBendCorrection=1;
        else dphiOverBendCorrection= fabs(tmp_stub_z)/tmp_stub_r;
        //float dphiOverBend = pitchOverSep*dphiOverBendCorrection;
        float stubBend = theStubs.at(istub)->getTriggerBend();
        if (!isBarrel && tmp_stub_z<0.0) stubBend=-stubBend;
        //float qOverPt = -(stubBend*dphiOverBend)/(tmp_stub_r*mMagneticFieldStrength*(3.0E8/2.0E11));//0.57 = B*c/(2E11*strip pitch)
        float trackBend = -(sensorSpacing*0.57*tmp_stub_r/10)/(pitch*tmp_trk_signedPt*dphiOverBendCorrection);
        //float trigDisplace = theStubs.at(istub)->getTriggerDisplacement();
        //float trigOffset = theStubs.at(istub)->getTriggerOffset();
        //float trigPos = theStubs.at(istub)->getTriggerPosition();

       // float qOverPtDiff = 1.0/tmp_trk_signedPt - qOverPt;
        float tmp_bend_diff = trackBend-stubBend;

        //float trackBend = -(1.8*0.57*tmp_stub_r/100)/(pitch*tmp_matchtrk_pt);
        //if (stub_tp_charge < 0 && stub_tp_charge!=-999) trackBend = -trackBend;
        if (fabs(tmp_trk_signedPt)<4){
          tmp_trk_bendchi2 += (tmp_bend_diff*tmp_bend_diff)/(sigma_bend1*sigma_bend1);
        }
        else {
          tmp_trk_bendchi2 += (tmp_bend_diff*tmp_bend_diff)/(sigma_bend2*sigma_bend2);
        }
        }
      tmp_trk_bendchi2=tmp_trk_bendchi2/(tmp_trk_nstub); 
  
/*      for (unsigned int istub=0; istub<(unsigned int)theStubs.size(); istub++) {
          bool isPS = false;
	  DetId detIdStub = theTrackerGeom->idToDet( ( theStubs.at(istub)->getClusterRef(0))->getDetId() )->geographicalId();	
	        const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detIdStub );
      const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );
	if (topol->nrows() == 960) isPS=true;	
       //DetId detId( theStubs.at(istub)->getDetId() );
      // if (detId.det() == DetId::Detector::Tracker) {
       //  if (detId.subdetId() == StripSubdetector::TOB && tTopo->tobLayer(detId) <= 3)  isPS = true;
        // else if (detId.subdetId() == StripSubdetector::TID && tTopo->tidRing(detId) <= 9)  isPS = true;
      // }
       if (isPS) nPS ++;
        }
*/
      //std::cout<<"n PS hits "<<nPS<<std::endl;
      int tmp_trk_genuine = 0;
      int tmp_trk_loose = 0;
      int tmp_trk_unknown = 0;
      int tmp_trk_combinatoric = 0;
      int tmp_eventid=-1;
	/*
      if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
      if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
      if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;
      
      if (DebugMode) {
	cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi 
	     << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
	if (tmp_trk_genuine) cout << " (is genuine)" << endl; 
	if (tmp_trk_unknown) cout << " (is unknown)" << endl; 
	if (tmp_trk_combinatoric) cout << " (is combinatoric)" << endl; 
      }
	*/
      //if(nPS>=NPS_minStubs && tmp_trk_pt>TP_minPt && fabs(m_pv_L1reco->at(0)-tmp_trk_z0)<DeltaZ0Cut && TrackQualityCuts(tmp_trk_pt,tmp_trk_nstub,tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar))){
	//&& (tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar)<CHI2MAX  || tmp_trk_pt<20) ){ 

      //if(tmp_trk_pt>TP_minPt && fabs(m_pv_L1reco->at(0)-tmp_trk_z0)<DeltaZ0Cut && (tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar)<CHI2MAX  || tmp_trk_pt<20) ){ 
      //if(nPS>=NPS_minStubs && tmp_trk_pt>TP_minPt && fabs(m_pv_MC->at(0)-tmp_trk_z0)<DeltaZ0Cut && (tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar)<CHI2MAX  || tmp_trk_pt<20) ){ 
     ///if(tmp_trk_pt>TP_minPt && fabs(m_pv_L1reco->at(0)-tmp_trk_z0)<DeltaZ0Cut && (tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar)<CHI2MAX  || tmp_trk_pt<20) 
//	&& tmp_trk_nstub>=TP_minNStub ){
      //if(tmp_trk_pt>200)tmp_trk_pt=200;
      //L1TwoLayerInputPtrs.push_back(l1track_ptr);
      //zbincount.push_back(0);
/*
      float tmp_trk_zchi2        = iterL1Track->getChi2(4); //YG hack: for 5-par tracks, use 4par to store rz fit chi2
      float tmp_trk_bconsistency = iterL1Track->getStubPtConsistency(L1Tk_nPar); //YG hack: consistency
    if ((fabs(tmp_trk_d0)>0.1 && tmp_trk_nstub>=5)||(tmp_trk_nstub==4 && fabs(tmp_trk_d0)>1.0)) tdtrk.push_back(1);
    else tdtrk.push_back(0);

    if ( (tmp_trk_chi2/(tmp_trk_nstub-3)) <3.5 && (tmp_trk_zchi2/(tmp_trk_nstub-2))<2 && tmp_trk_nstub>=5 && tmp_trk_bconsistency<4) ttrk.push_back(1);
    else ttrk.push_back(0);
    if ( (tmp_trk_chi2/(tmp_trk_nstub-3)) <3.5 && (tmp_trk_zchi2/(tmp_trk_nstub-2))<2 && tmp_trk_nstub>=5 && tmp_trk_bconsistency<4 && ((fabs(tmp_trk_d0)>0.1) ))ttdtrk.push_back(1);
    else ttdtrk.push_back(0);
*/  
    if(fabs(m_pv_L1reco->at(0)-tmp_trk_z0)<DeltaZ0Cut){
      m_trk_p ->push_back(tmp_trk_p); 
      m_trk_pt ->push_back(tmp_trk_pt);
      m_trk_eta->push_back(tmp_trk_eta);
      m_trk_phi->push_back(tmp_trk_phi);
      m_trk_z0 ->push_back(tmp_trk_z0);
      m_trk_sconsist ->push_back(tmp_trk_stubPtConsistency);
      m_trk_bconsist ->push_back(tmp_trk_bendchi2);
      if (L1Tk_nPar==5) m_trk_d0->push_back(tmp_trk_d0);
      else m_trk_d0->push_back(999.);
      m_trk_chi2 ->push_back(tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar));
      m_trk_psnstub->push_back(nPS);
      m_trk_nstub->push_back(tmp_trk_nstub);
      edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);
      if(!my_tp.isNull()) tmp_eventid = my_tp->eventId().event();
      if (MCTruthTTTrackHandle->isGenuine(l1track_ptr) && tmp_eventid==0) tmp_trk_genuine = 1;
      m_trk_genuine->push_back(tmp_trk_genuine);
      m_trk_loose->push_back(tmp_trk_loose);
      m_trk_unknown->push_back(tmp_trk_unknown);
      m_trk_combinatoric->push_back(tmp_trk_combinatoric);
	}
     // }

      // ----------------------------------------------------------------------------------------------
      // for studying the fake rate
      // ----------------------------------------------------------------------------------------------

      Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

      int myFake = 0;

      int myTP_pdgid = -999;
      float myTP_pt = -999;
      float myTP_eta = -999;
      float myTP_phi = -999;
      float myTP_z0 = -999;
      float myTP_dxy = -999;


      if (my_tp.isNull()) myFake = 0;
      else {
	int tmp_eventid = my_tp->eventId().event();

	if (tmp_eventid > 0) myFake = 2;
	else myFake = 1;

	myTP_pdgid = my_tp->pdgId();
	myTP_pt = my_tp->p4().pt();
	myTP_eta = my_tp->p4().eta();
	myTP_phi = my_tp->p4().phi();
	myTP_z0 = my_tp->vertex().z();
	
	float myTP_x0 = my_tp->vertex().x();
	float myTP_y0 = my_tp->vertex().y();
	myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

	if (DebugMode) {
	  cout << "TP matched to track has pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta() 
	       << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z()
	       << " pdgid = " <<  my_tp->pdgId() << " dxy = " << myTP_dxy << endl;
	}
      }

      m_trk_fake->push_back(myFake);

      m_trk_matchtp_pdgid->push_back(myTP_pdgid);
      m_trk_matchtp_pt->push_back(myTP_pt);
      m_trk_matchtp_eta->push_back(myTP_eta);
      m_trk_matchtp_phi->push_back(myTP_phi);
      m_trk_matchtp_z0->push_back(myTP_z0);
      m_trk_matchtp_dxy->push_back(myTP_dxy);

    }//end track loop
JetOutputs_.clear();
fjConstituents_.clear();
JetVz_.clear();
JetNtracks_.clear();
JetEventID_.clear();
}//end if SaveAllTracks

  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------


  if (DebugMode) cout << endl << "Loop over tracking particles!" << endl;
  int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;

  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
 
    Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
   ++this_tp;
    std::vector< Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    //if(tp_ptr.isNull())std::cout<<"TP ptr is Null "<<std::endl;
    int nStubTP = (int) theStubRefs.size();
    //std::cout<<"Nstubs for TP "<<nStubTP<<std::endl;
    // how many layers/disks have stubs?
    int hasStubInLayer[11] = {0};
    for (unsigned int is=0; is<theStubRefs.size(); is++) {
	  
	  //detID of stub
	//  DetId detIdStub = theTrackerGeom->idToDet( (theStubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();
	DetId detIdStub( (theStubRefs.at(is)->getDetId() ));
	  
         int layer = -1;
	  if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
	    layer  = static_cast<int>(tTopo->layer(detIdStub));
	    layer=layer-1;
	  }
	
	  else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
	    layer  =static_cast<int>(tTopo->layer(detIdStub));
	    layer=layer+5;
	  }	  
	  else if ( detIdStub.subdetId()==StripSubdetector::TIB ) {
	    layer  =static_cast<int>(tTopo->layer(detIdStub));

	  }	  
	  else if ( detIdStub.subdetId()==StripSubdetector::TEC ) {
	    layer  =static_cast<int>(tTopo->layer(detIdStub));
	  }
	else{
		std::cout<<"I am a stub I don't know where i am; Layer "<<static_cast<int>(tTopo->layer(detIdStub))<<std::endl;
	}	  
         if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer]<2)hasStubInLayer[layer] = 1;
         else hasStubInLayer[layer] = 2;
}
    //                                                          //treat genuine stubs separately (==2 is genuine, ==1 is not)
    int nStubLayerTP = 0;
    int nStubLayerTP_g = 0;
    for (int isum=0; isum<11; isum++) {
        if ( hasStubInLayer[isum] >= 1) nStubLayerTP   += 1;
        if ( hasStubInLayer[isum] == 2) nStubLayerTP_g += 1;
   }
    int tmp_eventid = iterTP->eventId().event();
    float  tmp_tp_vx=tp_ptr->vx();
    float  tmp_tp_vy=tp_ptr->vy();
    float  tmp_tp_vz=tp_ptr->vz();
    float tmp_tp_eta=tp_ptr->eta();
    float tmp_tp_phi=tp_ptr->phi();
    float tmp_tp_pt=tp_ptr->pt();
    float tmp_tp_p=tp_ptr->p();
    float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));

    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;
	
    float A = 0.01*0.5696;
    float Kmagnitude = A / tmp_tp_pt;
    
    float tmp_tp_charge = tp_ptr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;
	
    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tmp_tp_phi));
    float tmp_tp_rp = sqrt(tmp_tp_x0p*tmp_tp_x0p + tmp_tp_y0p*tmp_tp_y0p);
    float tmp_tp_d0 = tmp_tp_charge*tmp_tp_rp - (1. / (2. * K));
	
    tmp_tp_d0 = tmp_tp_d0*(-1); //fix d0 sign

    static double pi = 4.0*atan(1.0);
    float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);

    int tmp_tp_pdgid = iterTP->pdgId();
    float tmp_tp_z0_prod = tmp_tp_vz;
    float tmp_tp_d0_prod = -tmp_tp_vx*sin(tmp_tp_phi) + tmp_tp_vy*cos(tmp_tp_phi);



    if (tmp_tp_pt < TP_minPt) continue;
    if (fabs(tmp_tp_eta) > TP_maxEta) continue;

    
    if (fabs(tmp_tp_z0) > TP_maxZ0) continue;


    // for pions in ttbar, only consider TPs coming from near the IP!
    float dxy = sqrt(tmp_tp_vx*tmp_tp_vx + tmp_tp_vy*tmp_tp_vy);
    float tmp_tp_dxy = dxy;

    if (DebugMode) cout << "Tracking particle, pt: " << tmp_tp_pt << " eta: " << tmp_tp_eta << " phi: " << tmp_tp_phi 
			<< " z0: " << tmp_tp_z0 << " d0: " << tmp_tp_d0 
			<< " z_prod: " << tmp_tp_z0_prod << " d_prod: " << tmp_tp_d0_prod 
			<< " pdgid: " << tmp_tp_pdgid << " eventID: " << iterTP->eventId().event()
			<< " ttstubs " << MCTruthTTStubHandle->findTTStubRefs(tp_ptr).size()
			<< " tttracks " << MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr).size() << endl;


    // ----------------------------------------------------------------------------------------------



    if (TP_minNStub > 0) {
      if (DebugMode) cout << "Only consider TPs with >= " << TP_minNStub << " stubs" << endl;
      if (nStubTP < TP_minNStub) {
	if (DebugMode) cout << "TP fails minimum nbr stubs requirement! Continuing..." << endl;
	continue;
      }
    }


    // ----------------------------------------------------------------------------------------------
    // look for L1 tracks matched to the tracking particle

    std::vector< Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > > matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);
    
    int nMatch = 0;
    int i_track = -1;
    float i_chi2dof = 99999;

    if (matchedTracks.size() > 0) { 
    
      if (DebugMode && (matchedTracks.size()>1)) cout << "TrackingParticle has more than one matched L1 track!" << endl;

      // ----------------------------------------------------------------------------------------------
      // loop over matched L1 tracks
      // here, "match" means tracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters 
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTTrack

      for (int it=0; it<(int)matchedTracks.size(); it++) {

	bool tmp_trk_genuine = false;
	if (MCTruthTTTrackHandle->isGenuine(matchedTracks.at(it))) tmp_trk_genuine = true;
	if (!tmp_trk_genuine) continue;


	if (DebugMode) {
	  if (MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it)).isNull()) {
	    cout << "track matched to TP is NOT uniquely matched to a TP" << endl;
	  }
	  else {
	    Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	    cout << "TP matched to track matched to TP ... tp pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta() 
		 << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z() << endl;
	  }
	  cout << "   ... matched L1 track has pt = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).perp() 
	       << " eta = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).eta()
	       << " phi = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).phi()
	       << " chi2 = " << matchedTracks.at(it)->getChi2(L1Tk_nPar) 
	       << " consistency = " << matchedTracks.at(it)->getStubPtConsistency(L1Tk_nPar) 
	       << " z0 = " << matchedTracks.at(it)->getPOCA(L1Tk_nPar).z() 
	       << " nstub = " << matchedTracks.at(it)->getStubRefs().size();
	  if (tmp_trk_genuine) cout << " (genuine!) " << endl;
	}


	// ----------------------------------------------------------------------------------------------
	// further require L1 track to be (loosely) genuine, that there is only one TP matched to the track
	// + have >= L1Tk_minNStub stubs for it to be a valid match (only relevant is your track collection
	// e.g. stores 3-stub tracks but at plot level you require >= 4 stubs (--> tracklet case)
	int tmp_trk_nstub = matchedTracks.at(it)->getStubRefs().size();
	if (tmp_trk_nstub < TP_minNStub) continue;
	
	float dmatch_pt  = 999;
	float dmatch_eta = 999;
	float dmatch_phi = 999;
	int match_id = 999;

	Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	dmatch_pt  = fabs(my_tp->p4().pt() - tmp_tp_pt);
	dmatch_eta = fabs(my_tp->p4().eta() - tmp_tp_eta);
	dmatch_phi = fabs(my_tp->p4().phi() - tmp_tp_phi);
	match_id = my_tp->pdgId();

	float tmp_trk_chi2dof = (matchedTracks.at(it)->getChi2(L1Tk_nPar)) / (2*tmp_trk_nstub - L1Tk_nPar);
	
	// ensure that track is uniquely matched to the TP we are looking at!
	if (dmatch_pt<0.1 && dmatch_eta<0.1 && dmatch_phi<0.1 && tmp_tp_pdgid==match_id) { 
	  nMatch++;
	  if (i_track < 0 || tmp_trk_chi2dof < i_chi2dof) {
	    i_track = it;
	    i_chi2dof = tmp_trk_chi2dof;
	  }
	}

      }// end loop over matched L1 tracks

    }// end has at least 1 matched L1 track
    // ----------------------------------------------------------------------------------------------

   /*
    float tmp_matchtrk_p   = -999;
    float tmp_matchtrk_pt   = -999;
    float tmp_matchtrk_eta  = -999;
    float tmp_matchtrk_phi  = -999;
    float tmp_matchtrk_z0   = -999;
    float tmp_matchtrk_d0   = -999;
    float tmp_matchtrk_chi2 = -999;
    //float tmp_matchtrk_p=-999;
    int tmp_matchtrk_nstub  = -999;


    if (nMatch > 1 && DebugMode) cout << "WARNING *** 2 or more matches to genuine L1 tracks ***" << endl;

    if (nMatch > 0) {
      tmp_matchtrk_pt   = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).perp();
      tmp_matchtrk_eta  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).eta();
      tmp_matchtrk_phi  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).phi();
      tmp_matchtrk_z0   = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).z();
      tmp_matchtrk_p = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).mag();
      if (L1Tk_nPar == 5) {
	float tmp_matchtrk_x0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).x();
	float tmp_matchtrk_y0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).y();
	tmp_matchtrk_d0 = -tmp_matchtrk_x0*sin(tmp_matchtrk_phi) + tmp_matchtrk_y0*cos(tmp_matchtrk_phi);
      }

      tmp_matchtrk_chi2 = matchedTracks.at(i_track)->getChi2(L1Tk_nPar);
      tmp_matchtrk_nstub  = (int) matchedTracks.at(i_track)->getStubRefs().size();
    }
   */
    //if( nMatch>0 &&  fabs(tmp_tp_eta)<TP_maxEta && fabs(m_pv_L1TP->at(0)-tmp_tp_z0)<DeltaZ0Cut && nStubLayerTP>=TP_minNStubLayer && tmp_tp_pt>TP_minPt&& abs(tmp_tp_pdgid)!=11){
    //if( fabs(tmp_tp_eta)<TP_maxEta && fabs(m_pv_L1TP->at(0)-tmp_tp_z0)<DeltaZ0Cut && nStubLayerTP>=TP_minNStubLayer && tmp_tp_pt>TP_minPt&& abs(tmp_tp_pdgid)!=11){
   if( tmp_tp_dxy<1.0 && fabs(tmp_tp_eta)<TP_maxEta && fabs(m_pv_MC->at(0)-tmp_tp_z0)<DeltaZ0Cut && nStubLayerTP>=TP_minNStubLayer && tmp_tp_pt>TP_minPt){
    m_tp_p->push_back(tmp_tp_p);
    m_tp_pt->push_back(tmp_tp_pt);
    m_tp_eta->push_back(tmp_tp_eta);
    m_tp_phi->push_back(tmp_tp_phi);
    m_tp_dxy->push_back(tmp_tp_dxy);
    m_tp_z0->push_back(tmp_tp_z0);
    m_tp_d0->push_back(tmp_tp_d0);
    m_tp_z0_prod->push_back(tmp_tp_z0_prod);
    m_tp_d0_prod->push_back(tmp_tp_d0_prod);
    m_tp_pdgid->push_back(tmp_tp_pdgid);
    m_tp_nmatch->push_back(nMatch);
    m_tp_nstub->push_back(nStubTP);
    m_tp_nstublayers->push_back(nStubLayerTP);
    if(tmp_eventid<=0)tmp_eventid=1;
    else tmp_eventid=0;
    m_tp_eventid->push_back(tmp_eventid);
    m_tp_charge->push_back(tmp_tp_charge);
/*
    m_matchtrk_p ->push_back(tmp_matchtrk_p);
    m_matchtrk_pt ->push_back(tmp_matchtrk_pt);
    m_matchtrk_eta->push_back(tmp_matchtrk_eta);
    m_matchtrk_phi->push_back(tmp_matchtrk_phi);
    m_matchtrk_z0 ->push_back(tmp_matchtrk_z0);
    m_matchtrk_d0 ->push_back(tmp_matchtrk_d0);
    m_matchtrk_chi2 ->push_back(tmp_matchtrk_chi2);
    m_matchtrk_nstub->push_back(tmp_matchtrk_nstub);
*/
     }
  } //end loop tracking particles
JetOutputs_.clear();
fjConstituents_.clear();
JetVz_.clear();
JetNtracks_.clear();
JetEventID_.clear();
FillFastJets(*m_tp_pt, *m_tp_eta, *m_tp_phi, *m_tp_p, *m_tp_z0, *m_tp_eventid, CONESize, JetOutputs_, JetNtracks_, JetVz_,JetEventID_);
for(unsigned int j=0; j<JetOutputs_.size(); ++j){
	TLorentzVector temp;
	temp.SetPtEtaPhiE(JetOutputs_[j].pt(),JetOutputs_[j].eta(),JetOutputs_[j].phi_std(),JetOutputs_[j].modp());
        m_tpjet_pt->push_back(temp.Pt());
        m_tpjet_p->push_back(temp.P());
        m_tpjet_eta->push_back(temp.Eta());
        m_tpjet_phi->push_back(temp.Phi());
        m_tpjet_vz->push_back(JetVz_[j]);
	m_tpjet_ntracks->push_back(JetNtracks_[j]);
	float sumpt=0;
	float truesumpt=0;
	for(unsigned int i=0; i<fjConstituents_.size(); ++i){
        	auto index =fjConstituents_[i].user_index();	
		int eventId=m_tp_eventid->at(index);
		if(eventId<=0)truesumpt=truesumpt+m_tp_pt->at(index);
		sumpt=sumpt+m_tp_pt->at(index);
	}
	m_tpjet_tp_sumpt->push_back(sumpt);
	m_tpjet_truetp_sumpt->push_back(truesumpt);
}
JetOutputs_.clear();
fjConstituents_.clear();
JetVz_.clear();
std::vector<L1TkJetParticle>::const_iterator jetIter;
for (jetIter = TkJetHandle->begin(); jetIter != TkJetHandle->end(); ++jetIter) {
        m_trkjet_vz->push_back(jetIter->getJetVtx());
        m_trkjet_ntracks->push_back(jetIter->getTrkPtrs().size());
        m_trkjet_phi->push_back(jetIter->phi());
        m_trkjet_eta->push_back(jetIter->eta());
        m_trkjet_pt->push_back(jetIter->pt());
        m_trkjet_p->push_back(jetIter->p());
}
/*
for (jetIter = TwoLayerTkJetHandle->begin(); jetIter != TwoLayerTkJetHandle->end(); ++jetIter) {
        m_2ltrkjet_vz->push_back(jetIter->getJetVtx());
        m_2ltrkjet_ntracks->push_back(jetIter->getTrkPtrs().size());
        m_2ltrkjet_phi->push_back(jetIter->phi());
        m_2ltrkjet_eta->push_back(jetIter->eta());
        m_2ltrkjet_pt->push_back(jetIter->pt());
        m_2ltrkjet_p->push_back(jetIter->p());
}
*/
/* 
  if(L1TwoLayerInputPtrs.size()>0){
    maxzbin mzb;

       L2_cluster(L1TwoLayerInputPtrs, ttrk, tdtrk,ttdtrk,mzb);
	 for(int k = 0; k < mzb.nclust; ++k){
	m_2ltrkjet_vz->push_back(mzb.zbincenter);
        m_2ltrkjet_nttrk->push_back(mzb.clusters[k].numttrks);
        m_2ltrkjet_ndtrk->push_back(mzb.clusters[k].numtdtrks);
        m_2ltrkjet_ntdtrk->push_back(mzb.clusters[k].numttdtrks);
     	m_2ltrkjet_ntracks->push_back(mzb.clusters[k].numtracks);
        m_2ltrkjet_phi->push_back(mzb.clusters[k].phi);
        m_2ltrkjet_eta->push_back(mzb.clusters[k].eta);
        m_2ltrkjet_pt->push_back(mzb.clusters[k].pTtot);

        //m_2ltrkjet_p->push_back(jetIter->p());
	}
}
*/
reco::PFJetCollection::const_iterator PFjetIter;
for(PFjetIter=PFJetHandle->begin(); PFjetIter!=PFJetHandle->end(); ++PFjetIter){
        m_PFjet_phi->push_back(PFjetIter->phi());
        m_PFjet_eta->push_back(PFjetIter->eta());
        m_PFjet_pt->push_back(PFjetIter->pt());
        m_PFjet_p->push_back(PFjetIter->p());
	m_PFjet_vz->push_back(PFjetIter->vz());
	//std::cout<<"PFJet status "<<PFjetIter->status()<<std::endl;
}
FillCaloJets(PFCaloJetHandle,CaloTkJetHandle);
eventTree->Fill();

} // end of analyze()
bool L1TrackJetFastProducer::TrackQualityCuts(float trk_pt,int trk_nstub, double trk_chi2){
bool PassQuality=false;
if(trk_nstub==4 && trk_chi2<15)PassQuality=true;
if(trk_nstub==5 && trk_chi2<15 && trk_pt<10)PassQuality=true;
if(trk_nstub==5 && trk_chi2<7 && trk_pt>=10 && trk_pt<40)PassQuality=true;
if(trk_nstub==5 && trk_chi2<5 && trk_pt>40)PassQuality=true;
if(trk_nstub==6 && trk_chi2<5 && trk_pt>50)PassQuality=true;
if(trk_nstub==6 && trk_pt<=50)PassQuality=true;
return PassQuality;

}
void L1TrackJetFastProducer::FillFastJets(std::vector<float>pt, std::vector<float>eta, std::vector<float>phi, std::vector<float>p, std::vector<float>z0, std::vector<int> TruthID, float conesize, std::vector<fastjet::PseudoJet> &JetOutputs_, std::vector<int>&JetNtracks, std::vector<float>&JetVz, std::vector<float>&TrueSumPt){
fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, conesize);
      JetInputs_.clear();
	for(unsigned int j=0; j<pt.size(); ++j){
	 TLorentzVector temp;
         temp.SetPtEtaPhiE(pt[j], eta[j], phi[j], p[j]);	
         fastjet::PseudoJet psuedoJet(temp.Px(), temp.Py(), temp.Pz(), temp.P());
	 JetInputs_.push_back(psuedoJet);
	 JetInputs_.back().set_user_index(j);
	}
	fastjet::ClusterSequence cs(JetInputs_,jet_def);
	JetOutputs_.clear();
	JetOutputs_=fastjet::sorted_by_pt(cs.inclusive_jets(0));
	fjConstituents_.clear();
	for(unsigned int ijet=0;ijet<JetOutputs_.size();++ijet) {
		fjConstituents_ =fastjet::sorted_by_pt(cs.constituents(JetOutputs_[ijet]));		
		float truesumpt=0;
		float sumpt=0;
		float avgZ=0;
		JetNtracks.push_back(fjConstituents_.size());
		for(unsigned int i=0; i<fjConstituents_.size(); ++i){
		    auto index =fjConstituents_[i].user_index();
		    sumpt=sumpt+pt[index];
		    avgZ=avgZ+z0[index]*pt[index];
		    if(TruthID[index]>0)truesumpt=truesumpt+pt[index];	    	
		}
		avgZ=avgZ/sumpt;
		JetVz.push_back(avgZ);
		TrueSumPt.push_back(truesumpt);
	}
}
void L1TrackJetFastProducer::FillCaloJets(Handle<reco::PFJetCollection> PFCaloJetHandle,Handle< L1TkJetParticleCollection> CaloTkJetHandle){
  std::vector<L1TkJetParticle>::const_iterator jetIter;
  m_calotrkjet_eta->clear();
  m_calotrkjet_pt->clear();
  m_calotrkjet_vz->clear();
  m_calotrkjet_phi->clear();
  m_calotrkjet_p->clear();
  float vz=999;
 if(CaloTkJetHandle->begin()!=CaloTkJetHandle->end())vz=CaloTkJetHandle->begin()->getJetVtx();
 for (jetIter = CaloTkJetHandle->begin(); jetIter != CaloTkJetHandle->end(); ++jetIter) {
		  int ibx = jetIter->bx(); // only consider jets from the central BX
      		  if (ibx != 0) continue;
		  if(fabs(vz-jetIter->getJetVtx())>DeltaZ0Cut)continue;
		  m_calotrkjet_p->push_back(jetIter->p());
		  m_calotrkjet_phi->push_back(jetIter->phi());
		  m_calotrkjet_eta->push_back(jetIter->eta());
		  m_calotrkjet_pt->push_back(jetIter->pt());
		  m_calotrkjet_vz->push_back(jetIter->getJetVtx());
     }
  m_calojet_eta->clear();
  m_calojet_pt->clear();
  m_calojet_phi->clear();
  m_calojet_p->clear();
  reco::PFJetCollection::const_iterator PFjetIter=PFCaloJetHandle->begin();
  for(PFjetIter=PFCaloJetHandle->begin(); PFjetIter!=PFCaloJetHandle->end(); ++PFjetIter){
	m_calojet_eta->push_back(PFjetIter->eta());
	m_calojet_phi->push_back(PFjetIter->phi());
	m_calojet_pt->push_back(PFjetIter->pt());
	m_calojet_p->push_back(PFjetIter->p());
	
   }
}
/*
void L1TrackJetFastProducer::L2_cluster(vector< Ptr< L1TTTrackType > > L1TrackPtrs, vector<int>ttrk, vector<int>tdtrk,vector<int>ttdtrk,maxzbin &mzb){
  const int nz = Zbins;
  maxzbin  all_zbins[nz];
  if(all_zbins==NULL) cout<<" all_zbins memory not assigned"<<endl;
  //int best_ind=0;
  //          
  float zmin = -1.0*maxz;
  float zmax = zmin + 2*zstep;
  //                //Create grid of phibins! 
  etaphibin epbins[nphibins][nphibins];
   float phi = -1.0 * M_PI;
  float eta;
  float etamin, etamax, phimin, phimax;
  for(int i = 0; i < nphibins; ++i){
      eta = -1.0 * maxeta;
            for(int j = 0; j < netabins; ++j){
    phimin = phi;
    phimax = phi + phistep;
    etamin = eta;
    eta = eta + etastep;
    etamax = eta;
    epbins[i][j].phi = (phimin + phimax) / 2;
    epbins[i][j].eta = (etamin + etamax) / 2;
       }//for each etabin
       phi = phi + phistep;
   } //for each phibin (finished creating epbins)
  mzb = all_zbins[0];

for(int zbin = 0; zbin < Zbins-1; ++zbin){
  
        //First initialize pT, numtracks, used to 0 (or false)
        for(int i = 0; i < nphibins; ++i){
             for(int j = 0; j < netabins; ++j){
                 epbins[i][j].pTtot = 0;
                 epbins[i][j].used = false;
                 epbins[i][j].numtracks = 0;
                 epbins[i][j].numttrks = 0;
                 epbins[i][j].numtdtrks = 0;
                 epbins[i][j].numttdtrks = 0;
                 }//for each etabin
           } //for each phibin

   for (unsigned int k=0; k<L1TrackPtrs.size(); ++k){
      float trketa=L1TrackPtrs[k]->getMomentum().eta();
      float trkphi=L1TrackPtrs[k]->getMomentum().phi();
      float trkZ=L1TrackPtrs[k]->getPOCA(5).z();
      for(int i = 0; i < nphibins; ++i){
        for(int j = 0; j < netabins; ++j){
          if((zmin <= trkZ && zmax >= trkZ) &&
            ((epbins[i][j].eta - etastep / 2 <= trketa && epbins[i][j].eta + etastep / 2 >= trketa) 
              && epbins[i][j].phi - phistep / 2 <= trkphi && epbins[i][j].phi + phistep / 2 >= trkphi && (zbincount[k] != 2))){
            zbincount.at(k)=zbincount.at(k)+1;
            if(L1TrackPtrs[k]->getMomentum().perp()<PTMAX)epbins[i][j].pTtot += L1TrackPtrs[k]->getMomentum().perp();
	    else epbins[i][j].pTtot +=PTMAX;
            epbins[i][j].numttrks += ttrk[k];
            epbins[i][j].numtdtrks += tdtrk[k];
            epbins[i][j].numttdtrks += ttdtrk[k];
            ++epbins[i][j].numtracks;
    //        cout << epbins[i][j].phi << "\t" << tracks[k].pT << endl;
             } //if right bin
       } //for each phibin: j loop
      }//for each phibin: i loop
     //new 
    }
    etaphibin ** L1clusters = (etaphibin**)malloc(nphibins*sizeof(etaphibin*));

                for(int phislice = 0; phislice < nphibins; ++phislice){
      L1clusters[phislice] = L1_cluster(epbins[phislice]);
      for(int ind = 0; L1clusters[phislice][ind].pTtot != 0; ++ind){
        L1clusters[phislice][ind].used = false;
	//cout<<"L1 Clusters "<< L1clusters[phislice][ind].eta<<", "<< L1clusters[phislice][ind].phi<<", "<< L1clusters[phislice][ind].pTtot<<std::endl;
      }
    }
  //Create clusters array to hold output cluster data for Layer2; can't have more clusters than tracks.
    int ntracks=L1TrackPtrs.size();

    //etaphibin L2cluster[ntracks];//= (etaphibin *)malloc(ntracks * sizeof(etaphibin));
    etaphibin * L2cluster = (etaphibin *)malloc(ntracks * sizeof(etaphibin));
    if(L2cluster==NULL) cout<<"L2cluster memory not assigned"<<endl;

  //Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
    float hipT = 0;
    int nclust = 0;
    int phibin = 0;
    int imax=-1;
       //index of clusters array for each phislice.
    int index1;
    float E1 =0;
    float E0 =0;
    float E2 =0;
    int trx1, trx2;
    int ttrk1, ttrk2;
    int tdtrk1, tdtrk2;
    int ttdtrk1, ttdtrk2;
    int used1, used2, used3, used4;

      //Find eta-phibin with highest pT.
    for(phibin = 0; phibin < nphibins; ++phibin){
        while(true){
      hipT = 0;
      for(index1 = 0; L1clusters[phibin][index1].pTtot > 0; ++index1){
        if(!L1clusters[phibin][index1].used && L1clusters[phibin][index1].pTtot >= hipT){
          hipT = L1clusters[phibin][index1].pTtot;
          imax = index1;
        }
      }//for each index within the phibin
          //If highest pT is 0, all bins are used.
      if(hipT == 0){
        break;
      }
      E0 = hipT;   //E0 is pT of first phibin of the cluster.
      E1 = 0;
      E2 = 0;
      trx1 = 0;
      trx2 = 0;
      ttrk1 = 0;
      ttrk2 = 0;
      tdtrk1 = 0;
      tdtrk2 = 0;
      ttdtrk1 = 0;
      ttdtrk2 = 0;
      L2cluster[nclust] = L1clusters[phibin][imax];
      L1clusters[phibin][imax].used = true;
    //Add pT of upper neighbor.
    //E1 is pT of the middle phibin (should be highest pT)
      if(phibin != nphibins-1){
        used1 = -1;
        used2 = -1;
        for (index1 = 0; L1clusters[phibin+1][index1].pTtot != 0; ++index1){
          if(L1clusters[phibin+1][index1].used){
            continue;
          }
          if(fabs(L1clusters[phibin+1][index1].eta - L1clusters[phibin][imax].eta) <= 1.5*etastep){
            E1 += L1clusters[phibin+1][index1].pTtot;
            trx1 += L1clusters[phibin+1][index1].numtracks;
            ttrk1 += L1clusters[phibin+1][index1].numttrks;
            tdtrk1 += L1clusters[phibin+1][index1].numtdtrks;
            ttdtrk1 += L1clusters[phibin+1][index1].numttdtrks;
            if(used1 < 0)
              used1 = index1;
            else
              used2 = index1;
          }//if cluster is within one phibin
        } //for each cluster in above phibin
      //if E1 isn't higher, E0 and E1 are their own cluster.
        if(E1 < E0){
          L2cluster[nclust].pTtot += E1;   
          L2cluster[nclust].numtracks += trx1;
          L2cluster[nclust].numttrks += ttrk1;
          L2cluster[nclust].numtdtrks += tdtrk1;
          L2cluster[nclust].numttdtrks += ttdtrk1;
          if(used1 >= 0)
            L1clusters[phibin+1][used1].used = true;
          if(used2 >= 0)
            L1clusters[phibin+1][used2].used = true;
          ++nclust;
          continue;
        }
        
        if(phibin != nphibins-2){
                                      //E2 will be the pT of the third phibin (should be lower than E1).
          used3 = -1;
          used4 = -1;
          for (index1 = 0; L1clusters[phibin+2][index1].pTtot != 0; ++index1){
            if(L1clusters[phibin+2][index1].used){
              continue;
            }
            if(fabs(L1clusters[phibin+2][index1].eta - L1clusters[phibin][imax].eta) <= 1.5*etastep){
              E2 += L1clusters[phibin+2][index1].pTtot;
              trx2 += L1clusters[phibin+2][index1].numtracks;
              ttrk2 += L1clusters[phibin+2][index1].numttrks;
              tdtrk2 += L1clusters[phibin+2][index1].numtdtrks;
              ttdtrk2 += L1clusters[phibin+2][index1].numttdtrks;
              if(used3 < 0)
                used3 = index1;
              else
                used4 = index1;
            }
    
          }
             //if indeed E2 < E1, add E1 and E2 to E0, they're all a cluster together.
             //  otherwise, E0 is its own cluster.
          if(E2 < E1){
            L2cluster[nclust].pTtot += E1 + E2;
            L2cluster[nclust].numtracks += trx1 + trx2;
            L2cluster[nclust].numttrks += ttrk1 + ttrk2;
            L2cluster[nclust].numtdtrks += tdtrk1 + tdtrk2;
            L2cluster[nclust].numttdtrks += ttdtrk1 + ttdtrk2;
            L2cluster[nclust].phi = L1clusters[phibin+1][used1].phi;  
            if(used1 >= 0)
              L1clusters[phibin+1][used1].used = true;
            if(used2 >= 0)
              L1clusters[phibin+1][used2].used = true;
            if(used3 >= 0)
              L1clusters[phibin+2][used3].used = true;
            if(used4 >= 0)
              L1clusters[phibin+2][used4].used = true;
          }
          ++nclust;
          continue;
        } // end Not nphibins-2
        else{
          L2cluster[nclust].pTtot += E1;
          L2cluster[nclust].numtracks += trx1;
          L2cluster[nclust].numttrks += ttrk1;
          L2cluster[nclust].numtdtrks += tdtrk1;
          L2cluster[nclust].numttdtrks += ttdtrk1;
          L2cluster[nclust].phi = L1clusters[phibin+1][used1].phi;
          if(used1 >= 0)
            L1clusters[phibin+1][used1].used = true;
          if(used2 >= 0)
            L1clusters[phibin+1][used2].used = true;
          ++nclust;
          continue;
        }
      }//End not last phibin(23)
      else { //if it is phibin 23
        L1clusters[phibin][imax].used = true;
        ++nclust;
      }
        }//while hipT not 0
    }//for each phibin
    //for(int db=0;db<nclust;++db)cout<<L2cluster[db].phi<<"\t"<<L2cluster[db].pTtot<<"\t"<<L2cluster[db].numtracks<<endl;  
  //Now merge clusters, if necessary
 for(int m = 0; m < nclust -1; ++m){
                     for(int n = m+1; n < nclust; ++n)
                        if(L2cluster[n].eta == L2cluster[m].eta && (fabs(L2cluster[n].phi - L2cluster[m].phi) < 1.5*phistep || fabs(L2cluster[n].phi - L2cluster[m].phi) > 6.0)){
                                if(L2cluster[n].pTtot > L2cluster[m].pTtot){
                                        L2cluster[m].phi = L2cluster[n].phi;
                                }
                                L2cluster[m].pTtot += L2cluster[n].pTtot;
                                L2cluster[m].numtracks += L2cluster[n].numtracks;
        L2cluster[m].numttrks += L2cluster[n].numttrks;
        L2cluster[m].numtdtrks += L2cluster[n].numtdtrks;
        L2cluster[m].numttdtrks += L2cluster[n].numttdtrks;
                                for(int m1 = n; m1 < nclust-1; ++m1){
                                        L2cluster[m1] = L2cluster[m1+1];
                                }
                                nclust--;
                                m = -1;
                                break; //?????
                        }//end if clusters neighbor in eta
                }//end for (m) loop     
          //sum up all pTs in this zbin to find ht.
    float ht = 0;
    for(int k = 0; k < nclust; ++k){
                        if(L2cluster[k].pTtot>50 && L2cluster[k].numtracks<2)continue;
                        if(L2cluster[k].pTtot>100 && L2cluster[k].numtracks<=4)continue;
                        if(L2cluster[k].pTtot>5){
      			ht += L2cluster[k].pTtot;
                }
	}
     //if ht is larger than previous max, this is the new vertex zbin.
      //all_zbins[zbin].mcd = mcd;
    all_zbins[zbin].znum = zbin;
    all_zbins[zbin].clusters = (etaphibin *)malloc(nclust*sizeof(etaphibin));
    all_zbins[zbin].nclust = nclust;
    all_zbins[zbin].zbincenter=(zmin+zmax)/2.0;
    for(int k = 0; k < nclust; ++k){
      all_zbins[zbin].clusters[k].phi = L2cluster[k].phi;                               
      all_zbins[zbin].clusters[k].eta = L2cluster[k].eta;                             
      all_zbins[zbin].clusters[k].pTtot = L2cluster[k].pTtot;
      all_zbins[zbin].clusters[k].numtracks = L2cluster[k].numtracks;
      all_zbins[zbin].clusters[k].numttrks = L2cluster[k].numttrks;
      all_zbins[zbin].clusters[k].numtdtrks = L2cluster[k].numtdtrks;
      all_zbins[zbin].clusters[k].numttdtrks = L2cluster[k].numttdtrks;
    }
  //  for(int db=0;db<nclust;++db)cout<<all_zbins[zbin].clusters[db].phi<<"\t"<<all_zbins[zbin].clusters[db].pTtot<<endl; 
    all_zbins[zbin].ht = ht;
    if(ht >= mzb.ht){
      mzb = all_zbins[zbin];
      //mzb.zbincenter=(zmin+zmax)/2.0;
     // best_ind=zbin;
    }
    //Prepare for next zbin!
    zmin = zmin + zstep;
    zmax = zmax + zstep;
      
    //   for(int phislice = 0; phislice < nphibins; ++phislice){
    //   free(L1clusters[phislice]);      
    // }
    free(L1clusters);
    //free(all_zbins);
    // free(L2cluster);


    } //for each zbin
   std::cout<<"Chosen Z -bin "<<mzb.zbincenter<<std::endl; 
    for(int k = 0; k < mzb.nclust; ++k)std::cout<<"L2 Eta, Phi "<<mzb.clusters[k].eta<<", "<<mzb.clusters[k].phi<<", "<<mzb.clusters[k].pTtot<<std::endl;
}

etaphibin * L1TrackJetFastProducer::L1_cluster(etaphibin *phislice){

    etaphibin * clusters = (etaphibin *)malloc(netabins/2 * sizeof(etaphibin));
    if(clusters==NULL) cout<<"clusters memory not assigned"<<endl;
  //Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
    float my_pt, left_pt, right_pt, right2pt;

    int nclust = 0;
    right2pt=0;
    for(int etabin = 0; etabin < netabins; ++etabin){
      //assign values for my pT and neighbors' pT
      if(phislice[etabin].used) continue;
      my_pt = phislice[etabin].pTtot;
      if(etabin > 0 && !phislice[etabin-1].used) {
        left_pt = phislice[etabin-1].pTtot;
        // if(etabin > 1 && !phislice[etabin-2].used) {
        //   left2pt = phislice[etabin-2].pTtot;
        // } else left2pt = 0;
      } else left_pt = 0;
      if(etabin < netabins - 1 && !phislice[etabin+1].used) {
        right_pt = phislice[etabin+1].pTtot;
        if(etabin < netabins - 2 && !phislice[etabin+2].used) {
          right2pt = phislice[etabin+2].pTtot;
        } else right2pt = 0;
      } else right_pt = 0;
    
    //if I'm not a cluster, move on.
      if(my_pt < left_pt || my_pt <= right_pt) {
         //if unused pT in the left neighbor, spit it out as a cluster.
              if(left_pt > 0) {
          clusters[nclust] = phislice[etabin-1];
          phislice[etabin-1].used = true;
          ++nclust;
        }
        continue;
      }

    //I guess I'm a cluster-- should I use my right neighbor?
    // Note: left neighbor will definitely be used because if it 
    //       didn't belong to me it would have been used already
      clusters[nclust] = phislice[etabin];
      phislice[etabin].used = true;
      if(left_pt > 0) {
        clusters[nclust].pTtot += left_pt;
        clusters[nclust].numtracks += phislice[etabin-1].numtracks;
        clusters[nclust].numttrks += phislice[etabin-1].numttrks;
        clusters[nclust].numtdtrks += phislice[etabin-1].numtdtrks;
        clusters[nclust].numttdtrks += phislice[etabin-1].numttdtrks;
      }
      if(my_pt >= right2pt && right_pt > 0) {
        clusters[nclust].pTtot += right_pt;
        clusters[nclust].numtracks += phislice[etabin+1].numtracks;
        clusters[nclust].numttrks += phislice[etabin+1].numttrks;
        clusters[nclust].numtdtrks += phislice[etabin+1].numtdtrks;
        clusters[nclust].numttdtrks += phislice[etabin+1].numttdtrks;
        phislice[etabin+1].used = true;
      }

      ++nclust;
    } //for each etabin                       
                           
  //Now merge clusters, if necessary
    for(int m = 0; m < nclust -1; ++m){
      if(fabs(clusters[m+1].eta - clusters[m].eta) < 1.5*etastep){
        if(clusters[m+1].pTtot > clusters[m].pTtot){
          clusters[m].eta = clusters[m+1].eta;
        }
        clusters[m].pTtot += clusters[m+1].pTtot;
        clusters[m].numtracks += clusters[m+1].numtracks;  //Previous version didn't add tracks when merging. 
        clusters[m].numttrks += clusters[m+1].numttrks;
        clusters[m].numtdtrks += clusters[m+1].numtdtrks;
        clusters[m].numttdtrks += clusters[m+1].numttdtrks;
        for(int m1 = m+1; m1 < nclust-1; ++m1){
          clusters[m1] = clusters[m1+1];
        }
        nclust--;
        m = -1;
      }//end if clusters neighbor in eta
    }//end for (m) loop
//  for(int i = 0; i < nclust; ++i) cout << clusters[i].phi << "\t" << clusters[i].pTtot << "\t" << clusters[i].numtracks << endl;
  //zero out remaining unused clusters.
  for(int i = nclust; i < netabins/2; ++i){
    clusters[i].pTtot = 0;
  }
  return clusters;
}
*/
///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackJetFastProducer);
