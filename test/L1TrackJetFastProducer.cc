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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "DataFormats/L1TVertex/interface/Vertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

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
using namespace std;
using namespace edm;


//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TrackJetFastProducer : public edm::EDAnalyzer
{
public:

  // Constructor/destructor
  explicit L1TrackJetFastProducer(const edm::ParameterSet& iConfig);
  virtual ~L1TrackJetFastProducer();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
 // virtual void FillJets(std::vector<fastjet::PseudoJet>  JetInputs_, bool Prompt, bool TrueTP, edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle);  
  virtual void FillFastJets(std::vector<float>pt, std::vector<float>eta, std::vector<float>phi, std::vector<float>p, std::vector<float>z0, float conesize, std::vector<fastjet::PseudoJet> &JetOutputs_, std::vector<fastjet::PseudoJet> &fjConstituents, std::vector<float>&JetVz);

protected:
  
private:
  
  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config; 
  
  int MyProcess;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;       // lots of debug printout statements
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  int L1Tk_nPar;        // use 4 or 5 parameter track fit? 
  int TP_minNStub;      // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer; // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;      // save TPs with pt > minPt 
  double TP_maxEta;     // save TPs with |eta| < maxEta 
  double TP_maxZ0;      // save TPs with |z0| < maxZ0 
  double CHI2MAX;       // save Chi2 per ndof <5
  double DeltaZ0Cut;    // save with |L1z-z0| < maxZ0
  double CONESize;      // Use anti-kt with this cone size
  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)
  
  edm::InputTag L1TrackInputTag;        // L1 track collection
  edm::InputTag MCTruthTrackInputTag;   // MC truth collection
  edm::InputTag L1StubInputTag;
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;
  
  edm::InputTag TrueVertexInputTag;
  edm::InputTag MCVertexInputTag;
  edm::InputTag RecoVertexInputTag;
  edm::InputTag  GenParticleInputTag;
  edm::InputTag GenJetAK4;

  edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;

  edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  edm::EDGetTokenT< vector<reco::GenParticle> > HEPMCVertexToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet> >GenJetCollectionToken_;
  edm::EDGetTokenT<L1TkPrimaryVertexCollection>TPVertexToken_;
  edm::EDGetTokenT<L1TkPrimaryVertexCollection>MCVertexToken_;
  edm::EDGetTokenT< l1t::VertexCollection>L1VertexToken_;

  std::vector<fastjet::PseudoJet>  RecoJetInputs_;
  std::vector<fastjet::PseudoJet>  GenJetInputs_;
  std::vector<fastjet::PseudoJet>  JetInputs_;
  std::vector<float>JetVz_;
  std::vector<fastjet::PseudoJet> JetOutputs_;
  std::vector<fastjet::PseudoJet> fjConstituents_;
  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;
  // primary vertex
  std::vector<float>* m_pv_L1recofakesumpt;
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
  std::vector<float>* m_matchtrk_p;
  std::vector<float>* m_matchtrk_pt;
  std::vector<float>* m_matchtrk_eta;
  std::vector<float>* m_matchtrk_phi;
  std::vector<float>* m_matchtrk_d0; //this variable is only filled if L1Tk_nPar==5
  std::vector<float>* m_matchtrk_z0;
  std::vector<float>* m_matchtrk_chi2; 
  std::vector<int>*   m_matchtrk_nstub;

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

  std::vector<float>* m_recojet_vz;
  std::vector<float>* m_recojet_p;
  std::vector<float>* m_recojet_phi;
  std::vector<float>* m_recojet_eta;
  std::vector<int>* m_recojet_ntracks;
  std::vector<float>* m_recojet_tp_sumpt;
  std::vector<float>* m_recojet_truetp_sumpt;
  std::vector<float>* m_recojet_pt;

  std::vector<float>* m_tpjet_vz;
  std::vector<float>* m_tpjet_p;
  std::vector<float>* m_tpjet_phi;
  std::vector<float>* m_tpjet_eta;
  std::vector<int>* m_tpjet_ntracks;
  std::vector<float>* m_tpjet_tp_sumpt;
  std::vector<float>* m_tpjet_truetp_sumpt;
  std::vector<float>* m_tpjet_pt;
  std::vector<float>* m_jet_matchtrk_sumpt;
  std::vector<float>* m_jet_loosematchtrk_sumpt;
  std::vector<float>* m_jet_trk_sumpt;
};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1TrackJetFastProducer::L1TrackJetFastProducer(edm::ParameterSet const& iConfig) : 
  config(iConfig)
{

  MyProcess        = iConfig.getParameter< int >("MyProcess");
  DebugMode        = iConfig.getParameter< bool >("DebugMode");
  SaveAllTracks    = iConfig.getParameter< bool >("SaveAllTracks");
  L1Tk_nPar        = iConfig.getParameter< int >("L1Tk_nPar");
  //Cuts on Tracks
  TP_minNStub      = iConfig.getParameter< int >("TP_minNStub");
  TP_minNStubLayer = iConfig.getParameter< int >("TP_minNStubLayer");
  TP_minPt         = iConfig.getParameter< double >("PTMINTRA");
  TP_maxEta        = 2.4;
  TP_maxZ0         = iConfig.getParameter< double >("ZMAX");
  DeltaZ0Cut        = iConfig.getParameter< double >("DeltaZ0Cut");
  CHI2MAX 	    = iConfig.getParameter< double >("CHI2MAX");
  CONESize	    =iConfig.getParameter<double>("CONESize");
  L1TrackInputTag      = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  MCVertexInputTag    = iConfig.getParameter<edm::InputTag>("MCVertexInputTag");
  TrueVertexInputTag    = iConfig.getParameter<edm::InputTag>("TrueVertexInputTag");
  RecoVertexInputTag    = iConfig.getParameter<edm::InputTag>("RecoVertexInputTag");
  GenParticleInputTag   = iConfig.getParameter<edm::InputTag >("GenParticleInputTag");
  MCTruthStubInputTag = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  GenJetAK4=iConfig.getParameter<edm::InputTag>("GenJetAK4");
  HEPMCVertexToken_=consumes< std::vector< reco::GenParticle> >(GenParticleInputTag);
  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);
  ttTrackMCTruthToken_ = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthTrackInputTag);
  ttStubMCTruthToken_ = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);
  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);
  GenJetCollectionToken_=consumes< std::vector<reco::GenJet > >(GenJetAK4);
  MCVertexToken_=consumes<L1TkPrimaryVertexCollection>(MCVertexInputTag);
  TPVertexToken_=consumes<L1TkPrimaryVertexCollection>(TrueVertexInputTag);
  L1VertexToken_=consumes<l1t::VertexCollection>(RecoVertexInputTag);

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
  edm::Service<TFileService> fs;


  // initilize
  m_trk_p    = new std::vector<float>;
  m_trk_pt    = new std::vector<float>;
  m_trk_eta   = new std::vector<float>;
  m_trk_phi   = new std::vector<float>;
  m_trk_z0    = new std::vector<float>;
  m_trk_d0    = new std::vector<float>;
  m_trk_chi2  = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
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

  m_matchtrk_p    = new std::vector<float>;
  m_matchtrk_pt    = new std::vector<float>;
  m_matchtrk_eta   = new std::vector<float>;
  m_matchtrk_phi   = new std::vector<float>;
  m_matchtrk_z0    = new std::vector<float>;
  m_matchtrk_d0    = new std::vector<float>;
  m_matchtrk_chi2  = new std::vector<float>;
  m_matchtrk_nstub = new std::vector<int>;
  
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

  m_recojet_eta = new std::vector<float>;
  m_recojet_vz = new std::vector<float>;
  m_recojet_phi = new std::vector<float>;
  m_recojet_p = new std::vector<float>;
  m_recojet_pt = new std::vector<float>;
  m_recojet_ntracks = new std::vector<int>;
  m_recojet_tp_sumpt = new std::vector<float>;
  m_recojet_truetp_sumpt = new std::vector<float>;


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
    eventTree->Branch("trk_pt",    &m_trk_pt);
    eventTree->Branch("trk_eta",   &m_trk_eta);
    eventTree->Branch("trk_phi",   &m_trk_phi);
    eventTree->Branch("trk_d0",    &m_trk_d0);
    eventTree->Branch("trk_z0",    &m_trk_z0);
    eventTree->Branch("trk_chi2",  &m_trk_chi2);
    eventTree->Branch("trk_nstub", &m_trk_nstub);

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

  eventTree->Branch("matchtrk_p",      &m_matchtrk_p);
  eventTree->Branch("matchtrk_pt",      &m_matchtrk_pt);
  eventTree->Branch("matchtrk_eta",     &m_matchtrk_eta);
  eventTree->Branch("matchtrk_phi",     &m_matchtrk_phi);
  eventTree->Branch("matchtrk_z0",      &m_matchtrk_z0);
  eventTree->Branch("matchtrk_d0",      &m_matchtrk_d0);
  eventTree->Branch("matchtrk_chi2",    &m_matchtrk_chi2);
  eventTree->Branch("matchtrk_nstub",   &m_matchtrk_nstub);
  }

    eventTree->Branch("pv_L1recofakesumpt", &m_pv_L1recofakesumpt);
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
    eventTree->Branch("recojet_eta", &m_recojet_eta);
    eventTree->Branch("recojet_vz", &m_recojet_vz);
    eventTree->Branch("recojet_p", &m_recojet_p);
    eventTree->Branch("recojet_pt", &m_recojet_pt);
    eventTree->Branch("recojet_phi", &m_recojet_phi);
    eventTree->Branch("recojet_ntracks", &m_recojet_ntracks);
    eventTree->Branch("recojet_truetp_sumpt", m_recojet_truetp_sumpt);
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
void L1TrackJetFastProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  m_matchtrk_p->clear();
  m_matchtrk_pt->clear();
  m_matchtrk_eta->clear();
  m_matchtrk_phi->clear();
  m_matchtrk_z0->clear();
  m_matchtrk_d0->clear();
  m_matchtrk_chi2->clear();
  m_matchtrk_nstub->clear();

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

  m_recojet_eta->clear();
  m_recojet_pt->clear();
  m_recojet_vz->clear();
  m_recojet_phi->clear();
  m_recojet_p->clear();
  m_recojet_ntracks->clear();
  m_recojet_truetp_sumpt->clear();
  m_recojet_tp_sumpt->clear();
  m_tpjet_eta->clear();
  m_tpjet_pt->clear();
  m_tpjet_vz->clear();
  m_tpjet_phi->clear();
  m_tpjet_p->clear();
  m_tpjet_ntracks->clear();
  m_tpjet_tp_sumpt->clear();
  m_tpjet_truetp_sumpt->clear();
  m_pv_L1recofakesumpt->clear();
  m_pv_L1recotruesumpt->clear();
  m_pv_L1recosumpt->clear();
  m_pv_L1reco->clear();
  m_pv_L1TPsumpt->clear();
  m_pv_L1TP->clear();
  m_pv_MC->clear();
  m_MC_lep->clear();
  m_pv_MCChgSumpT->clear();
  // -----------------------------------------------------------------------------------------------
  // retrieve various containers
  // -----------------------------------------------------------------------------------------------

  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);
    edm::Handle<std::vector<reco::GenJet> >GenJetsAK4Handle;
  iEvent.getByToken(GenJetCollectionToken_,GenJetsAK4Handle); 


   edm::Handle< std::vector< reco::GenParticle> > GenParticleHandle;
   iEvent.getByToken(HEPMCVertexToken_,GenParticleHandle);
        vector<reco::GenParticle>::const_iterator genpartIter ;
    	int leptonicCount=0;
        for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
		if (abs(genpartIter ->pdgId())!=24)continue;
		//std::cout<<"W-boson mother "<<genpartIter ->mother(0)->pdgId()<<std::endl;

		if(abs(genpartIter ->mother(0)->pdgId())!=6 && abs(genpartIter ->mother(0)->pdgId())!=24 )continue;
		if( abs(genpartIter ->daughter(0)->pdgId())==11  ||  abs(genpartIter ->daughter(0)->pdgId())==13 ||  abs(genpartIter ->daughter(0)->pdgId())==15){++leptonicCount;}
	}
  m_MC_lep->push_back(leptonicCount);
  // MC truth association maps
   edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
   iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);

  edm::Handle<L1TkPrimaryVertexCollection >L1TkPrimaryVertexTPHandle;
  iEvent.getByToken(TPVertexToken_, L1TkPrimaryVertexTPHandle);

  edm::Handle<L1TkPrimaryVertexCollection >L1TkPrimaryVertexMCHandle;
  iEvent.getByToken(MCVertexToken_, L1TkPrimaryVertexMCHandle);

  if(L1TkPrimaryVertexMCHandle.isValid()){
	m_pv_MC->push_back(L1TkPrimaryVertexMCHandle->begin()->getZvertex());  
	m_pv_MCChgSumpT->push_back(L1TkPrimaryVertexMCHandle->begin()->getSum());
  }
  if(L1TkPrimaryVertexTPHandle.isValid()){
	m_pv_L1TPsumpt->push_back(L1TkPrimaryVertexTPHandle->begin()->getSum());
	m_pv_L1TP->push_back(L1TkPrimaryVertexTPHandle->begin()->getZvertex());
  }
  edm::Handle< l1t::VertexCollection >L1TkPrimaryVertexHandle;
  iEvent.getByToken(L1VertexToken_, L1TkPrimaryVertexHandle);
 if(L1TkPrimaryVertexHandle.isValid()){
	m_pv_L1reco->push_back(L1TkPrimaryVertexHandle->begin()->z0());
	//add Vertex True quality, Vertex Fake Content, total sumpT, and number of tracks 
	std::vector<edm::Ptr< TTTrack<Ref_Phase2TrackerDigi_> > >Vtxtracks=L1TkPrimaryVertexHandle->begin()->tracks();
	//std::cout<<"Ntracks in Vertex "<<Vtxtracks.size()<<std::endl;
	float sumpt=0;
	float trueContent=0;
	float fakeContent=0;
	for(unsigned int t=0; 	t<Vtxtracks.size(); ++t){
		sumpt=Vtxtracks[t]->getMomentum(L1Tk_nPar).perp()+sumpt;
		//if(Vtxtracks[t].isAvailable())std::cout<<"Has a collection in memory "<<Vtxtracks[t].id()<<std::endl;
		//edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(Vtxtracks[t]);	
		//if(MCTruthTTTrackHandle->isGenuine(l1track_ptr))std::cout<<"MC Handle is valid "<<std::endl;//trueContent=trueContent+Vtxtracks[t]->getMomentum(L1Tk_nPar).perp();
		if(MCTruthTTTrackHandle->isGenuine(Vtxtracks[t]))trueContent=trueContent+Vtxtracks[t]->getMomentum(L1Tk_nPar).perp();
		else fakeContent=fakeContent+Vtxtracks[t]->getMomentum(L1Tk_nPar).perp();
	}
	//std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track=Vtxtracks->begin();;		
	//std::cout<<"Vtx Sum pT "<<trueContent<<std::endl;
	m_pv_L1recofakesumpt->push_back(fakeContent);
	m_pv_L1recotruesumpt->push_back(trueContent);
	m_pv_L1recosumpt->push_back(sumpt);	
}
  // -----------------------------------------------------------------------------------------------
  // more for TTStubs
  edm::ESHandle<TrackerGeometry> geometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

  edm::ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  //const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();
//Gen Jet branches
std::vector<float>genConstpt;
std::vector<float>genConsteta;
std::vector<float>genConstphi;
std::vector<float>genConstp;
std::vector<float>genConstz;
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
	}
		else NeuEnergy=NeuEnergy+genJet.getGenConstituent(g)->energy();
        	}
		m_genak4jet_chgfrac->push_back(ChgEnergy);
	        m_genak4jet_neufrac->push_back(NeuEnergy);
    }
  }
//fill chg particle fastjets
JetOutputs_.clear();
fjConstituents_.clear();
JetVz_.clear();
FillFastJets(genConstpt, genConsteta, genConstphi, genConstp, genConstz, CONESize, JetOutputs_, fjConstituents_, JetVz_);
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
	RecoJetInputs_.clear();

    for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);

      this_l1track++;
      float tmp_trk_p   = iterL1Track->getMomentum(L1Tk_nPar).mag();
      float tmp_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
      float tmp_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
      float tmp_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
      float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm
      
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

      int tmp_trk_genuine = 0;
      int tmp_trk_loose = 0;
      int tmp_trk_unknown = 0;
      int tmp_trk_combinatoric = 0;
      if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
      if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
      if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
      if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;
      
      if (DebugMode) {
	cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi 
	     << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
	if (tmp_trk_genuine) cout << " (is genuine)" << endl; 
	if (tmp_trk_unknown) cout << " (is unknown)" << endl; 
	if (tmp_trk_combinatoric) cout << " (is combinatoric)" << endl; 
      }
      if(tmp_trk_pt>TP_minPt && fabs(m_pv_L1reco->at(0)-tmp_trk_z0)<DeltaZ0Cut && (tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar)<CHI2MAX  || tmp_trk_pt<20) && tmp_trk_nstub>=TP_minNStub ){
      m_trk_p ->push_back(tmp_trk_p); 
      m_trk_pt ->push_back(tmp_trk_pt);
      m_trk_eta->push_back(tmp_trk_eta);
      m_trk_phi->push_back(tmp_trk_phi);
      m_trk_z0 ->push_back(tmp_trk_z0);
      if (L1Tk_nPar==5) m_trk_d0->push_back(tmp_trk_d0);
      else m_trk_d0->push_back(999.);
      m_trk_chi2 ->push_back(tmp_trk_chi2/(2*tmp_trk_nstub - L1Tk_nPar));
      m_trk_nstub->push_back(tmp_trk_nstub);
      m_trk_genuine->push_back(tmp_trk_genuine);
      m_trk_loose->push_back(tmp_trk_loose);
      m_trk_unknown->push_back(tmp_trk_unknown);
      m_trk_combinatoric->push_back(tmp_trk_combinatoric);
      }

      // ----------------------------------------------------------------------------------------------
      // for studying the fake rate
      // ----------------------------------------------------------------------------------------------

      edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

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

  }//end if SaveAllTracks

  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------


  if (DebugMode) cout << endl << "Loop over tracking particles!" << endl;
  int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;

  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
 
    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
   ++this_tp;
    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
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

    std::vector< edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > > matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);
    
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
	    edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
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

	edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
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

    if( fabs(tmp_tp_eta)<TP_maxEta && fabs(m_pv_L1TP->at(0)-tmp_tp_z0)<DeltaZ0Cut && nStubLayerTP>=TP_minNStubLayer && tmp_tp_pt>TP_minPt){
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
    m_tp_eventid->push_back(tmp_eventid);
    m_tp_charge->push_back(tmp_tp_charge);

    m_matchtrk_p ->push_back(tmp_matchtrk_p);
    m_matchtrk_pt ->push_back(tmp_matchtrk_pt);
    m_matchtrk_eta->push_back(tmp_matchtrk_eta);
    m_matchtrk_phi->push_back(tmp_matchtrk_phi);
    m_matchtrk_z0 ->push_back(tmp_matchtrk_z0);
    m_matchtrk_d0 ->push_back(tmp_matchtrk_d0);
    m_matchtrk_chi2 ->push_back(tmp_matchtrk_chi2);
    m_matchtrk_nstub->push_back(tmp_matchtrk_nstub);
     }
  } //end loop tracking particles
JetOutputs_.clear();
fjConstituents_.clear();
JetVz_.clear();
FillFastJets(*m_tp_pt, *m_tp_eta, *m_tp_phi, *m_tp_p, *m_tp_z0, CONESize, JetOutputs_, fjConstituents_, JetVz_);
for(unsigned int j=0; j<JetOutputs_.size(); ++j){
	TLorentzVector temp;
	temp.SetPtEtaPhiE(JetOutputs_[j].pt(),JetOutputs_[j].eta(),JetOutputs_[j].phi_std(),JetOutputs_[j].modp());
        m_tpjet_pt->push_back(temp.Pt());
        m_tpjet_p->push_back(temp.P());
        m_tpjet_eta->push_back(temp.Eta());
        m_tpjet_phi->push_back(temp.Phi());
        m_tpjet_vz->push_back(JetVz_[j]);
	m_tpjet_ntracks->push_back(fjConstituents_.size());
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

FillFastJets(*m_trk_pt, *m_trk_eta, *m_trk_phi, *m_trk_p, *m_trk_z0, CONESize, JetOutputs_, fjConstituents_, JetVz_);
for(unsigned int j=0; j<JetOutputs_.size(); ++j){
        TLorentzVector temp;
        temp.SetPtEtaPhiE(JetOutputs_[j].pt(),JetOutputs_[j].eta(),JetOutputs_[j].phi_std(),JetOutputs_[j].modp());
	m_recojet_vz->push_back(JetVz_[j]);
	m_recojet_ntracks->push_back(fjConstituents_.size());
  	m_recojet_phi->push_back(JetOutputs_[j].phi_std());
  	m_recojet_eta->push_back(JetOutputs_[j].eta());
  	m_recojet_pt->push_back(JetOutputs_[j].pt());
  	m_recojet_p->push_back(JetOutputs_[j].modp());
        float sumpt=0;
        float truesumpt=0;
	for(unsigned int i=0; i<fjConstituents_.size(); ++i){
        	auto index =fjConstituents_[i].user_index();	
		int eventId=m_trk_genuine->at(index);
		if(eventId>0)truesumpt=truesumpt+m_trk_pt->at(index);
		sumpt=sumpt+m_trk_pt->at(index);
	}
	m_recojet_tp_sumpt->push_back(sumpt);
	m_recojet_truetp_sumpt->push_back(truesumpt);
}
JetOutputs_.clear();
fjConstituents_.clear();
JetVz_.clear();

//FillFastJets(*m_matchtrk_pt, *m_matchtrk_eta, *m_matchtrk_phi, *m_matchtrk_p, *m_matchtrk_z0, CONESize, JetOutputs_, fjConstituents_, JetVz_);
eventTree->Fill();

} // end of analyze()

void L1TrackJetFastProducer::FillFastJets(std::vector<float>pt, std::vector<float>eta, std::vector<float>phi, std::vector<float>p, std::vector<float>z0, float conesize, std::vector<fastjet::PseudoJet> &JetOutputs_, std::vector<fastjet::PseudoJet> &fjConstituents, std::vector<float>&JetVz){
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
	for(unsigned int ijet=0;ijet<JetOutputs_.size();++ijet) {
		std::vector<fastjet::PseudoJet> fjConstituents =fastjet::sorted_by_pt(cs.constituents(JetOutputs_[ijet]));		
		float sumpt=0;
		float avgZ=0;
		for(unsigned int i=0; i<fjConstituents.size(); ++i){
		    auto index =fjConstituents[i].user_index();
		    sumpt=sumpt+pt[index];
		    avgZ=avgZ+z0[index]*pt[index];
		    	
		}
		avgZ=avgZ/sumpt;
		JetVz.push_back(avgZ);
	}
}
///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackJetFastProducer);
