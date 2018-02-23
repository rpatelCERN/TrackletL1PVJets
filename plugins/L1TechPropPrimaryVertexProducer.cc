// -*- C++ -*-
//
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Tue Nov 12 17:03:19 CET 2013
// $Id$
//
//

// -------------------------------------------------------------------------------------------------------
//
//	********  OLD CODE   ********
//
//	********  The latest producer for the primary vertex is  L1TkFastVertexProducer.cc      ********
//
// --------------------------------------------------------------------------------------------------------



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
#include "DataFormats/L1TVertex/interface/Vertex.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "TH1F.h"

//
// class declaration
//

class L1TechPropPrimaryVertexProducer : public edm::EDProducer {
   public:

   typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
   typedef std::vector< L1TTTrackType >  L1TTTrackCollectionType;

      explicit L1TechPropPrimaryVertexProducer(const edm::ParameterSet&);
      ~L1TechPropPrimaryVertexProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      void TrackingParticleVertex(edm::Event& iEvent,edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle, edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle, const TrackerTopology* topol);
      float MaxPtVertex(const edm::Handle<L1TTTrackCollectionType> & L1TTTrackHandle,
                float& sum,
                int nmin, int nPSmin, float ptmin, int imode,
		const TrackerTopology* topol) ;

      float SumPtVertex(const edm::Handle<L1TTTrackCollectionType> & L1TTTrackHandle,
                float z, int nmin, int nPSmin, float ptmin, int imode,
		const TrackerTopology* topol, std::vector< edm::Ptr< L1TTTrackType > >  &l1Tracks);
      float WeightedZVertex(const edm::Handle<L1TTTrackCollectionType> & L1TTTrackHandle,
                float z, int nmin, int nPSmin, float ptmin, int imode,
		const TrackerTopology* topol, std::vector< edm::Ptr< L1TTTrackType > >  &l1Tracks);


   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      //virtual void endRun(edm::Run&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

	float ZMAX;	// in cm
	float DeltaZ;	// in cm
	float CHI2MAX;
	float PTMINTRA ; 	// in GeV

	unsigned int nStubsmin ;		// minimum number of stubs 
	int nStubsPSmin ;	// minimum number of stubs in PS modules 
	bool SumPtSquared;
	bool MC;
  	edm::InputTag TrackingParticleInputTag;
  	edm::InputTag GenInputTag;
        edm::InputTag MCTruthStubInputTag;

  	edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  	edm::EDGetTokenT< std::vector<reco::GenParticle> > HEPMCVertexToken_;
        
        edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;

        const edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken;
	std::vector<edm::Ptr<L1TTTrackType> > tracks_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
L1TechPropPrimaryVertexProducer::L1TechPropPrimaryVertexProducer(const edm::ParameterSet& iConfig) :
  trackToken(consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag")))
{
   //register your products
   //now do what ever other initialization is needed
  
  ZMAX = (float)iConfig.getParameter<double>("ZMAX");
  DeltaZ = (float)iConfig.getParameter<double>("DeltaZ");
  CHI2MAX = (float)iConfig.getParameter<double>("CHI2MAX");
  PTMINTRA = (float)iConfig.getParameter<double>("PTMINTRA");

  nStubsmin = iConfig.getParameter<int>("nStubsmin");
  nStubsPSmin = iConfig.getParameter<int>("nStubsPSmin");
  MC=iConfig.getParameter<bool>("MonteCarloVertex");
  SumPtSquared = iConfig.getParameter<bool>("SumPtSquared");
  GenInputTag = iConfig.getParameter<edm::InputTag>("GenParticleInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  MCTruthStubInputTag = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");
  HEPMCVertexToken_=consumes< std::vector<reco::GenParticle> >(GenInputTag);
  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);  
  ttStubMCTruthToken_ = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);
  produces<L1TkPrimaryVertexCollection>();
  produces<L1TkPrimaryVertexCollection>("TrackingParticleVtx");
  produces< l1t::VertexCollection >( "l1vertices" );
}


L1TechPropPrimaryVertexProducer::~L1TechPropPrimaryVertexProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1TechPropPrimaryVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 std::unique_ptr<L1TkPrimaryVertexCollection> result(new L1TkPrimaryVertexCollection);

  
  ////////////////////////
  // GET MAGNETIC FIELD //
  edm::ESHandle<MagneticField> magneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
  const MagneticField* theMagneticField = magneticFieldHandle.product();
  double mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
  if ( mMagneticFieldStrength < 0) std::cout << "mMagneticFieldStrength < 0 " << std::endl;  // for compil when not used

  //////////////////////
  // Tracker Topology //
  edm::ESHandle<TrackerTopology> tTopoHandle_;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle_);
  const TrackerTopology* tTopo = tTopoHandle_.product();


  edm::Handle<L1TTTrackCollectionType> L1TTTrackHandle;
  iEvent.getByToken(trackToken, L1TTTrackHandle);   

  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);

   edm::Handle< std::vector< reco::GenParticle> > GenParticleHandle;
   iEvent.getByToken(HEPMCVertexToken_,GenParticleHandle);

  edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);

if(MC ){
    std::vector<reco::GenParticle>::const_iterator genpartIter ;
	float zvtx_gen = -999;
	float gensumpt=0;
        for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
           int status = genpartIter -> status() ;
	   if (status< 21 || status>29) continue;
	   //if (status!=1) continue;
	   if ( genpartIter -> numberOfMothers() == 0) continue;
	   zvtx_gen = genpartIter -> vz() ;
	   
	   break;
	} 
   for (genpartIter = GenParticleHandle->begin(); genpartIter != GenParticleHandle->end(); ++genpartIter) {
 	int status = genpartIter -> status() ;
	if(status!=1)continue;
	if(genpartIter ->charge()==0)continue;
	if(genpartIter->pt()<2 || fabs(genpartIter->eta())>2.4)continue;
	gensumpt=gensumpt+genpartIter->pt();
   }

   //sum over gen charged Particles	
   L1TkPrimaryVertex vtx1( zvtx_gen, gensumpt );
   result -> push_back( vtx1 );
   iEvent.put( std::move(result));
   return;	
}
   
 if( !L1TTTrackHandle.isValid() )
        {
          LogError("L1TechPropPrimaryVertexProducer")
            << "\nWarning: LTTkTrackCollection not found in the event. Exit"
            << std::endl;
 	    return;
        }


   float sum1 = -999;
   int nmin = nStubsmin;
   int nPSmin = nStubsPSmin ;
   float ptmin = PTMINTRA ;
   int imode = 2;	// max(Sum PT2)
   if (! SumPtSquared)  imode = 1;   // max(Sum PT)
   tracks_.clear();
   float z1 = WeightedZVertex( L1TTTrackHandle, sum1, nmin, nPSmin, ptmin, imode, tTopo,tracks_);
   L1TkPrimaryVertex vtx1( z1, sum1 );
   result -> push_back( vtx1 );
   iEvent.put( std::move(result) );
   //std::cout<<"z1 "<<z1<<" Num of tracks "<<tracks_.size()<<std::endl;
   //loop over tracks and tracking particles in the chosen Z-window and get the input tracks that pass the cut
   std::unique_ptr<l1t::VertexCollection> lProduct(new std::vector<l1t::Vertex>());
   lProduct->push_back(l1t::Vertex(z1, tracks_));
   iEvent.put(std::move(lProduct), "l1vertices"); 

   //loop over tracking particles for the True MCVtx 
   TrackingParticleVertex(iEvent,TrackingParticleHandle,MCTruthTTStubHandle,tTopo); 
}


float L1TechPropPrimaryVertexProducer::WeightedZVertex(const edm::Handle<L1TTTrackCollectionType> & L1TTTrackHandle,
                float z, int nmin, int nPSmin, float ptmin, int imode,
                const TrackerTopology* topol, std::vector< edm::Ptr< L1TTTrackType > >  &l1Tracks){
float vtxZ=-999;
 int nbins=2*ZMAX/DeltaZ;
 float xmin = -1*ZMAX ;
 float xmax = ZMAX ;
  TH1F*htmp = new TH1F("htmp",";z (cm); Tracks",nbins,xmin,xmax);
  TH1F*htmp_weight = new TH1F("htmp_weight",";z (cm); Tracks",nbins,xmin,xmax);
  l1Tracks.clear();
  std::vector< edm::Ptr< L1TTTrackType > > SelectedTracks;
  L1TTTrackCollectionType::const_iterator trackIter;
    int this_l1track = 0;
  for (trackIter = L1TTTrackHandle->begin(); trackIter != L1TTTrackHandle->end(); ++trackIter) {
	  edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(L1TTTrackHandle, this_l1track);
      this_l1track++;

    float pt = trackIter->getMomentum().perp();
    float chi2 = trackIter->getChi2();
    float ztr  = trackIter->getPOCA().z();
    if(ztr!=ztr)continue; //just if there is nan
    if (pt < ptmin) continue;
    if (fabs(ztr) > ZMAX ) continue;
    if (chi2 > CHI2MAX) continue;
   

	// get the number of stubs and the number of stubs in PS layers
    int nPS = 0.;     // number of stubs in PS modules
    int nstubs = 0;

      // get pointers to stubs associated to the L1 track
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > >  theStubs = trackIter -> getStubRefs() ;

     int tmp_trk_nstub = (int) theStubs.size();
      if ( tmp_trk_nstub < 0) {
	std::cout << " ... could not retrieve the vector of stubs in L1TechPropPrimaryVertexProducer::SumPtVertex " << std::endl;
	continue;
      }
	nstubs=theStubs.size();
      // loop over the stubs
      for (unsigned int istub=0; istub<(unsigned int)theStubs.size(); istub++) {
        //bool genuine = theStubs.at(istub)->isGenuine();
        //if (genuine) {
	   bool isPS = false;
	   DetId detId( theStubs.at(istub)->getDetId() );
	   if (detId.det() == DetId::Detector::Tracker) {
	     if (detId.subdetId() == StripSubdetector::TOB && topol->tobLayer(detId) <= 3)  isPS = true;
	     else if (detId.subdetId() == StripSubdetector::TID && topol->tidRing(detId) <= 9)  isPS = true;
	   }
	   if (isPS) nPS ++;
	   //if (isPS) cout << " this is a stub in a PS module " << endl;
           //if (isPS) nPS ++;
	//} // endif genuine
       } // end loop over stubs

        if (imode == 1 || imode == 2 ) {
            if (nPS < nPSmin) continue;
        }
	if ( nstubs < nmin) continue;
	     SelectedTracks.push_back(l1track_ptr);
	     htmp -> Fill( ztr );
     	     htmp_weight -> Fill( ztr, pt );

	}
  float zvtx_sliding = -999;
  float sigma_max = -999;
  int nb = htmp -> GetNbinsX();
  for (int i=2; i <= nb-1; i++) {
     float a0 = htmp -> GetBinContent(i-1);
     float a1 = htmp -> GetBinContent(i);
     float a2 = htmp -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp -> GetBinCenter(i-1);
        float z1 = htmp -> GetBinCenter(i);
        float z2 = htmp -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
zvtx_sliding = -999;
 sigma_max = -999;
 for (int i=2; i <= nb-1; i++) {
     float a0 = htmp_weight -> GetBinContent(i-1);
     float a1 = htmp_weight -> GetBinContent(i);
     float a2 = htmp_weight -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp_weight -> GetBinCenter(i-1);
        float z1 = htmp_weight -> GetBinCenter(i);
        float z2 = htmp_weight -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
  vtxZ=zvtx_sliding;
	
  for(unsigned int t=0; t<SelectedTracks.size(); ++t){
	float z0=SelectedTracks[t]->getPOCA().z();
	if(fabs(vtxZ-z0) >DeltaZ)continue;
	l1Tracks.push_back(SelectedTracks[t]);
}

  return vtxZ;
}
float L1TechPropPrimaryVertexProducer::MaxPtVertex(const edm::Handle<L1TTTrackCollectionType> & L1TTTrackHandle,
 		float& Sum,
		int nmin, int nPSmin, float ptmin, int imode,
		const TrackerTopology* topol){
        // return the zvtx corresponding to the max(SumPT)
        // of tracks with at least nPSmin stubs in PS modules
   
      float sumMax = 0;
      float zvtxmax = -999;
      int nIter = (int)(ZMAX* 2.)/DeltaZ ;

      for (int itest = 0; itest <= nIter; itest ++) {
	
        //float z = -100 + itest;         // z in cm
	float z = -ZMAX + itest *DeltaZ;  	// z in cm
        //z = z/10.  ;   // z in cm
   	std::vector< edm::Ptr< L1TTTrackType > > l1Tracks;
        float sum = SumPtVertex(L1TTTrackHandle, z, nmin, nPSmin, ptmin, imode, topol,l1Tracks);
        if (sumMax >0 && sum == sumMax) {
          //cout << " Note: Several vertices have the same sum " << zvtxmax << " " << z << " " << sumMax << endl;
        }
   
        if (sum > sumMax) {
           sumMax = sum;
           zvtxmax = z;
	   tracks_=l1Tracks;
        }
       }  // end loop over tested z 
   
 Sum = sumMax;
 return zvtxmax;
}  


float L1TechPropPrimaryVertexProducer::SumPtVertex(const edm::Handle<L1TTTrackCollectionType> & L1TTTrackHandle,
		float z, int nmin, int nPSmin, float ptmin, int imode,
		const TrackerTopology* topol,std::vector< edm::Ptr< L1TTTrackType > >  &l1Tracks) {

        // sumPT of tracks with >= nPSmin stubs in PS modules
        // z in cm
 float sumpt = 0;

  l1Tracks.clear();
  L1TTTrackCollectionType::const_iterator trackIter;
    int this_l1track = 0;
  for (trackIter = L1TTTrackHandle->begin(); trackIter != L1TTTrackHandle->end(); ++trackIter) {
	  edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(L1TTTrackHandle, this_l1track);
      this_l1track++;

    float pt = trackIter->getMomentum().perp();
    float chi2 = trackIter->getChi2();
    float ztr  = trackIter->getPOCA().z();

    if (pt < ptmin) continue;
    if (fabs(ztr) > ZMAX ) continue;
    if (chi2 > CHI2MAX) continue;
    if ( fabs(ztr - z) > DeltaZ) continue;   // eg DeltaZ = 1 mm


	// get the number of stubs and the number of stubs in PS layers
    int nPS = 0.;     // number of stubs in PS modules
    int nstubs = 0;

      // get pointers to stubs associated to the L1 track
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > >  theStubs = trackIter -> getStubRefs() ;

      int tmp_trk_nstub = (int) theStubs.size();
      if ( tmp_trk_nstub < 0) {
	std::cout << " ... could not retrieve the vector of stubs in L1TechPropPrimaryVertexProducer::SumPtVertex " << std::endl;
	continue;
      }
	nstubs=theStubs.size();
      // loop over the stubs
      for (unsigned int istub=0; istub<(unsigned int)theStubs.size(); istub++) {
        //bool genuine = theStubs.at(istub)->isGenuine();
        //if (genuine) {
	   bool isPS = false;
	   DetId detId( theStubs.at(istub)->getDetId() );
	   if (detId.det() == DetId::Detector::Tracker) {
	     if (detId.subdetId() == StripSubdetector::TOB && topol->tobLayer(detId) <= 3)  isPS = true;
	     else if (detId.subdetId() == StripSubdetector::TID && topol->tidRing(detId) <= 9)  isPS = true;
	   }
	   if (isPS) nPS ++;
	   //if (isPS) cout << " this is a stub in a PS module " << endl;
           //if (isPS) nPS ++;
	//} // endif genuine
       } // end loop over stubs

        if (imode == 1 || imode == 2 ) {
            if (nPS < nPSmin) continue;
        }
	if ( nstubs < nmin) continue;
	l1Tracks.push_back(l1track_ptr);
        if (imode == 2) sumpt += pt*pt;
        if (imode == 1) sumpt += pt;

  } // end loop over the tracks

 return sumpt;

}



// ------------ method called once each job just before starting event loop  ------------
void 
L1TechPropPrimaryVertexProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TechPropPrimaryVertexProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
L1TechPropPrimaryVertexProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{


}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TechPropPrimaryVertexProducer::TrackingParticleVertex(edm::Event& iEvent,edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle,
edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle,
const TrackerTopology* topol){
double vtxZ=0;
 int nbins=2*ZMAX/DeltaZ;
 float xmin = -1*ZMAX ;
 float xmax = ZMAX ;
 std::unique_ptr<L1TkPrimaryVertexCollection> result(new L1TkPrimaryVertexCollection);
  TH1F*htmp = new TH1F("htmp",";z (cm); Tracks",nbins,xmin,xmax);
  TH1F*htmp_weight = new TH1F("htmp_weight",";z (cm); Tracks",nbins,xmin,xmax);

int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;
  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
	++this_tp;
//   if(tp_ptr->pt()>200)continue;
   if(tp_ptr->pt()<PTMINTRA)continue;
   if(fabs(tp_ptr->vz())>ZMAX)continue;

    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
//    std::cout<<"Nstubs "<<theStubRefs.size()<<std::endl;

    if(theStubRefs.size()<nStubsmin)continue;
   
    // how many layers/disks have stubs?
    int hasStubInLayer[11] = {0};
    int nPS=0;
    for (unsigned int is=0; is<theStubRefs.size(); is++) {
        //bool genuine = theStubs.at(istub)->isGenuine();
        //if (genuine) {
	   bool isPS = false;
	   DetId detId( (theStubRefs.at(is)->getDetId() ));
	   if (detId.det() == DetId::Detector::Tracker) {
	     if (detId.subdetId() == StripSubdetector::TOB && topol->tobLayer(detId) <= 3)  isPS = true;
	     else if (detId.subdetId() == StripSubdetector::TID && topol->tidRing(detId) <= 9)  isPS = true;
	   }
	   int layer = -1;
	  if ( detId.subdetId()==StripSubdetector::TOB ) {

	    layer  = static_cast<int>(topol->layer(detId));
	    layer=layer-1;
	  }

	  else if ( detId.subdetId()==StripSubdetector::TID ) {
	    layer  =static_cast<int>(topol->layer(detId));
	    layer=layer+5;
	  }	  
 	   if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer]<2)hasStubInLayer[layer] = 1; //doesn't have a truth matched stub
	   else  hasStubInLayer[layer] = 2; //has 2 hits from the same TP
	   //if (isPS) nPS ++;
	   //if (isPS) cout << " this is a stub in a PS module " << endl;
           if (isPS) nPS ++;
	//} // endif genuine
        // end loop over stubs
	}
    unsigned int nStubLayerTP = 0;
    for (int isum=0; isum<11; isum++) {
        if ( hasStubInLayer[isum] >= 1) nStubLayerTP   += 1;
	}
        if(nStubLayerTP<nStubsmin)continue;
        if (nStubsPSmin > nPS) continue;
     	   htmp -> Fill( tp_ptr->vz());
           if(SumPtSquared)htmp_weight -> Fill( tp_ptr->vz(), tp_ptr->pt()* tp_ptr->pt());
	   else htmp_weight -> Fill( tp_ptr->vz(), tp_ptr->pt() );
    }
  float zvtx_sliding = -999;
  float sigma_max = -999;
  int nb = htmp -> GetNbinsX();

  for (int i=2; i <= nb-1; i++) {
     float a0 = htmp -> GetBinContent(i-1);
     float a1 = htmp -> GetBinContent(i);
     float a2 = htmp -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp -> GetBinCenter(i-1);
        float z1 = htmp -> GetBinCenter(i);
        float z2 = htmp -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
 zvtx_sliding = -999;
 sigma_max = -999;
 for (int i=2; i <= nb-1; i++) {
     float a0 = htmp_weight -> GetBinContent(i-1);
     float a1 = htmp_weight -> GetBinContent(i);
     float a2 = htmp_weight -> GetBinContent(i+1);
     float sigma = a0 + a1 + a2;
     if (sigma > sigma_max) {
        sigma_max = sigma;
        float z0 = htmp_weight -> GetBinCenter(i-1);
        float z1 = htmp_weight -> GetBinCenter(i);
        float z2 = htmp_weight -> GetBinCenter(i+1);
        zvtx_sliding =  ( a0 * z0 + a1 * z1 + a2 * z2 ) / sigma;
     }
  }
   vtxZ=zvtx_sliding;
   float gensumpt=htmp_weight->GetBinContent(htmp_weight->GetXaxis()->FindBin(vtxZ));
   L1TkPrimaryVertex vtx1( vtxZ, gensumpt );
   result -> push_back( vtx1 );
   iEvent.put( std::move(result),"TrackingParticleVtx");

}
void
L1TechPropPrimaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TechPropPrimaryVertexProducer);
