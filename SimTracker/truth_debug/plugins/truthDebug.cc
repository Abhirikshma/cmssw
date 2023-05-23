// -*- C++ -*-
//
// Package:    SimTracker/truthDebug
// Class:      truthDebug
//
/**\class truthDebug truthDebug.cc SimTracker/truthDebug/plugins/truthDebug.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Abhirikshma Nandi
//         Created:  Fri, 14 Apr 2023 14:55:21 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class truthDebug : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit truthDebug(const edm::ParameterSet&);
  ~truthDebug() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> assocToken_;

  const TrackerTopology *ttopo;
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
truthDebug::truthDebug(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes(iConfig.getParameter<edm::InputTag>("tracks"))),
    tpToken_(consumes<std::vector<TrackingParticle>>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
    assocToken_(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getParameter<edm::InputTag>("tkToTpAssociator")))
     {
  
  //now do what ever initialization is needed
}

truthDebug::~truthDebug() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void truthDebug::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::View<reco::Track>> trackCollectionH;
  iEvent.getByToken(tracksToken_, trackCollectionH);
  const edm::View<reco::Track> &trackCollection = *(trackCollectionH.product());

  edm::Handle<std::vector<TrackingParticle>> tpCollectionH;
  iEvent.getByToken(tpToken_, tpCollectionH);
  const std::vector<TrackingParticle> &tpCollection = *(tpCollectionH.product());

  std::cout << "# tps in collection " << tpCollection.size() << std::endl;
  std::cout << "# reco tracks in collection " << trackCollection.size() << std::endl; 

  const auto& associator = iEvent.get(assocToken_);
  const reco::RecoToSimCollection &tkToTpMap = associator.associateRecoToSim(trackCollectionH, tpCollectionH);

  // Find match to track in reco to sim collection
  for (unsigned i = 0; i < trackCollection.size(); ++i) {
    edm::RefToBase<reco::Track> tk_i(trackCollectionH, i);

    auto foundHere = tkToTpMap.find(tk_i); // check if track is in the map
    if (foundHere != tkToTpMap.end()) {

      auto foundVal = tkToTpMap[tk_i];
      if (foundVal.size() != 0) {
        // Reading the matched TPs
        std::cout << "reco track " << i << " pT: " << tk_i->pt() << " matched to " << foundVal.size() << " tracking particles\n";
        for (auto val_i : foundVal) {
          auto tp_i = val_i.first;
          double qual = val_i.second;
          std::cout << "\tTP pT: " << tp_i->pt() << " quality " << qual << " genpart status " << tp_i->status() << " pdg id: " << tp_i->pdgId() << std::endl;
          std::cout << "\t\t parent TrackingVertex position (x, y, z): " << tp_i->vx() << ", " << tp_i->vy() << ", " << tp_i->vz() << std::endl;

        }
      }
      
    }
    else std::cout << "reco track " << i << " pT: " << tk_i->pt() << " matched to 0 tracking particles\n";
  }
  
  // const auto& tkToTpMap = iEvent.get(tkToTpCollToken_);
}

// ------------ method called once each job just before starting event loop  ------------
void truthDebug::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void truthDebug::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void truthDebug::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("tkToTpAssociator", edm::InputTag("trackAssociatorByHits"));
  descriptions.add("truthDebug", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(truthDebug);
