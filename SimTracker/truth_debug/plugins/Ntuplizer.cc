#include "TTree.h"
#include "TFile.h"

#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "RecoTracker/FinalTrackSelectors/interface/getBestVertex.h"

class Ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
    explicit Ntuplizer(const edm::ParameterSet&);
    ~Ntuplizer() override;
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void beginJob() override;
    // void beginRun(const edm::Run&, const edm::EventSetup&) override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    // void endRun(const edm::Run& iEvent, const edm::EventSetup&) override;
    void endJob() override;

    // tokens
    edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
    edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
    edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> assocToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

    using TrackingParticleRefKeyToIndex = std::unordered_map<reco::RecoToSimCollection::index_type, size_t>;

    // output trees
    TTree* track_tree_;
    TTree* assoc_tree_;
    TTree* tp_tree_;

    void clearVariables();

    unsigned int event_index;

    // variables for branches

    // tracks    
    std::vector<int> track_ev;
    std::vector<unsigned int> track_id;
    std::vector<float> track_pt;
    std::vector<float> track_eta;
    std::vector<float> track_phi;
    std::vector<int> track_nHits;
    std::vector<int> track_nHits_pixel;
    std::vector<int> track_nHits_strip;
    std::vector<int> track_nLost;
    std::vector<int> track_charge;
    std::vector<float> track_inner_px;
    std::vector<float> track_inner_py;
    std::vector<float> track_inner_pz;
    std::vector<float> track_inner_pt;
    std::vector<float> track_outer_px;
    std::vector<float> track_outer_py;
    std::vector<float> track_outer_pz;
    std::vector<float> track_outer_pt;
    std::vector<float> track_dxy;
    std::vector<float> track_dz;
    std::vector<float> track_dxy_err;
    std::vector<float> track_dz_err;

    // association
    std::vector<std::vector<unsigned int>> track_TP_recoToSim; // vector of TPs associated to every reco track
    std::vector<std::vector<int>> track_TP_recoToSim_qual; 

    // tracking particles
    std::vector<int> tp_ev;
    std::vector<unsigned int> tp_id;
    std::vector<float> tp_pt;
    std::vector<float> tp_eta;
    std::vector<float> tp_phi;
    std::vector<int> tp_pdgID;
    std::vector<double> tp_parentVertex_x;
    std::vector<double> tp_parentVertex_y;
    std::vector<double> tp_parentVertex_z;

};

Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) 
    : tracksToken_(consumes(iConfig.getParameter<edm::InputTag>("tracks"))),
    tpToken_(consumes<std::vector<TrackingParticle>>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
    assocToken_(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getParameter<edm::InputTag>("tkToTpAssociator"))),
    beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))) {
    };

Ntuplizer::~Ntuplizer() {
    clearVariables();
};

void Ntuplizer::clearVariables() {
    track_ev.clear();
    track_id.clear();
    track_pt.clear();
    track_eta.clear();
    track_phi.clear();
    track_nHits.clear();
    track_nHits_pixel.clear();
    track_nHits_strip.clear();
    track_nLost.clear();
    track_charge.clear();
    track_inner_px.clear();
    track_inner_py.clear();
    track_inner_pz.clear();
    track_inner_pt.clear();
    track_outer_px.clear();
    track_outer_py.clear();
    track_outer_pz.clear();
    track_outer_pt.clear();
    track_dxy.clear();
    track_dz.clear();
    track_dxy_err.clear();
    track_dz_err.clear();

    track_TP_recoToSim.clear();
    track_TP_recoToSim_qual.clear();

    tp_ev.clear();
    tp_id.clear();
    tp_pt.clear();
    tp_eta.clear();
    tp_phi.clear();
    tp_pdgID.clear();
    tp_parentVertex_x.clear();
    tp_parentVertex_y.clear();
    tp_parentVertex_z.clear();
};

void Ntuplizer::beginJob() {
    // usesResource(TFileService::kSharedResource);
    edm::Service<TFileService> fs;
//    std::cout << "making trees\n";

    // tracks
    track_tree_ = fs->make<TTree>("tracks","RECO tracks");
    
    track_tree_->Branch("event", &track_ev);
    track_tree_->Branch("id", &track_id);
    track_tree_->Branch("pt", &track_pt);
    track_tree_->Branch("eta", &track_eta);
    track_tree_->Branch("phi", &track_phi);
    track_tree_->Branch("nHits", &track_nHits);
    track_tree_->Branch("nHits_pixel", &track_nHits_pixel);
    track_tree_->Branch("nHits_strip", &track_nHits_strip);
    track_tree_->Branch("nLost", &track_nLost);
    track_tree_->Branch("charge", &track_charge);
    track_tree_->Branch("inner_px", &track_inner_px);
    track_tree_->Branch("inner_py", &track_inner_py);
    track_tree_->Branch("inner_pz", &track_inner_pz);
    track_tree_->Branch("inner_pt", &track_inner_pt);
    track_tree_->Branch("outer_px", &track_outer_px);
    track_tree_->Branch("outer_py", &track_outer_py);
    track_tree_->Branch("outer_pz", &track_outer_pz);
    track_tree_->Branch("outer_pt", &track_outer_pt);
    track_tree_->Branch("dxy", &track_dxy);
    track_tree_->Branch("dz", &track_dz);
    track_tree_->Branch("dxy_error", &track_dxy_err);
    track_tree_->Branch("dz_err", &track_dz_err);

    // associator
    assoc_tree_ = fs->make<TTree>("associations","track to TP");

    assoc_tree_->Branch("recoToSimMap", &track_TP_recoToSim);
    assoc_tree_->Branch("recoToSimQual", &track_TP_recoToSim_qual);

    // tracking particles
    tp_tree_ = fs->make<TTree>("trackingParticles", "tracking particles");

    tp_tree_->Branch("event", &tp_ev);
    tp_tree_->Branch("id", &tp_id);
    tp_tree_->Branch("pt", &tp_pt);
    tp_tree_->Branch("eta", &tp_eta);
    tp_tree_->Branch("phi", &tp_phi);
    tp_tree_->Branch("pdgID", &tp_pdgID);
    tp_tree_->Branch("parentVertex_x", &tp_parentVertex_x);
    tp_tree_->Branch("parentVertex_y", &tp_parentVertex_y);
    tp_tree_->Branch("parentVertex_z", &tp_parentVertex_z);

    // std::cout << "made trees\n";
    
};

void Ntuplizer::endJob() { };

void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& setup) {
    event_index++;
    clearVariables();
    // std::cout << "getting collections...\n";
    edm::Handle<edm::View<reco::Track>> trackCollectionH;
    iEvent.getByToken(tracksToken_, trackCollectionH);
    const auto& tracks = *(trackCollectionH.product());

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexToken_, vertices);

    const auto& beamSpot = iEvent.get(beamSpotToken_);

    edm::Handle<TrackingParticleCollection> tpCollectionH;
    iEvent.getByToken(tpToken_, tpCollectionH);
    const TrackingParticleCollection& tpCollection = *(tpCollectionH.product());

    const auto& associator = iEvent.get(assocToken_);
    const reco::RecoToSimCollection tkToTpMap = associator.associateRecoToSim(trackCollectionH, tpCollectionH);

    track_TP_recoToSim.resize(tracks.size());
    track_TP_recoToSim_qual.resize(tracks.size());

    TrackingParticleRefKeyToIndex tpKeyToIndex; // map from Ref::key() to index for TPs
    // tracking particle loop
    for (size_t i=0; i < tpCollection.size(); ++i) {
        tpKeyToIndex[TrackingParticleRef(tpCollectionH, i).key()] = i;
        const auto& tp_i = tpCollection[i];
        tp_ev.push_back(event_index);
        tp_id.push_back(i);
        tp_pt.push_back(tp_i.pt());
        tp_eta.push_back(tp_i.eta());
        tp_phi.push_back(tp_i.phi());
        tp_pdgID.push_back(tp_i.pdgId());
        tp_parentVertex_x.push_back(tp_i.vx());
        tp_parentVertex_y.push_back(tp_i.vy());
        tp_parentVertex_z.push_back(tp_i.vz());
    }

    // track loop
    // std::cout << "entering track loop\n";
    for (unsigned int i=0; i < tracks.size(); ++i) {
        // std::cout << "track " << i << std::endl;
        const auto& tk_i = tracks[i];
        const auto& hp = tk_i.hitPattern();
        // std::cout << "got track and hit pattern\n";
        track_ev.push_back(event_index);
        track_id.push_back(i);
        track_pt.push_back(tk_i.pt());
        track_eta.push_back(tk_i.eta());
        track_phi.push_back(tk_i.phi());
        track_nHits.push_back(hp.numberOfValidHits());
        track_nHits_pixel.push_back(hp.numberOfValidPixelHits());
        track_nHits_strip.push_back(hp.numberOfValidStripHits());
        track_nLost.push_back(hp.numberOfLostHits(reco::HitPattern::TRACK_HITS));
        track_charge.push_back(tk_i.charge());
        track_inner_px.push_back(tk_i.innerMomentum().x());
        track_inner_py.push_back(tk_i.innerMomentum().y());
        track_inner_pz.push_back(tk_i.innerMomentum().z());
        track_inner_pt.push_back(tk_i.innerMomentum().rho());
        track_outer_px.push_back(tk_i.outerMomentum().x());
        track_outer_py.push_back(tk_i.outerMomentum().y());
        track_outer_pz.push_back(tk_i.outerMomentum().z());
        track_outer_pt.push_back(tk_i.outerMomentum().rho());
        track_dxy.push_back(tk_i.dxy(beamSpot.position())); // would probably work without the beam spot position
        track_dz.push_back(tk_i.dz(beamSpot.position()));
        track_dxy_err.push_back(tk_i.dxyError());
        track_dz_err.push_back(tk_i.dzError());
        
        // Point bestPV = getBestVertex(trk_i, vertices); // closest PV

        // associator
        // std::cout << "associator\n";
        edm::RefToBase<reco::Track> tkRef_i(trackCollectionH, i);
        auto foundHere = tkToTpMap.find(tkRef_i);
        if (foundHere != tkToTpMap.end()) {
            auto foundVal = tkToTpMap[tkRef_i];
            if (foundVal.size() != 0) {
                // std::cout << "found sim match\n";
                for (const auto& val_i : foundVal) {
                    // track_TP_recoToSim[i].push_back(val_i.first.index());
                    track_TP_recoToSim[i].push_back(tpKeyToIndex.at(val_i.first.key()));
                    track_TP_recoToSim_qual[i].push_back(val_i.second);
                }
            }
        }

    } // track loop

    // fill trees
    track_tree_->Fill();
    assoc_tree_->Fill();
    tp_tree_->Fill();
}

void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("tkToTpAssociator", edm::InputTag("trackAssociatorByHits"));
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices")); 
  descriptions.add("Ntuplizer", desc);
}

DEFINE_FWK_MODULE(Ntuplizer);
