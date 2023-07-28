#include "TTree.h"
#include "TFile.h"

#include <iostream>

#include "HepMC/GenVertex.h"

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
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"

#include "RecoTracker/FinalTrackSelectors/interface/getBestVertex.h"
#include "SimTracker/TrackHistory/interface/VertexClassifier.h"

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
    edm::EDGetTokenT<reco::VertexToTrackingVertexAssociator> vertAssocToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex>> secondaryVertexToken_;
    edm::EDGetTokenT<TrackingVertexCollection> tvToken_;
    // edm::EDGetTokenT<reco::SecondaryVertexTagInfo> secondaryVertexToken_;

    // edm::InputTag svTagInfoProducer_;

    VertexClassifier classifier_;

    bool debug_; // enables printing information for debugging

    using RefKeyToIndex = std::unordered_map<reco::RecoToSimCollection::index_type, size_t>;

    // output trees
    TTree* track_tree_;
    TTree* vert_tree_;
    TTree* assoc_tree_;
    TTree* v_assoc_tree_;
    TTree* tv_tree_;
    TTree* tp_tree_;

    void clearVariables();

    unsigned int event_index;

    // variables for branches

    // tracks    
    std::vector<float> track_pt;
    std::vector<float> track_eta;
    std::vector<float> track_phi;
    std::vector<int> track_nHits;
    std::vector<int> track_nHits_pixel;
    std::vector<int> track_nHits_strip;
    std::vector<int> track_nLost;
    std::vector<int> track_charge;
    std::vector<float> track_inner_x;
    std::vector<float> track_inner_y;
    std::vector<float> track_inner_z;
    std::vector<float> track_outer_x;
    std::vector<float> track_outer_y;
    std::vector<float> track_outer_z;
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

    // IVF vertices
    std::vector<double> vert_x;
    std::vector<double> vert_y;
    std::vector<double> vert_z;
    std::vector<double> vert_x_err;
    std::vector<double> vert_y_err;
    std::vector<double> vert_z_err;
    std::vector<int> vert_nTracks;
    std::vector<std::vector<unsigned int>> vert_tracks;
    std::vector<std::vector<float>> vert_track_weights;
    std::vector<double> vert_chi2;
    std::vector<double> vert_ndof;
    std::vector<bool> vert_valid;
    std::vector<bool> vert_fake;

    // track association
    std::vector<std::vector<unsigned int>> track_TP_recoToSim; // vector of TPs associated to every reco track
    std::vector<std::vector<double>> track_TP_recoToSim_qual; 

    // vertex (IVF) association
    std::vector<std::vector<unsigned int>> vert_TV_recoToSim;
    std::vector<std::vector<double>> vert_TV_recoToSim_qual;
    std::vector<std::vector<unsigned int>> vert_TV_simToReco;
    std::vector<std::vector<unsigned int>> vert_TV_simToReco_qual;

    // tracking vertices
    std::vector<double> tv_x;
    std::vector<double> tv_y;
    std::vector<double> tv_z;
    std::vector<bool> tv_inVol;
    std::vector<bool> tv_isSecondary;
    std::vector<bool> tv_isTertiary;
    std::vector<bool> tv_isBWeak;
    std::vector<bool> tv_isCWeak;

    // tracking particles
    std::vector<float> tp_pt;
    std::vector<float> tp_eta;
    std::vector<float> tp_phi;
    std::vector<int> tp_pdgID;
    std::vector<unsigned int> tp_parentVertex;
    std::vector<std::vector<unsigned int>> tp_daughterVertices;
    std::vector<double> tp_parentVertex_x;
    std::vector<double> tp_parentVertex_y;
    std::vector<double> tp_parentVertex_z;

};

Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) 
    : tracksToken_(consumes(iConfig.getParameter<edm::InputTag>("tracks"))),
    tpToken_(consumes<std::vector<TrackingParticle>>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
    assocToken_(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getParameter<edm::InputTag>("tkToTpAssociator"))),
    vertAssocToken_(consumes<reco::VertexToTrackingVertexAssociator>(iConfig.getParameter<edm::InputTag>("vertToTvAssociator"))),
    beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
    vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    secondaryVertexToken_(consumes(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
    tvToken_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("trackingVertices"))),
    classifier_(iConfig, consumesCollector()),
    debug_(iConfig.getParameter<bool>("enableDebug")) {
        event_index = 0;
        // svTagInfoProducer_ = iConfig.getParameter<edm::InputTag>("secondaryVertices");
        // consumes<reco::SecondaryVertexTagInfoCollection>(iConfig.getParameter<edm::InputTag>("secondaryVertices"));
    };

Ntuplizer::~Ntuplizer() {
    clearVariables();
};

void Ntuplizer::clearVariables() {
    track_pt.clear();
    track_eta.clear();
    track_phi.clear();
    track_nHits.clear();
    track_nHits_pixel.clear();
    track_nHits_strip.clear();
    track_nLost.clear();
    track_charge.clear();
    track_inner_x.clear();
    track_inner_y.clear();
    track_inner_z.clear();
    track_outer_x.clear();
    track_outer_y.clear();
    track_outer_z.clear();
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

    vert_x.clear();
    vert_y.clear();
    vert_z.clear();
    vert_x_err.clear();
    vert_y_err.clear();
    vert_z_err.clear();
    vert_nTracks.clear();
    vert_tracks.clear();
    vert_track_weights.clear();
    vert_chi2.clear();
    vert_ndof.clear();
    vert_valid.clear();
    vert_fake.clear();

    track_TP_recoToSim.clear();
    track_TP_recoToSim_qual.clear();

    vert_TV_recoToSim.clear();
    vert_TV_recoToSim_qual.clear();
    vert_TV_simToReco.clear();
    vert_TV_simToReco_qual.clear();

    tv_x.clear();
    tv_y.clear();
    tv_z.clear();
    tv_inVol.clear();
    tv_isSecondary.clear();
    tv_isTertiary.clear();
    tv_isBWeak.clear();
    tv_isCWeak.clear();

    tp_pt.clear();
    tp_eta.clear();
    tp_phi.clear();
    tp_pdgID.clear();
    tp_parentVertex.clear();
    tp_daughterVertices.clear();
    tp_parentVertex_x.clear();
    tp_parentVertex_y.clear();
    tp_parentVertex_z.clear();
};

void Ntuplizer::beginJob() {
    // usesResource(TFileService::kSharedResource);
    edm::Service<TFileService> fs;

    // tracks
    track_tree_ = fs->make<TTree>("tracks","RECO tracks");
    
    track_tree_->Branch("pt", &track_pt);
    track_tree_->Branch("eta", &track_eta);
    track_tree_->Branch("phi", &track_phi);
    track_tree_->Branch("nHits", &track_nHits);
    track_tree_->Branch("nHits_pixel", &track_nHits_pixel);
    track_tree_->Branch("nHits_strip", &track_nHits_strip);
    track_tree_->Branch("nLost", &track_nLost);
    track_tree_->Branch("charge", &track_charge);
    track_tree_->Branch("inner_x", &track_inner_x);
    track_tree_->Branch("inner_y", &track_inner_y);
    track_tree_->Branch("inner_z", &track_inner_z);
    track_tree_->Branch("outer_x", &track_outer_x);
    track_tree_->Branch("outer_y", &track_outer_y);
    track_tree_->Branch("outer_z", &track_outer_z);
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

    // IVF vertices
    vert_tree_ = fs->make<TTree>("ivf_verts", "IVF vertices");

    vert_tree_->Branch("x", &vert_x);
    vert_tree_->Branch("y", &vert_y);
    vert_tree_->Branch("z", &vert_z);
    vert_tree_->Branch("x_error", &vert_x_err);
    vert_tree_->Branch("y_error", &vert_y_err);
    vert_tree_->Branch("z_error", &vert_z_err);
    vert_tree_->Branch("nTracks", &vert_nTracks);
    vert_tree_->Branch("tracks", &vert_tracks);
    vert_tree_->Branch("track_weights", &vert_track_weights);
    vert_tree_->Branch("chi2", &vert_chi2);
    vert_tree_->Branch("ndof", &vert_ndof);
    vert_tree_->Branch("isValid", &vert_valid);
    vert_tree_->Branch("isFake", &vert_fake);

    // track associator
    assoc_tree_ = fs->make<TTree>("track_associations","track to TP");

    assoc_tree_->Branch("recoToSimMap", &track_TP_recoToSim);
    assoc_tree_->Branch("recoToSimQual", &track_TP_recoToSim_qual);

    // vertex associator
    v_assoc_tree_ = fs->make<TTree>("vertex_associations", "IVF vertices to TV");

    v_assoc_tree_->Branch("recoToSimMap", &vert_TV_recoToSim);
    v_assoc_tree_->Branch("recoToSimQual", &vert_TV_recoToSim_qual);
    v_assoc_tree_->Branch("simToRecoMap", &vert_TV_simToReco);
    v_assoc_tree_->Branch("simToRecoQual", &vert_TV_simToReco_qual);

    // tracking vertices
    tv_tree_ = fs->make<TTree>("trackingVertices", "tracking vertices");

    tv_tree_->Branch("x", &tv_x);
    tv_tree_->Branch("y", &tv_y);
    tv_tree_->Branch("z", &tv_z);
    tv_tree_->Branch("inVol", &tv_inVol);
    tv_tree_->Branch("isSecondary", &tv_isSecondary);
    tv_tree_->Branch("isTertiary", &tv_isTertiary);
    tv_tree_->Branch("isBWeak", &tv_isBWeak);
    tv_tree_->Branch("isCWeak", &tv_isCWeak);
    
    // tracking particles
    tp_tree_ = fs->make<TTree>("trackingParticles", "tracking particles");

    tp_tree_->Branch("pt", &tp_pt);
    tp_tree_->Branch("eta", &tp_eta);
    tp_tree_->Branch("phi", &tp_phi);
    tp_tree_->Branch("pdgID", &tp_pdgID);
    tp_tree_->Branch("parentVertex", &tp_parentVertex);
    tp_tree_->Branch("daughterVertices", &tp_daughterVertices);
    tp_tree_->Branch("parentVertex_x", &tp_parentVertex_x);
    tp_tree_->Branch("parentVertex_y", &tp_parentVertex_y);
    tp_tree_->Branch("parentVertex_z", &tp_parentVertex_z);    
};

void Ntuplizer::endJob() { };

void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& setup) {
    event_index++;
    clearVariables();
    classifier_.newEvent(iEvent, setup);

    edm::Handle<edm::View<reco::Track>> trackCollectionH;
    iEvent.getByToken(tracksToken_, trackCollectionH);
    const auto& tracks = *(trackCollectionH.product());

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexToken_, vertices);

    const auto& beamSpot = iEvent.get(beamSpotToken_);

    edm::Handle<TrackingParticleCollection> tpCollectionH;
    iEvent.getByToken(tpToken_, tpCollectionH);
    const TrackingParticleCollection& tpCollection = *(tpCollectionH.product());

    edm::Handle<TrackingVertexCollection> tvCollectionH;
    iEvent.getByToken(tvToken_, tvCollectionH);
    const TrackingVertexCollection& tvCollection = *(tvCollectionH.product());

    // edm::Handle<reco::SecondaryVertexTagInfoCollection> svTagInfoCollection;
    // iEvent.getByLabel(svTagInfoProducer_, svTagInfoCollection);
    edm::Handle<edm::View<reco::Vertex>> svH;
    iEvent.getByToken(secondaryVertexToken_, svH);
    const auto& secondaryVertices = *(svH.product());
    
    if (debug_) {
        std::cout << "Event " << event_index << "  N (IVF secondary vertices) " << secondaryVertices.size() << std::endl;
        std::cout << "secondary vertices in collection - positions: " << std::endl; 
        for (const auto& secV : secondaryVertices) {
            std::cout << "\t" << secV.position().x() << ", " << secV.position().y() << ", " << secV.position().z() << std::endl;
        }
    }

    // Associations: reco tracks to TrackingParticles
    const auto& associator = iEvent.get(assocToken_);
    const reco::RecoToSimCollection tkToTpMap = associator.associateRecoToSim(trackCollectionH, tpCollectionH);

    // Associations: reco IVF vertices -> TrackingVertices
    const auto& vert_associator = iEvent.get(vertAssocToken_);
    reco::VertexRecoToSimCollection vertToTvMap = vert_associator.associateRecoToSim(svH, tvCollectionH);
    reco::VertexSimToRecoCollection tvToVertMap = vert_associator.associateSimToReco(svH, tvCollectionH);
    if (debug_) {
        std::cout << "VertexRecoToSim size " << vertToTvMap.size() << std::endl;
        std::cout << "VertexSimToReco size " << tvToVertMap.size() << std::endl;
        std::cout << "Number of TrackingVertices " << tvCollection.size() << std::endl;
    }

    // making Ref::key() to index maps for secondary vertices
    // used while saving the VertexSimToReco maps 
    RefKeyToIndex svKeyToIndex;
    for (unsigned int i=0; i < secondaryVertices.size(); ++i) {
        reco::VertexBaseRef v_ref (svH, i);
        svKeyToIndex[v_ref.key()] = i;
    }

    RefKeyToIndex tvKeyToIndex; // map from Ref::key() to index for TVs
    vert_TV_simToReco.resize(tvCollection.size());
    vert_TV_simToReco_qual.resize(tvCollection.size());
    // tracking vertex loop
    for (unsigned int i=0; i < tvCollection.size(); ++i) {
        TrackingVertexRef tv_ref(tvCollectionH, i);
        tvKeyToIndex[tv_ref.key()] = i;
        const auto& tv_i = tvCollection[i];
        tv_x.push_back(tv_i.position().x());
        tv_y.push_back(tv_i.position().y());
        tv_z.push_back(tv_i.position().z());
        tv_inVol.push_back(tv_i.inVolume()); // is it inside tracker volume
        // vertex classification flags
        classifier_.evaluate(tv_ref);
        if (classifier_.is(VertexCategories::SecondaryVertex)) tv_isSecondary.push_back(true); else tv_isSecondary.push_back(false);
        if (classifier_.is(VertexCategories::TertiaryVertex)) tv_isTertiary.push_back(true); else tv_isTertiary.push_back(false);
        if (classifier_.is(VertexCategories::BWeakDecay)) tv_isBWeak.push_back(true); else tv_isBWeak.push_back(false);
        if (classifier_.is(VertexCategories::CWeakDecay)) tv_isCWeak.push_back(true); else tv_isCWeak.push_back(false);

        // vertex sim-to-reco
        auto foundHere = tvToVertMap.find(tv_ref);
        if (foundHere != tvToVertMap.end()) {
            auto foundVal = tvToVertMap[tv_ref];
            if (foundVal.size() != 0) {
                for (const auto& val_i : foundVal) {
                    unsigned int id = svKeyToIndex.at(val_i.first.key());
                    vert_TV_simToReco[i].push_back(id);
                    vert_TV_simToReco_qual[i].push_back(val_i.second);
                }
            }
        }
    }

    track_TP_recoToSim.resize(tracks.size());
    track_TP_recoToSim_qual.resize(tracks.size());
    tp_daughterVertices.resize(tpCollection.size());

    RefKeyToIndex tpKeyToIndex; // map from Ref::key() to index for TPs
    // tracking particle loop
    for (unsigned int i=0; i < tpCollection.size(); ++i) {
        tpKeyToIndex[TrackingParticleRef(tpCollectionH, i).key()] = i;
        const auto& tp_i = tpCollection[i];
        tp_pt.push_back(tp_i.pt());
        tp_eta.push_back(tp_i.eta());
        tp_phi.push_back(tp_i.phi());
        tp_pdgID.push_back(tp_i.pdgId());
        tp_parentVertex_x.push_back(tp_i.vx());
        tp_parentVertex_y.push_back(tp_i.vy());
        tp_parentVertex_z.push_back(tp_i.vz());

        tp_parentVertex.push_back(tvKeyToIndex.at(tp_i.parentVertex().key())); // at -> [] ??
        for (auto d = tp_i.decayVertices_begin(); d != tp_i.decayVertices_end(); ++d) {
            tp_daughterVertices[i].push_back(tvKeyToIndex.at(d.key()));
        }
    }

    RefKeyToIndex tkKeyToIndex; // Ref::key() to index map for tracks
    // track loop
    for (unsigned int i=0; i < tracks.size(); ++i) {
        const auto& tk_i = tracks[i];
        const auto& hp = tk_i.hitPattern();
        track_pt.push_back(tk_i.pt());
        track_eta.push_back(tk_i.eta());
        track_phi.push_back(tk_i.phi());
        track_nHits.push_back(hp.numberOfValidHits());
        track_nHits_pixel.push_back(hp.numberOfValidPixelHits());
        track_nHits_strip.push_back(hp.numberOfValidStripHits());
        track_nLost.push_back(hp.numberOfLostHits(reco::HitPattern::TRACK_HITS));
        track_charge.push_back(tk_i.charge());
        track_inner_x.push_back(tk_i.innerPosition().x());
        track_inner_y.push_back(tk_i.innerPosition().y());
        track_inner_z.push_back(tk_i.innerPosition().z());
        track_outer_x.push_back(tk_i.outerPosition().x());
        track_outer_y.push_back(tk_i.outerPosition().y());
        track_outer_z.push_back(tk_i.outerPosition().z());
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
        edm::RefToBase<reco::Track> tkRef_i(trackCollectionH, i);
        tkKeyToIndex[tkRef_i.key()] = i;
        auto foundHere = tkToTpMap.find(tkRef_i);
        if (foundHere != tkToTpMap.end()) {
            auto foundVal = tkToTpMap[tkRef_i];
            if (foundVal.size() != 0) {
                for (const auto& val_i : foundVal) {
                    track_TP_recoToSim[i].push_back(tpKeyToIndex.at(val_i.first.key())); // at -> [] ??
                    track_TP_recoToSim_qual[i].push_back(val_i.second);
                }
            }
        }

    } // track loop

    vert_TV_recoToSim.resize(secondaryVertices.size());
    vert_TV_recoToSim_qual.resize(secondaryVertices.size());
    vert_tracks.resize(secondaryVertices.size());
    vert_track_weights.resize(secondaryVertices.size());
    // secondary vertices loop
    for (unsigned int i=0; i<secondaryVertices.size(); ++i) {
        const auto& v = secondaryVertices[i];
        vert_x.push_back(v.x());
        vert_y.push_back(v.y());
        vert_z.push_back(v.z());
        vert_x_err.push_back(v.xError());
        vert_y_err.push_back(v.yError());
        vert_z_err.push_back(v.zError());
        vert_nTracks.push_back(v.nTracks()); // number of tracks which contibuted to vertex with weight > 0.5
        vert_ndof.push_back(v.ndof());
        vert_chi2.push_back(v.chi2());
        vert_valid.push_back(v.isValid());
        vert_fake.push_back(v.isFake());

        reco::VertexBaseRef v_ref (svH, i);

        // tracks in vertex
        if (debug_) std::cout << "tracks in vertex " << i << ":" << std::endl; 
        for (auto tk_i = v.tracks_begin(); tk_i != v.tracks_end(); ++tk_i) {
            vert_tracks[i].push_back(tkKeyToIndex.at(tk_i->key()));
            vert_track_weights[i].push_back(v.trackWeight(*tk_i));
            if (debug_)
            std::cout << "tk " <<tkKeyToIndex.at(tk_i->key()) << " wt " << v.trackWeight(*tk_i) << "  ";
        }

        if (debug_) {
            std::cout << "\nreco vert posn. " << v.x() << " " << v.y() << " " << v.z() << "  matches: "<< std::endl;
            if (v.isValid()) std::cout << "valid!\n";
            if (v.isFake()) std::cout << "fake!\n"; 
        }
        
        // reco vert -> find in map -> TVs (first), qualities (second)
        auto foundHere = vertToTvMap.find(v_ref);
        if (foundHere != vertToTvMap.end()) {
            auto foundVal = vertToTvMap[v_ref];
            if (foundVal.size() != 0) {
                for (const auto& val_i : foundVal) {
                    unsigned int id = tvKeyToIndex.at(val_i.first.key()); // ask in map for index of TVs
                    if (debug_) std::cout << "index " << id << std::endl;
                    vert_TV_recoToSim[i].push_back(id); 
                    vert_TV_recoToSim_qual[i].push_back(val_i.second);
                    
                    if (debug_) {
                        TrackingVertexRef tv_ref (tvCollectionH, id);
                        const auto& tv = tvCollection[id];
                        std::cout << "\tTV posn. " << tv.position().x() << " " << tv.position().y() << " " << tv.position().z() << "  qual " << val_i.second << std::endl;
                        classifier_.evaluate(tv_ref);
                        if (classifier_.is(VertexCategories::SecondaryVertex)) std::cout << "\tSecondary vertex" << std::endl;
                        if (classifier_.is(VertexCategories::BWeakDecay)) std::cout << "\tB weak decay" << std::endl;
                        if (classifier_.is(VertexCategories::CWeakDecay)) std::cout << "\tC weak decay" << std::endl;
                    }
                }
            }
        }
    }


    // fill trees
    track_tree_->Fill();
    vert_tree_->Fill();
    assoc_tree_->Fill();
    v_assoc_tree_->Fill();
    tv_tree_->Fill();
    tp_tree_->Fill();
}

void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("tkToTpAssociator", edm::InputTag("trackAssociatorByHits"));
  desc.add<edm::InputTag>("vertToTvAssociator", edm::InputTag("VertexAssociator"));
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("secondaryVertices", edm::InputTag("inclusiveSecondaryVertices"));
  desc.add<edm::InputTag>("trackingVertices", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<bool>("enableDebug", false); 
  descriptions.add("Ntuplizer", desc);
}

DEFINE_FWK_MODULE(Ntuplizer);
