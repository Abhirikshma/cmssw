//
// Original Author:  Marco Rovere
//         Created:  Fri May  1 07:21:02 CEST 2020
//
//
//
// system include files
#include <memory>
#include <iostream>
#include <numeric>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
//
// class declaration
//

class TiclDebugger : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit TiclDebugger(const edm::ParameterSet&);
  ~TiclDebugger() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  const edm::InputTag trackstersMerge_;
  const edm::InputTag trackstersClue3d_;
  const edm::InputTag ticlGraph_;
  const edm::InputTag ticlCandidates_;
  const edm::InputTag tracks_;
  const edm::InputTag caloParticles_;
  const edm::InputTag layerClusters_;
  const edm::InputTag associatorRecoSim_;
  const edm::InputTag associatorSimReco_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_token_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersMergeToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_clue3d_token_;
  edm::EDGetTokenT<TICLGraph> ticlGraph_token_;
  edm::EDGetTokenT<std::vector<TICLCandidate>> ticlCandidateToken_; 
  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClustersToken_;
  edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimToken_;

  TTree *tree_;

  unsigned nTracksters;
  std::vector<unsigned> indices;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> eta;
  std::vector<double> phi;
  std::vector<double> raw_E;
  std::vector<double> raw_em_E;
  std::vector<double> time;
  std::vector<std::vector<unsigned>> linked_inners;
  std::vector<std::vector<unsigned>> linked_outers;
  std::vector<std::vector<unsigned>> tracksters_in_candidate;

  std::vector<std::vector<unsigned>> tracksters_recoToSim;
  std::vector<std::vector<unsigned>> tracksters_recoToSim_score;

};

TiclDebugger::TiclDebugger(const edm::ParameterSet& iConfig)
    : trackstersMerge_(iConfig.getParameter<edm::InputTag>("trackstersMerge")),
      trackstersClue3d_(iConfig.getParameter<edm::InputTag>("trackstersClue3d")),
      ticlGraph_(iConfig.getParameter<edm::InputTag>("ticlgraph")),
      ticlCandidates_(iConfig.getParameter<edm::InputTag>("ticlcandidates")),
      tracks_(iConfig.getParameter<edm::InputTag>("tracks")),
      caloParticles_(iConfig.getParameter<edm::InputTag>("caloParticles")),
      layerClusters_(iConfig.getParameter<edm::InputTag>("layerClusters")),
      associatorRecoSim_(iConfig.getParameter<edm::InputTag>("recoToSimAssociator")),
      caloGeometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()) {
  edm::ConsumesCollector&& iC = consumesCollector();
  trackstersMergeToken_ = iC.consumes<std::vector<ticl::Trackster>>(trackstersMerge_);
  tracksters_clue3d_token_ = iC.consumes<std::vector<ticl::Trackster>>(trackstersClue3d_);
  ticlGraph_token_ = iC.consumes<TICLGraph>(ticlGraph_);
  ticlCandidateToken_ = iC.consumes<std::vector<TICLCandidate>>(ticlCandidates_);
  tracksToken_ = iC.consumes<std::vector<reco::Track>>(tracks_);
  caloParticlesToken_ = iC.consumes<std::vector<CaloParticle>>(caloParticles_);
  layerClustersToken_ = iC.consumes<std::vector<reco::CaloCluster>>(layerClusters_);
  tsRecoToSimToken_ = iC.consumes<hgcal::RecoToSimCollectionSimTracksters>(associatorRecoSim_);
  
}

TiclDebugger::~TiclDebugger() {}

void TiclDebugger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "begin analyze" <<std::endl;
  static const char* particle_kind[] = {"gam", "e", "mu", "pi0", "h", "h0", "?", "!"};
  using namespace edm;
  using std::begin;
  using std::end;
  using std::iota;
  using std::sort;

  edm::Handle<std::vector<ticl::Trackster>> trackstersMergeH;

  iEvent.getByToken(trackstersMergeToken_, trackstersMergeH);
  auto const& tracksters = *trackstersMergeH.product();
  std::vector<int> sorted_tracksters_idx(tracksters.size());
  iota(begin(sorted_tracksters_idx), end(sorted_tracksters_idx), 0);
  sort(begin(sorted_tracksters_idx), end(sorted_tracksters_idx), [&tracksters](int i, int j) {
    return tracksters[i].raw_energy() > tracksters[j].raw_energy();
  });

  edm::Handle<std::vector<ticl::Trackster>> trackstersclue3dH;
  iEvent.getByToken(tracksters_clue3d_token_, trackstersclue3dH);
  auto const& trackstersclue3d = *trackstersclue3dH.product();

  edm::Handle<TICLGraph> ticlGraphH;
  iEvent.getByToken(ticlGraph_token_, ticlGraphH);
  const auto& ticlGraph = *ticlGraphH.product();

  edm::Handle<std::vector<TICLCandidate>> ticlCandidateH;
  iEvent.getByToken(ticlCandidateToken_, ticlCandidateH);
  const auto& ticlCandidates = *ticlCandidateH.product();

  edm::Handle<std::vector<reco::CaloCluster>> layerClustersH;
  iEvent.getByToken(layerClustersToken_, layerClustersH);
  auto const& layerClusters = *layerClustersH.product();

  edm::Handle<std::vector<reco::Track>> tracksH;
  iEvent.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH.product();

  edm::Handle<std::vector<CaloParticle>> caloParticlesH;
  iEvent.getByToken(caloParticlesToken_, caloParticlesH);
  auto const& caloParticles = *caloParticlesH.product();
  std::vector<std::pair<int, float>> bestCPMatches;

  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimH;
  iEvent.getByToken(tsRecoToSimToken_, tsRecoToSimH);
  auto const& tsRecoSimMap = *tsRecoToSimH;

  nTracksters = 0;
  indices.clear();
  x.clear();
  y.clear();
  z.clear();
  eta.clear();
  phi.clear();
  raw_E.clear();
  raw_em_E.clear();
  time.clear();
  linked_inners.clear();
  linked_outers.clear();
  tracksters_in_candidate.clear();

  nTracksters = trackstersclue3d.size();
  linked_inners.resize(nTracksters);
  linked_outers.resize(nTracksters);
  tracksters_recoToSim.resize(nTracksters);
  tracksters_recoToSim_score.resize(nTracksters, -1.);
  for (unsigned i = 0; i < nTracksters; ++i) {
    std::cout << i << std::endl;
    const auto &t = trackstersclue3d[i];
    indices.push_back(i);
    const auto &barycenter = t.barycenter();
    x.push_back(barycenter.x());
    y.push_back(barycenter.y());
    z.push_back(barycenter.z());
    eta.push_back(barycenter.eta());
    phi.push_back(barycenter.phi());
    raw_E.push_back(t.raw_energy());
    raw_em_E.push_back(t.raw_em_energy());
    time.push_back(t.time());

    // TICLGraph
    const auto& node = ticlGraph.getNode((int) i);
    auto this_inners = node.getInner();
    auto this_outers = node.getOuter();
    linked_inners[i].insert(linked_inners[i].end(), this_inners.begin(), this_inners.end());
    linked_outers[i].insert(linked_outers[i].end(), this_outers.begin(), this_outers.end());

    // Associations
    const edm::Ref<ticl::Trackster> tsRef(trackstersclue3dH, i);
    const auto& stsAssociated = tsRecoSimMap.find(tsRef);

    if (stsAssociated == tsRecoSimMap.end()) continue; // no matches
    for (auto &sts : stsAssociated) {
      std::cout << sts.first << " " << sts.second << "  ";
      tracksters_recoToSim[i].push_back(sts.first);
      tracksters_recoToSim_score[i].push_back(sts.second);  
    }
  }

  unsigned nCandidates = ticlCandidates.size();
  tracksters_in_candidate.resize(nCandidates);
  for (unsigned i = 0; i < nCandidates; ++i) {
    const auto &c = ticlCandidates[i];
    auto trackster_ptrs = c.tracksters();
    for (auto ts_ptr : trackster_ptrs) {
      unsigned ts_idx = ts_ptr.get() - (edm::Ptr<ticl::Trackster>(trackstersclue3dH, 0)).get();
      tracksters_in_candidate[i].push_back(ts_idx);
    }
  }

  tree_->Fill();

  auto bestCaloParticleMatches = [&](const ticl::Trackster& t) -> void {
    bestCPMatches.clear();
    auto idx = 0;
    auto separation = 0.;
    for (auto const& cp : caloParticles) {
      separation = reco::deltaR2(t.barycenter(), cp.momentum());
      if (separation < 0.05) {
        bestCPMatches.push_back(std::make_pair(idx, separation));
      }
      ++idx;
    }
  };

  std::stringstream prob_id_str;

  for (auto const& t : sorted_tracksters_idx) {
    auto const& trackster = tracksters[t];
    auto const& probs = trackster.id_probabilities();
    // Sort probs in descending order
    std::vector<int> sorted_probs_idx(probs.size());
    iota(begin(sorted_probs_idx), end(sorted_probs_idx), 0);
    sort(begin(sorted_probs_idx), end(sorted_probs_idx), [&probs](int i, int j) { return probs[i] > probs[j]; });
    // Sort edges in ascending order
    std::vector<int> sorted_edges_idx(trackster.edges().size());
    iota(begin(sorted_edges_idx), end(sorted_edges_idx), 0);
    sort(begin(sorted_edges_idx), end(sorted_edges_idx), [&](int i, int j) {
      int layers = rhtools_.lastLayer();
      auto const& ed_i = trackster.edges()[i];
      auto const& ed_j = trackster.edges()[j];
      auto const& cl_i_in = layerClusters[ed_i[0]].hitsAndFractions()[0].first;
      auto const& cl_i_out = layerClusters[ed_i[1]].hitsAndFractions()[0].first;
      auto const& cl_j_in = layerClusters[ed_j[0]].hitsAndFractions()[0].first;
      auto const& cl_j_out = layerClusters[ed_j[1]].hitsAndFractions()[0].first;
      auto const layer_i_in = rhtools_.getLayerWithOffset(cl_i_in) + layers * ((rhtools_.zside(cl_i_in) + 1) >> 1) - 1;
      auto const layer_i_out =
          rhtools_.getLayerWithOffset(cl_i_out) + layers * ((rhtools_.zside(cl_i_out) + 1) >> 1) - 1;
      auto const layer_j_in = rhtools_.getLayerWithOffset(cl_j_in) + layers * ((rhtools_.zside(cl_j_in) + 1) >> 1) - 1;
      auto const layer_j_out =
          rhtools_.getLayerWithOffset(cl_j_out) + layers * ((rhtools_.zside(cl_j_out) + 1) >> 1) - 1;
      if (layer_i_in != layer_j_in)
        return layer_i_in < layer_j_in;
      else
        return layer_i_out < layer_j_out;
    });

    for (auto p_idx : sorted_probs_idx) {
      prob_id_str << "(" << particle_kind[p_idx] << "):" << std::fixed << std::setprecision(4) << probs[p_idx] << " ";
    }
    LogVerbatim("TICLDebugger") << "\nTrksIdx: " << t << "\n bary: " << trackster.barycenter()
                                << " baryEta: " << trackster.barycenter().eta()
                                << " baryPhi: " << trackster.barycenter().phi()
                                << "\n raw_energy: " << trackster.raw_energy()
                                << " raw_em_energy: " << trackster.raw_em_energy()
                                << "\n raw_pt: " << trackster.raw_pt() << " raw_em_pt: " << trackster.raw_em_pt()
                                << "\n seedIdx: " << trackster.seedIndex() << "\n Probs: " << prob_id_str.str();
    prob_id_str.str("");
    prob_id_str.clear();
    LogVerbatim("TICLDebugger") << "\n time: " << trackster.time() << "+/-" << trackster.timeError() << std::endl
                                << " vertices: " << trackster.vertices().size() << " average usage: "
                                << std::accumulate(std::begin(trackster.vertex_multiplicity()),
                                                   std::end(trackster.vertex_multiplicity()),
                                                   0.) /
                                       trackster.vertex_multiplicity().size()
                                << std::endl;
    LogVerbatim("TICLDebugger") << " link connections: " << trackster.edges().size() << std::endl;
    auto dumpLayerCluster = [&layerClusters](hgcal::RecHitTools const& rhtools, int cluster_idx) {
      auto const& cluster = layerClusters[cluster_idx];
      const auto firstHitDetId = cluster.hitsAndFractions()[0].first;
      int layers = rhtools.lastLayer();
      int lcLayerId =
          rhtools.getLayerWithOffset(firstHitDetId) + layers * ((rhtools.zside(firstHitDetId) + 1) >> 1) - 1;

      LogVerbatim("TICLDebugger") << "Idx: " << cluster_idx << "(" << lcLayerId << ", "
                                  << cluster.hitsAndFractions().size() << ", " << cluster.position() << ") ";
    };
    for (auto link : sorted_edges_idx) {
      LogVerbatim("TICLDebugger") << "(" << trackster.edges()[link][0] << ", " << trackster.edges()[link][1] << ")  ";
      dumpLayerCluster(rhtools_, trackster.edges()[link][0]);
      dumpLayerCluster(rhtools_, trackster.edges()[link][1]);
      LogVerbatim("TICLDebugger") << std::endl;
    }
    if (trackster.seedID().id() != 0) {
      auto const& track = tracks[trackster.seedIndex()];
      LogVerbatim("TICLDebugger") << " Seeding Track:" << std::endl;
      LogVerbatim("TICLDebugger") << "   p: " << track.p() << " pt: " << track.pt() << " charge: " << track.charge()
                                  << " eta: " << track.eta() << " outerEta: " << track.outerEta()
                                  << " phi: " << track.phi() << " outerPhi: " << track.outerPhi() << std::endl;
    }
    bestCaloParticleMatches(trackster);
    if (!bestCPMatches.empty()) {
      LogVerbatim("TICLDebugger") << " Best CaloParticles Matches:" << std::endl;
      ;
      for (auto const& i : bestCPMatches) {
        auto const& cp = caloParticles[i.first];
        LogVerbatim("TICLDebugger") << "   " << i.first << "(" << i.second << "):" << cp.pdgId()
                                    << " simCl size:" << cp.simClusters().size() << " energy:" << cp.energy()
                                    << " pt:" << cp.pt() << " momentum:" << cp.momentum() << std::endl;
      }
      LogVerbatim("TICLDebugger") << std::endl;
    }
  }
}

void TiclDebugger::beginRun(edm::Run const&, edm::EventSetup const& es) {
  const CaloGeometry& geom = es.getData(caloGeometry_token_);
  rhtools_.setGeometry(geom);
}

void TiclDebugger::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("EventTree", "Event Data");
  std::cout << "made tree" <<std::endl;

  tree_->Branch("Trackster_N", &nTracksters);
  tree_->Branch("indices", &indices);
  tree_->Branch("x", &x);
  tree_->Branch("y", &y);
  tree_->Branch("z", &z);
  tree_->Branch("eta", &eta);
  tree_->Branch("phi", &phi);
  tree_->Branch("raw_energy", &raw_E);
  tree_->Branch("raw_em_energy", &raw_em_E);
  tree_->Branch("time", &time);
  tree_->Branch("linked_inners", &linked_inners);
  tree_->Branch("linked_outers", &linked_outers);
  tree_->Branch("tracksters_in_candidate", &tracksters_in_candidate);
}

void TiclDebugger::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TiclDebugger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersMerge", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("trackstersClue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("ticlgraph", edm::InputTag("ticlGraph"));
  desc.add<edm::InputTag>("ticlcandidates", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("recoToSimAssociator", edm::InputTag("tracksterSimTracksterAssociationPR"));
  descriptions.add("ticlDebugger", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(TiclDebugger);
