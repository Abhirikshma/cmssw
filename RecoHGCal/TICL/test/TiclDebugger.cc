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
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/TICLSeedingRegion.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

//#include "RecoHGCal/TICL/interface/GlobalCache.h"
#include "RecoHGCal/TICL/plugins/SeedingRegionAlgoBase.h"
#include "RecoHGCal/TICL/plugins/SeedingRegionByTracks.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH2D.h"

//
// class declaration
//

class TiclDebugger : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit TiclDebugger(const edm::ParameterSet&);
  ~TiclDebugger() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  typedef math::XYZVector Vector;
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  void linkTracksters(const edm::Event &, const edm::EventSetup &, const std::vector<ticl::Trackster> &, edm::Handle<std::vector<CaloParticle>>);
  void buildFirstLayers();
  bool isCP(const ticl::Trackster& t, edm::Handle<std::vector<CaloParticle>>);

  const edm::InputTag trackstersMerge_;
  const edm::InputTag tracks_;
  const edm::InputTag track_col_;
  const edm::InputTag caloParticles_;
  const edm::InputTag layerClusters_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_token_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersMergeToken_;
  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> track_col_token_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClustersToken_;

  const std::string detector_;
  const std::string propName_;

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;

  std::once_flag initializeGeometry_;

  const HGCalDDDConstants* hgcons_;
  const StringCutObjectSelector<reco::Track> cutTk_;

  std::unique_ptr<GeomDet> firstDisk_[2];

  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
  MagneticField const* bFieldProd_;

  TTree* tree_;
  TH1D* h_distAll;
  TH1D* h_distCaloP;
  TH1D* h_NTracks;
  TH1D* h_NCaloParticles;
  TH1D* h_binDiff_track_eta;
  TH1D* h_binDiff_track_phi;
  TH1D* h_binDiff_track_global;
  TH1D* h_binDiff_trackster_eta;
  TH1D* h_binDiff_trackster_phi;
  TH1D* h_binDiff_trackster_global;
  TH1D* h_tksters_linked_eta;
  TH1D* h_tksters_tot_eta;
  TH1D* h_eff_eta;
  TH1D* h_eff_phi;
  TH1D* h_tksters_tot_phi;
  TH1D* h_eff_pT;
  TH1D* h_tksters_tot_pT;
  TH1D* h_eff_E;
  TH1D* h_tksters_tot_E;

  int totTracksters;
  int Ntsos_notvalid;
  int noTracks;
};

TiclDebugger::TiclDebugger(const edm::ParameterSet& iConfig)
    : trackstersMerge_(iConfig.getParameter<edm::InputTag>("trackstersMerge")),
      tracks_(iConfig.getParameter<edm::InputTag>("tracks")),
      track_col_(iConfig.getParameter<edm::InputTag>("tracks")),
      caloParticles_(iConfig.getParameter<edm::InputTag>("caloParticles")),
      layerClusters_(iConfig.getParameter<edm::InputTag>("layerClusters")),
      caloGeometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()), 
      detector_(iConfig.getParameter<std::string>("detector")),
      propName_(iConfig.getParameter<std::string>("propagator")),
      bfield_token_(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagator_token_(esConsumes<Propagator, TrackingComponentsRecord>(
        edm::ESInputTag("",propName_))),
      cutTk_(iConfig.getParameter<std::string>("cutTk")) {
  edm::ConsumesCollector&& iC = consumesCollector();
  trackstersMergeToken_ = iC.consumes<std::vector<ticl::Trackster>>(trackstersMerge_);
  tracksToken_ = iC.consumes<std::vector<reco::Track>>(tracks_);
  track_col_token_ = iC.consumes<reco::TrackCollection>(track_col_);
  caloParticlesToken_ = iC.consumes<std::vector<CaloParticle>>(caloParticles_);
  layerClustersToken_ = iC.consumes<std::vector<reco::CaloCluster>>(layerClusters_);


  std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
  hdc_token_ = iC.esConsumes<HGCalDDDConstants, IdealGeometryRecord>(
    edm::ESInputTag("",detectorName_));

  totTracksters = 0;
  Ntsos_notvalid = 0;
  noTracks = 0;
}

TiclDebugger::~TiclDebugger() {}

bool TiclDebugger::isCP(const ticl::Trackster& t, edm::Handle<std::vector<CaloParticle>> caloParticlesH) {
  //auto const& caloParticles = *caloParticlesH.product();
  edm::ProductID cp_id = caloParticlesH.id();
  //std::cout <<"in isCP, cp_id = " << cp_id;
  auto tracksterID = t.seedID();
  //std::cout <<" trackster seed ID = " << tracksterID << std::endl;
  auto index = t.seedIndex();
  bool result = false;
  if (tracksterID == cp_id) {
    result = true;
    auto cp = (*caloParticlesH)[index];
  }
  return result;
}

void TiclDebugger::buildFirstLayers() {
  float zVal = hgcons_->waferZ(1, true);
  std::pair<double, double> rMinMax = hgcons_->rangeR(zVal, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());
  }
}

void TiclDebugger::linkTracksters(const edm::Event &evt, 
                                  const edm::EventSetup &es, 
                                  const std::vector<ticl::Trackster> &tracksters, 
                                  edm::Handle<std::vector<CaloParticle>> caloParticlesH) {
  edm::ESHandle<HGCalDDDConstants> hdc  = es.getHandle(hdc_token_);
  //edm::ESHandle<MagneticField> bfield = es.getHandle(bfield_token_);
  edm::ESHandle<Propagator> propagator = es.getHandle(propagator_token_);

  hgcons_ = hdc.product();

  buildFirstLayers();

  edm::Handle<reco::TrackCollection> tracks_h;
  evt.getByToken(track_col_token_, tracks_h);
  // edm::ProductID trkId = tracks_h.id();
  //auto bFieldProd = bfield.product();
  const Propagator &prop = (*propagator);

  std::vector<std::pair<std::pair<Vector,Vector>, int>> trackPcol; // propagated track points, errors and index of track in collection
  std::vector<std::pair<Vector, int>> tracksterPcol; // same but for Tracksters

  ticl::TICLLayerTile tracksterProp_tile;
  ticl::TICLLayerTile trackProp_tile;

  // distance in cm
  // Propagate tracks to the HGCal first surface

  unsigned nTracks = tracks_h->size();
  std::cout << std::endl << std::endl << "---NEW EVENT---" << std::endl;
  //std::cout << "nTracks = " << nTracks << std::endl;
  h_NTracks->Fill(nTracks);

  for (unsigned i = 0; i < nTracks; ++i) {
    const reco::Track &tk = (*tracks_h)[i];
    if (!cutTk_((tk))) {
      ++noTracks;
      continue;
    }

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd_);
  int iSide = int(tk.eta() > 0);
  TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
  std::cout << "propagate!" << std::endl;
  if (tsos.isValid()) {
    // do something with the propagated trajectories
    std::cout << "x: " << tsos.globalPosition().x() <<" y: "<< tsos.globalPosition().y() << " z: " << tsos.globalPosition().z() <<std::endl;
    std::cout << "error : " << tsos.localError().positionError() << std::endl;

    Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
    double x_err = pow(tsos.localError().positionError().xx(), 0.5);
    double y_err = pow(tsos.localError().positionError().yy(), 0.5);
    Vector track_err(x_err, y_err, 0.0);

    h_binDiff_track_eta->Fill(abs(trackProp_tile.etaBin((trackP-track_err).Eta()) - trackProp_tile.etaBin((trackP+track_err).Eta())));
    h_binDiff_track_phi->Fill(abs(trackProp_tile.phiBin((trackP-track_err).Phi()) - trackProp_tile.phiBin((trackP+track_err).Phi())));
    h_binDiff_track_global->Fill(abs(trackProp_tile.globalBin((trackP-track_err).Eta(), (trackP-track_err).Phi()) - trackProp_tile.globalBin((trackP+track_err).Eta(), (trackP+track_err).Phi())));
    
    trackPcol.push_back(std::make_pair(std::make_pair(trackP,track_err), i));
    trackProp_tile.fill(trackP.Eta(), trackP.Phi(), i);
  }
  else {
    ++Ntsos_notvalid;
  }
  }

  // Propagate tracksters (barycenters) to the same plane

  const auto &trackstersLink = tracksters;
  std::cout << "# tracksters = " << trackstersLink.size() << std::endl;

  for (unsigned i = 0; i < trackstersLink.size(); ++i) {
  //for (auto &t : tracksters) {
    std::cout << "tracksters loop" << std::endl;
    const auto &t = trackstersLink[i];
    Vector baryc = t.barycenter();
    Vector directnv = t.eigenvectors(0);
    Vector trackster_err(pow(t.sigmas()[0],0.5), pow(t.sigmas()[1], 0.5), 0.0);
    h_tksters_tot_eta->Fill(t.barycenter().Eta());
    h_tksters_tot_phi->Fill(t.barycenter().Phi());
    h_tksters_tot_pT->Fill(t.raw_pt());
    h_tksters_tot_E->Fill(t.raw_energy());
    
    
    float zVal = hgcons_->waferZ(1, true);
    zVal *= (baryc.Z() >= 0) ? 1 : -1;

    double par = (zVal - baryc.Z())/directnv.Z();

    double xOnSurface = par*directnv.X() + baryc.X();
    double yOnSurface = par*directnv.Y() + baryc.Y();

    std::cout << "x = " << xOnSurface << std::endl;
    std::cout << "y = " << yOnSurface << std::endl;
    std::cout << "z = " << zVal << std::endl;
    std::cout << "error :  " << t.sigmas()[0] << "  " << t.sigmas()[1] << "  " << t.sigmas()[2] << std::endl;
    //std::cout << "PCA err :  " << t.sigmasPCA()[0] << "  " << t.sigmasPCA()[1] << "  " << t.sigmasPCA()[2] << std::endl;

    Vector tracksterP(xOnSurface, yOnSurface, zVal);
    tracksterPcol.push_back(std::make_pair(tracksterP, i));
    tracksterProp_tile.fill(tracksterP.Eta(), tracksterP.Phi(), i);

    h_binDiff_trackster_eta->Fill(abs(tracksterProp_tile.etaBin((tracksterP-trackster_err).Eta()) - tracksterProp_tile.etaBin((tracksterP+trackster_err).Eta())));
    h_binDiff_trackster_phi->Fill(abs(tracksterProp_tile.phiBin((tracksterP-trackster_err).Phi()) - tracksterProp_tile.phiBin((tracksterP+trackster_err).Phi())));
    h_binDiff_trackster_global->Fill(abs(tracksterProp_tile.globalBin((tracksterP-trackster_err).Eta(), (tracksterP-trackster_err).Phi()) - 
                                                          tracksterProp_tile.globalBin((tracksterP+trackster_err).Eta(), (tracksterP+trackster_err).Phi())));
  }

  // finding closest tracks to tracksters

  //std::vector<bool> tracks_mask (trackPcol.size(), 0);
  totTracksters+=tracksterPcol.size();
  std::vector<int> closestTracks (tracksterPcol.size(), -1); // index of closest track for every trackster, arranged as in tracksterPcol
  int nth = 0;
  auto const& caloParticles = *caloParticlesH.product();
  std::cout << "n caloparticles = " << caloParticles.size() << std::endl;
  h_NCaloParticles->Fill(caloParticles.size());
  /*std::cout << "pdg id: ";
  for (const auto &cp : caloParticles) {
    std::cout << cp.pdgId() <<"  ";
  }*/
  std::cout << std::endl;

  for (auto const i : tracksterPcol) {
    double minSep = 99999.;
    for (auto const j : trackPcol) {
      double sep = pow((i.first - j.first.first).Mag2(),0.5);
      //std::cout << "sep= " << sep << std::endl;
      if (sep < minSep) {
        minSep = sep;
        closestTracks[nth] = j.second;
      }
    }
    if (minSep > 0 && minSep < 1000) {
      h_distAll->Fill(minSep);
    }
    else std::cout << "sep=" << minSep;

    /*const auto &t = trackstersLink[i.second];
    if (isCP(t, caloParticlesH)) {
      //h_distCaloP->Fill(minSep);
      std::cout << i.second << "from Caloparticle" << std::endl;
    }
    if (caloParticles[0].g4Tracks()[0].crossedBoundary() || caloParticles[1].g4Tracks()[0].crossedBoundary()) {
      h_distCaloP->Fill(minSep);
      std::cout << i.second << "from Caloparticle" << std::endl;
    }
    //2 tracks, both propagated 
    if (nTracks == 2 && trackPcol.size() == 2) {
      h_distCaloP->Fill(minSep);
    }*/
    //2 tracks
    if (nTracks == 2) {
      h_distCaloP->Fill(minSep);
    }
    ++nth;
  } // trackster propagated point collection loop end

  std::cout << "closest tracks to tracksters: "; 
  for (int i : closestTracks) {
    std::cout << i <<" ";
  }
  std::cout << std::endl;


  // Search box over tiles + preliminary linking

  //std::vector<int> tracksters_near[trackPcol.second.max()+1]
  std::vector<unsigned> tracksters_near[trackPcol.size()]; // vector element for every track contains the indices of tracksters near to that track

  nth = 0;
  for (auto i : trackPcol) {
    auto trackP = i.first.first;
    double tk_eta = trackP.Eta();
    double tk_phi = trackP.Phi();
    auto track_err = i.first.second;
    //double delta = 0.05;

    int etaBin = tracksterProp_tile.etaBin(tk_eta);
    int phiBin = tracksterProp_tile.phiBin(tk_phi);
    //std::cout << etaBin <<"  "<< phiBin << std::endl;

    int delta_eta = abs(trackProp_tile.etaBin((trackP-track_err).Eta()) - trackProp_tile.etaBin((trackP+track_err).Eta()));
    int delta_phi = abs(trackProp_tile.phiBin((trackP-track_err).Phi()) - trackProp_tile.phiBin((trackP+track_err).Phi()));

    //std::array<int, 4> search_box = tracksterProp_tile.searchBoxEtaPhi(tk_eta - delta, tk_eta + delta, tk_phi - delta, tk_phi + delta);
    std::array<int, 4> search_box = {{etaBin - delta_eta, etaBin + delta_eta, phiBin - delta_phi, phiBin + delta_phi}};

    for (int eta_i = search_box[0]; eta_i < search_box[1] + 1; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i < search_box[3] + 1; ++phi_i) {
        auto tracksters_in_tile = tracksterProp_tile[tracksterProp_tile.globalBin(eta_i,phi_i)]; 
        //tracksters_near[nth].push_back(trackster_id);
        tracksters_near[nth].insert(std::end(tracksters_near[nth]), std::begin(tracksters_in_tile), std::end(tracksters_in_tile));
      }
    } // search in box ends
    ++nth;
  } // track loop ends

  std::vector<unsigned> linked_unique;
  std::cout << std::endl << "box search results (track : tracksters): " << std::endl;
  int blah  = trackPcol.size();
  for (int i = 0; i < blah; ++i) {
    std::cout << std::endl << "track " << i <<" :";   
    for (auto j : tracksters_near[i]) {
      std::cout << " " << j << " ";
      bool already_in = (std::find(linked_unique.begin(), linked_unique.end(), j) != linked_unique.end());
      if (!already_in) {
        linked_unique.push_back(j);
        const auto &t = trackstersLink[j];
        h_tksters_linked_eta->Fill(t.barycenter().Eta());
        h_eff_phi->Fill(t.barycenter().Phi()); // divided later by h_tksters_tot_phi
        h_eff_pT->Fill(t.raw_pt());
        h_eff_E->Fill(t.raw_energy());
      }
    }
  } // track loop ends

  
}

void TiclDebugger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  static const char* particle_kind[] = {"gam", "e", "mu", "pi0", "h", "h0", "?", "!"};
  using namespace edm;
  using std::begin;
  using std::end;
  using std::iota;
  using std::sort;

  edm::Handle<std::vector<ticl::Trackster>> trackstersMergeH;

  iEvent.getByToken(trackstersMergeToken_, trackstersMergeH);
  auto const& tracksters = *trackstersMergeH.product();

  edm::Handle<std::vector<CaloParticle>> caloParticlesH;
  iEvent.getByToken(caloParticlesToken_, caloParticlesH);
  auto const& caloParticles = *caloParticlesH.product();

  // Linking
  linkTracksters(iEvent, iSetup, tracksters, caloParticlesH);

  std::vector<int> sorted_tracksters_idx(tracksters.size());
  iota(begin(sorted_tracksters_idx), end(sorted_tracksters_idx), 0);
  sort(begin(sorted_tracksters_idx), end(sorted_tracksters_idx), [&tracksters](int i, int j) {
    return tracksters[i].raw_energy() > tracksters[j].raw_energy();
  });

  edm::Handle<std::vector<reco::CaloCluster>> layerClustersH;
  iEvent.getByToken(layerClustersToken_, layerClustersH);
  auto const& layerClusters = *layerClustersH.product();

  edm::Handle<std::vector<reco::Track>> tracksH;
  iEvent.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH.product();

  std::vector<std::pair<int, float>> bestCPMatches;

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
  tree_->Fill();
}

void TiclDebugger::beginRun(edm::Run const&, edm::EventSetup const& es) {
  const CaloGeometry& geom = es.getData(caloGeometry_token_);
  bFieldProd_ = &es.getData(bfield_token_);
  //edm::ESHandle<MagneticField> bfield = es.getHandle(bfield_token_);
  //bFieldProd_ = bfield.product();
  rhtools_.setGeometry(geom);
}

void TiclDebugger::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");
  h_distAll = fs->make<TH1D>("allMinSep", "Trackster-Track closest distance, all", 100, 0., 100.);
  h_distAll->GetXaxis()->SetTitle("cm");
  h_distCaloP = fs->make<TH1D>("CaloPMinSep", "Trackster-Track closest distance, only unconverted", 100, 0., 40.);
  h_distCaloP->GetXaxis()->SetTitle("cm");
  h_NTracks = fs->make<TH1D>("NTracks", "Number of Tracks", 40, 0., 10.);
  //h_NTracks->GetXaxis()->SetTitle("#");
  h_NCaloParticles = fs->make<TH1D>("NCaloParticles", "Number of Caloparticles", 40, 0., 10.);

  // efficiencies
  h_binDiff_track_eta = fs->make<TH1D>("diff_Track_eta", "Bin difference - eta", 20, 0., 5.);
  h_binDiff_track_phi = fs->make<TH1D>("diff_Track_phi", "Bin difference - phi", 20, 0., 5.);
  h_binDiff_track_global = fs->make<TH1D>("diff_Track_glob", "Bin difference - global", 20, 0., 5.);

  h_binDiff_trackster_eta = fs->make<TH1D>("diff_Trackster_eta", "Bin difference - eta", 20, 0., 5.);
  h_binDiff_trackster_phi = fs->make<TH1D>("diff_Trackster_phi", "Bin difference - phi", 20, 0., 5.);
  h_binDiff_trackster_global = fs->make<TH1D>("diff_Trackster_glob", "Bin difference - global", 20, 0., 5.);

  h_tksters_linked_eta = fs->make<TH1D>("trackster_linked_eta", "Trackster linked - eta", 40, -4.0, 4.0);
  h_tksters_tot_eta = fs->make<TH1D>("trackster_eta", "Trackster - eta", 40, -4.0, 4.0);
  h_eff_eta = fs->make<TH1D>("eff_eta", "Efficiency - eta", 40, -4.0, 4.0);
  h_eff_eta->Sumw2();
  h_eff_eta->GetXaxis()->SetTitle("#eta of baryc");

  h_tksters_tot_phi = fs->make<TH1D>("trackster_phi", "Trackster - phi", 40, -1*M_PI, M_PI);
  h_eff_phi = fs->make<TH1D>("eff_phi", "Efficiency - phi", 40, -1*M_PI, M_PI);
  h_eff_phi->Sumw2();
  h_eff_phi->GetXaxis()->SetTitle("#phi of baryc");

  h_tksters_tot_pT = fs->make<TH1D>("trackster_pT", "Trackster - pT", 50, 0, 50);
  h_eff_pT = fs->make<TH1D>("eff_pT", "Efficiency - pT", 50, 0, 50);
  h_eff_pT->Sumw2();
  h_eff_pT->GetXaxis()->SetTitle("pT [GeV]");

  h_tksters_tot_E = fs->make<TH1D>("trackster_E", "Trackster - energy", 50, 0, 500);
  h_eff_E = fs->make<TH1D>("eff_E", "Efficiency - energy", 50, 0, 500);
  h_eff_E->Sumw2();
  h_eff_E->GetXaxis()->SetTitle("Energy [GeV]");

}

void TiclDebugger::endJob() {
  std::cout << std::endl;
  std::cout << std::endl << "Total number of Tracksters = " << totTracksters << std::endl;
  std::cout << "Number of tsos not valid = " << Ntsos_notvalid << std::endl;
  std::cout << "tracks not satisfying cutTk = " << noTracks << std::endl;

  h_eff_eta->Divide(h_tksters_linked_eta,h_tksters_tot_eta);
  h_eff_phi->Divide(h_tksters_tot_phi);
  h_eff_pT->Divide(h_tksters_tot_pT);
  h_eff_E->Divide(h_tksters_tot_E);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TiclDebugger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersMerge", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<std::string>("detector","HGCAL");
  desc.add<std::string>("propagator","PropagatorWithMaterial");
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  descriptions.add("ticlDebugger", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(TiclDebugger);
