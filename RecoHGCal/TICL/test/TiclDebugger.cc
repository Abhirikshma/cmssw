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
#include <array>

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
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

//#include "RecoHGCal/TICL/interface/GlobalCache.h"
#include "RecoHGCal/TICL/plugins/SeedingRegionAlgoBase.h"
#include "RecoHGCal/TICL/plugins/SeedingRegionByTracks.h"
//#include "RecoHGCal/TICL/plugins/TrackstersPCA.h"
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
#include "THStack.h"

#include "PCA_mod.h"

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
  typedef ticl::Trackster::IterationIndex TracksterIterIndex;
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  void linkTracksters(const edm::Event&,
                      const edm::EventSetup&,
                      std::vector<ticl::Trackster>&,
                      std::vector<ticl::Trackster>&,
                      edm::Handle<std::vector<CaloParticle>>);
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
  edm::EDGetTokenT<std::vector<ticl::Trackster>> simTSToken_;
  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> track_col_token_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClustersToken_;
  edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;

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
  TH1D* h_st_distAll;
  TH1D* h_distTS_notLinked;
  TH1D* h_distSTS_notLinked;
  TH1D* h_distCaloP;

  TH1D* h_NTracks;
  TH1D* h_NCaloParticles;

  // TS efficiencies
  TH1D* h_tksters_linked_eta;
  TH1D* h_tksters_tot_eta;
  TH1D* h_eff_eta;
  TH1D* h_eff_phi;
  TH1D* h_tksters_tot_phi;
  TH1D* h_eff_pT;
  TH1D* h_tksters_tot_pT;
  TH1D* h_eff_E;
  TH1D* h_tksters_tot_E;
  // STS
  TH1D* h_stksters_linked_eta;
  TH1D* h_stksters_tot_eta;
  TH1D* h_st_eff_eta;
  TH1D* h_st_eff_phi;
  TH1D* h_stksters_tot_phi;
  TH1D* h_st_eff_pT;
  TH1D* h_stksters_tot_pT;
  TH1D* h_st_eff_E;
  TH1D* h_stksters_tot_E;

  TH1D* h_theta_tk_TS;
  TH1D* h_theta_tk_STS;

  TH1D* h_EfracLinked_TS;
  TH1D* h_EfracLinked_STS;

  TH1D* h_EfracLayer_TS_EM;
  TH1D* h_EfracLayer_TS_HAD;
  TH1D* h_EfracLayer_STS;
  TH1D* h_stdDevLayer_TS_EM;
  TH1D* h_stdDevLayer_TS_HAD;

  TH1D* h_layerOfLargeE_TS_EM;
  TH1D* h_layerOfLargeE_TS_HAD;
  TH1D* h_layerOfLargeE_STS;

  TH1D* h_lowEmaxLayer_STS_E;

  TH2D* h_barycenters_EM;
  TH2D* h_barycenters_HAD;

  THStack* hs_EfracLayer;
  THStack* hs_stdDevLayer_TS;

  int NEvent;

  int totTracksters;
  int totTrackstersPropagated;
  int totsimTracksters;
  int totsimTrackstersPropagated;
  int totTrackstersLinked;
  int totsimTrackstersLinked;
  int totsimTSfromCP;
  int Ntsos_notvalid;
  int noTracks;
  int NSimTracksCrossedBoundary;
  int NsimTS_notValid;

  double totEnergyLinked_TS;
  double totEnergyLinked_STS;

  double layer_Efrac_mean_TS_EM[100];
  double layer_Efrac_mean_TS_HAD[100];
  double layer_Efrac_mean_STS[100];

  double stdDev_by_layer_EM[100];
  double stdDev_by_layer_HAD[100];

  std::vector<int> badPCAs;

  // tracksters selected for propagation
  std::vector<ticl::Trackster> selectedTS;
  std::vector<ticl::Trackster> selectedSTS;
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
      propagator_token_(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", propName_))),
      cutTk_(iConfig.getParameter<std::string>("cutTk")) {
  edm::ConsumesCollector&& iC = consumesCollector();
  trackstersMergeToken_ = iC.consumes<std::vector<ticl::Trackster>>(trackstersMerge_);
  simTSToken_ = iC.consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("label_simTst"));
  tracksToken_ = iC.consumes<std::vector<reco::Track>>(tracks_);
  track_col_token_ = iC.consumes<reco::TrackCollection>(track_col_);
  caloParticlesToken_ = iC.consumes<std::vector<CaloParticle>>(caloParticles_);
  layerClustersToken_ = iC.consumes<std::vector<reco::CaloCluster>>(layerClusters_);
  clustersTime_token_ =
      iC.consumes<edm::ValueMap<std::pair<float, float>>>(iConfig.getParameter<edm::InputTag>("layer_clustersTime"));

  std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
  hdc_token_ = iC.esConsumes<HGCalDDDConstants, IdealGeometryRecord>(edm::ESInputTag("", detectorName_));

  NEvent = 0;

  totTracksters = 0;
  totTrackstersPropagated = 0;
  totsimTracksters = 0;
  totsimTrackstersPropagated = 0;
  totTrackstersLinked = 0;
  totsimTrackstersLinked = 0;
  totsimTSfromCP = 0;
  Ntsos_notvalid = 0;
  noTracks = 0;
  NSimTracksCrossedBoundary = 0;
  NsimTS_notValid = 0;

  for (int i = 0; i < 100; ++i) {
    layer_Efrac_mean_TS_EM[i] = 0;
    layer_Efrac_mean_TS_HAD[i] = 0;
    layer_Efrac_mean_STS[i] = 0;
    stdDev_by_layer_EM[i] = 0;
    stdDev_by_layer_HAD[i] = 0;
  }
}

TiclDebugger::~TiclDebugger() {}

bool TiclDebugger::isCP(const ticl::Trackster& t, edm::Handle<std::vector<CaloParticle>> caloParticlesH) {
  edm::ProductID cp_id = caloParticlesH.id();
  auto tracksterID = t.seedID();
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

void TiclDebugger::linkTracksters(const edm::Event& evt,
                                  const edm::EventSetup& es,
                                  std::vector<ticl::Trackster>& trackstersInput,
                                  std::vector<ticl::Trackster>& simTrackstersInput,
                                  edm::Handle<std::vector<CaloParticle>> caloParticlesH) {
  // distance in cm

  std::cout << "Linking begin" << std::endl;
  edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
  //edm::ESHandle<MagneticField> bfield = es.getHandle(bfield_token_);
  edm::ESHandle<Propagator> propagator = es.getHandle(propagator_token_);

  auto& trackstersLink = trackstersInput;
  auto& simTracksters = simTrackstersInput;

  /*sort(trackstersLink.begin(), trackstersLink.end(), [](const ticl::Trackster &a, const ticl::Trackster &b){
    return a.raw_energy() > b.raw_energy();
  });

  sort(simTracksters.begin(), simTracksters.end(), [](const ticl::Trackster &a, const ticl::Trackster &b){
    return a.raw_energy() > b.raw_energy();
  });*/

  hgcons_ = hdc.product();

  buildFirstLayers();

  edm::Handle<reco::TrackCollection> tracks_h;
  evt.getByToken(track_col_token_, tracks_h);
  //auto bFieldProd = bfield.product();
  const Propagator& prop = (*propagator);

  // propagated point collections
  std::vector<std::pair<std::pair<Vector, Vector>, int>>
      trackPcol;  // propagated track points, errors and index of track in collection
  std::vector<std::pair<Vector, int>> tracksterPcol;     // same but for Tracksters, errors not included here yet
  std::vector<std::pair<Vector, int>> simtracksterPcol;  // simTracksters

  // tiles
  ticl::TICLLayerTile tracksterProp_tile_fw;
  ticl::TICLLayerTile tracksterProp_tile_bw;
  ticl::TICLLayerTile simtracksterProp_tile_fw;
  ticl::TICLLayerTile simtracksterProp_tile_bw;
  ticl::TICLLayerTile trackProp_tile_fw;
  ticl::TICLLayerTile trackProp_tile_bw;

  bool singleSC_zplus = false;
  bool singleSC_zminus = false;
  bool unconverted = false;  //set true if all STS originate from CPs

  totEnergyLinked_TS = 0.;
  totEnergyLinked_STS = 0.;

  // Caloparticle things

  auto const& caloParticles = *caloParticlesH.product();
  std::cout << std::endl << "n caloparticles = " << caloParticles.size() << std::endl;
  h_NCaloParticles->Fill(caloParticles.size());

  double CP_zlpusE = 0.;
  double CP_zminusE = 0.;
  for (auto const& cp : caloParticles) {

    if (cp.g4Tracks()[0].crossedBoundary()) {
      if (cp.g4Tracks()[0].getPositionAtBoundary().Z() > 0) {
        singleSC_zplus = true;
        CP_zlpusE = cp.energy();
      } else if (cp.g4Tracks()[0].getPositionAtBoundary().Z() < 0) {
        singleSC_zminus = true;
        CP_zminusE = cp.energy();
      }
    }

    for (auto const& SimTk : cp.g4Tracks()) {
      if (SimTk.crossedBoundary()) {
        ++NSimTracksCrossedBoundary;
      }
    }
  }
  std::cout << "z - : " << singleSC_zminus << "  z + : " << singleSC_zplus << std::endl;

  // Propagate tracks to the HGCal front

  std::vector<Vector> trackP_mom_fw;
  std::vector<Vector> trackP_mom_bw;
  unsigned nTracks = tracks_h->size();

  h_NTracks->Fill(nTracks);

  for (unsigned i = 0; i < nTracks; ++i) {
    const reco::Track& tk = (*tracks_h)[i];
    if (!cutTk_((tk))) {
      ++noTracks;
      continue;
    }

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd_);
    int iSide = int(tk.eta() > 0);
    TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    std::cout << "propagate!" << std::endl;
    if (tsos.isValid()) {
      std::cout << "x: " << tsos.globalPosition().x() << " y: " << tsos.globalPosition().y()
                << " z: " << tsos.globalPosition().z() << std::endl;
      std::cout << "error : " << tsos.localError().positionError() << std::endl;

      Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
      Vector tkP_momentum(tsos.globalMomentum().x(), tsos.globalMomentum().y(), tsos.globalMomentum().z());
      double x_err = pow(tsos.localError().positionError().xx(), 0.5);
      double y_err = pow(tsos.localError().positionError().yy(), 0.5);
      Vector track_err(x_err, y_err, 0.0);

      trackPcol.push_back(std::make_pair(std::make_pair(trackP, track_err), i));

      if (trackP.Eta() > 0 && singleSC_zplus) {
        trackProp_tile_fw.fill(trackP.Eta(), trackP.Phi(), i);
        trackP_mom_fw.push_back(tkP_momentum);
      }

      else if (trackP.Eta() < 0 && singleSC_zminus) {
        trackProp_tile_bw.fill(trackP.Eta(), trackP.Phi(), i);
        trackP_mom_bw.push_back(tkP_momentum);
      }
    }

    else
      ++Ntsos_notvalid;
  }

  std::cout << "# tracksters = " << trackstersLink.size() << std::endl;
  std::cout << "# simTracksters = " << simTracksters.size() << std::endl;

  // simTrackster things

  int N_from_CP = 0;  // number of simTracksters from CaloParticle
  std::vector<Vector> sts_eigenv_fw;
  std::vector<Vector> sts_eigenv_bw;
  std::vector<unsigned> selectedSTS_idx;
  selectedSTS.clear();

  for (unsigned i = 0; i < simTracksters.size(); ++i) {
    std::cout << "simtracksters loop" << std::endl;

    const auto& st = simTracksters[i];

    if (isCP(st, caloParticlesH))
      ++N_from_CP;

    Vector baryc = st.barycenter();
    Vector directnv = st.eigenvectors(0);

    if (abs(directnv.Z()) < 0.00001) {
      std::cout << "eigenvec = " << directnv << std::endl;
      continue;
    }

    bool zplus_reject = true;
    bool zminus_reject = true;

    if (baryc.Z() > 0 && singleSC_zplus) {
      if (st.raw_energy() > 0. * CP_zlpusE)
        zplus_reject = false;
    }

    else if (baryc.Z() < 0 && singleSC_zminus) {
      if (st.raw_energy() > 0. * CP_zminusE)
        zminus_reject = false;
    }

    // Propagate simTracksters to the HGCal front

    //if ((baryc.Z() > 0 && singleSC_zplus) || (baryc.Z() < 0 && singleSC_zminus)) {
    if (!(zplus_reject && zminus_reject)) {
      //if (zplus_reject && zminus_reject)  continue;
      selectedSTS.push_back(st);

      float zVal = hgcons_->waferZ(1, true);
      zVal *= (baryc.Z() > 0) ? 1 : -1;

      double par = (zVal - baryc.Z()) / directnv.Z();

      double xOnSurface = par * directnv.X() + baryc.X();
      double yOnSurface = par * directnv.Y() + baryc.Y();

      std::cout << "x = " << xOnSurface << std::endl;
      std::cout << "y = " << yOnSurface << std::endl;
      std::cout << "z = " << zVal << std::endl;

      Vector simtracksterP(xOnSurface, yOnSurface, zVal);
      selectedSTS_idx.push_back(i);
      simtracksterPcol.push_back(std::make_pair(simtracksterP, i));

      //simtracksters that can be linked
      h_stksters_tot_eta->Fill(st.barycenter().Eta());
      h_stksters_tot_phi->Fill(st.barycenter().Phi());
      h_stksters_tot_pT->Fill(st.raw_pt());
      h_stksters_tot_E->Fill(st.raw_energy());
      if (!isCP(st, caloParticlesH))
        std::cout << "wrong selection" << std::endl;

      if (simtracksterP.Eta() > 0) {
        simtracksterProp_tile_fw.fill(simtracksterP.Eta(), simtracksterP.Phi(), i);
        sts_eigenv_fw.push_back(directnv);
      }

      else if (simtracksterP.Eta() < 0) {
        simtracksterProp_tile_bw.fill(simtracksterP.Eta(), simtracksterP.Phi(), i);
        sts_eigenv_bw.push_back(directnv);
      }

      //if (st.vertices().size() > 0) totsimTracksters++;
    }

  }  // STS loop end

  totsimTracksters += simTracksters.size();
  totsimTrackstersPropagated += simtracksterPcol.size();
  totsimTSfromCP += N_from_CP;

  std::cout << std::endl << N_from_CP << " from CPs " << std::endl;
  if (N_from_CP == static_cast<int>(simTracksters.size()))
    unconverted = true;

  // Propagate tracksters (barycenters) to the HGCal front
  selectedTS.clear();

  std::vector<Vector> ts_eigenv_fw;
  std::vector<Vector> ts_eigenv_bw;
  std::vector<unsigned> selectedTS_idx;
  for (unsigned i = 0; i < trackstersLink.size(); ++i) {
    std::cout << "tracksters loop" << std::endl;
    const auto& t = trackstersLink[i];

    std::cout << "E = " << t.raw_energy() << std::endl;
  
    Vector baryc = t.barycenter();
    Vector directnv = t.eigenvectors(0);
    Vector trackster_err(pow(t.sigmas()[0], 0.5), pow(t.sigmas()[1], 0.5), 0.0);

    bool zplus_reject = true;
    bool zminus_reject = true;

    if (baryc.Z() > 0 && singleSC_zplus) {
      if (t.raw_energy() > 0.1 * CP_zlpusE)
        zplus_reject = false;
    }

    else if (baryc.Z() < 0 && singleSC_zminus) {
      if (t.raw_energy() > 0.1 * CP_zminusE)
        zminus_reject = false;
    }

    //if ((baryc.Z() > 0 && singleSC_zplus) || (baryc.Z() < 0 && singleSC_zminus)) {
    if (!(zplus_reject && zminus_reject)) {
      selectedTS.push_back(t);

      float zVal = hgcons_->waferZ(1, true);
      zVal *= (baryc.Z() > 0) ? 1 : -1;

      double par = (zVal - baryc.Z()) / directnv.Z();

      double xOnSurface = par * directnv.X() + baryc.X();
      double yOnSurface = par * directnv.Y() + baryc.Y();

      std::cout << "x = " << xOnSurface << std::endl;
      std::cout << "y = " << yOnSurface << std::endl;
      std::cout << "z = " << zVal << std::endl;
      std::cout << "error :  " << t.sigmas()[0] << "  " << t.sigmas()[1] << "  " << t.sigmas()[2] << std::endl;
      //std::cout << "PCA err :  " << t.sigmasPCA()[0] << "  " << t.sigmasPCA()[1] << "  " << t.sigmasPCA()[2] << std::endl;

      Vector tracksterP(xOnSurface, yOnSurface, zVal);
      selectedTS_idx.push_back(i);
      tracksterPcol.push_back(std::make_pair(tracksterP, i));

      //tracksters than can be linked
      h_tksters_tot_eta->Fill(t.barycenter().Eta());
      h_tksters_tot_phi->Fill(t.barycenter().Phi());
      h_tksters_tot_pT->Fill(t.raw_pt());
      h_tksters_tot_E->Fill(t.raw_energy());

      if (tracksterP.Eta() > 0) {
        tracksterProp_tile_fw.fill(tracksterP.Eta(), tracksterP.Phi(), i);
        ts_eigenv_fw.push_back(directnv);
      }

      else if (tracksterP.Eta() < 0) {
        tracksterProp_tile_bw.fill(tracksterP.Eta(), tracksterP.Phi(), i);
        ts_eigenv_bw.push_back(directnv);
      }
    }
  }

  totTracksters += trackstersLink.size();
  totTrackstersPropagated += tracksterPcol.size();

  // Search box over tiles + preliminary linking

  std::vector<unsigned> tracksters_near
      [trackPcol.size()];  // vector element for every track contains the indices of tracksters near to that track
  std::vector<unsigned> simtracksters_near[trackPcol.size()];
  int nth = 0;
  for (auto i : trackPcol) {
    auto trackP = i.first.first;
    //auto track_err = i.first.second;

    double tk_eta = trackP.Eta();
    double tk_phi = trackP.Phi();

    double delta = 0.02;

    if (tk_eta > 0) {
      /*int etaBin = tracksterProp_tile_fw.etaBin(tk_eta);
      int phiBin = tracksterProp_tile_fw.phiBin(tk_phi);

      int etaBin_st = simtracksterProp_tile_fw.etaBin(tk_eta); // same as bins from TS tile
      int phiBin_st = simtracksterProp_tile_fw.phiBin(tk_phi);


      // deltas from the uncertainties in track propagation
      int delta_eta = abs(trackProp_tile_fw.etaBin((trackP-track_err).Eta()) - trackProp_tile_fw.etaBin((trackP+track_err).Eta()));
      int delta_phi = abs(trackProp_tile_fw.phiBin((trackP-track_err).Phi()) - trackProp_tile_fw.phiBin((trackP+track_err).Phi()));*/

      double eta_min = ((tk_eta - delta) > 0) ? (tk_eta - delta) : 0;

      std::array<int, 4> search_box =
          tracksterProp_tile_fw.searchBoxEtaPhi(eta_min, tk_eta + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }

      //std::array<int, 4> search_box = {{etaBin - 1, etaBin + 1, phiBin - 1, phiBin + 1}}; // 3x3
      //std::array<int, 4> search_box = {{etaBin - delta_eta, etaBin + delta_eta, phiBin - delta_phi, phiBin + delta_phi}}; // f(uncertainties)

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto tracksters_in_box = tracksterProp_tile_fw[tracksterProp_tile_fw.globalBin(eta_i, phi_i)];
          tracksters_near[nth].insert(
              std::end(tracksters_near[nth]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // search in box ends, TS

      std::array<int, 4> search_box_st =
          simtracksterProp_tile_fw.searchBoxEtaPhi(eta_min, tk_eta + delta, tk_phi - delta, tk_phi + delta);
      if (search_box_st[2] > search_box_st[3]) {
        double temp = search_box_st[3];
        search_box_st[3] = search_box_st[2];
        search_box_st[2] = temp;
      }
      //std::array<int, 4> search_box_st = {{etaBin_st - 1, etaBin_st + 1, phiBin_st - 1, phiBin_st + 1}}; // 3x3
      //std::array<int, 4> search_box_st = {{etaBin_st - delta_eta, etaBin_st + delta_eta, phiBin_st - delta_phi, phiBin_st + delta_phi}}; // f(uncertainties)

      for (int eta_i = search_box_st[0]; eta_i <= search_box_st[1]; ++eta_i) {
        for (int phi_i = search_box_st[2]; phi_i <= search_box_st[3]; ++phi_i) {
          auto simtracksters_in_box = simtracksterProp_tile_fw[simtracksterProp_tile_fw.globalBin(eta_i, phi_i)];
          simtracksters_near[nth].insert(
              std::end(simtracksters_near[nth]), std::begin(simtracksters_in_box), std::end(simtracksters_in_box));
        }
      }  // search in box ends, STS
    }

    if (tk_eta < 0) {
      double eta_min = ((abs(tk_eta) - delta) > 0) ? (abs(tk_eta) - delta) : 0;

      std::array<int, 4> search_box =
          tracksterProp_tile_bw.searchBoxEtaPhi(eta_min, abs(tk_eta) + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }
      
      //std::array<int, 4> search_box = {{etaBin - 1, etaBin + 1, phiBin - 1, phiBin + 1}}; // 3x3
      //std::array<int, 4> search_box = {{etaBin - delta_eta, etaBin + delta_eta, phiBin - delta_phi, phiBin + delta_phi}}; // f(uncertainties)

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto tracksters_in_box = tracksterProp_tile_bw[tracksterProp_tile_bw.globalBin(eta_i, phi_i)];
          tracksters_near[nth].insert(
              std::end(tracksters_near[nth]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // search in box ends, TS

      std::array<int, 4> search_box_st =
          simtracksterProp_tile_bw.searchBoxEtaPhi(eta_min, abs(tk_eta) + delta, tk_phi - delta, tk_phi + delta);
      if (search_box_st[2] > search_box_st[3]) {
        double temp = search_box_st[3];
        search_box_st[3] = search_box_st[2];
        search_box_st[2] = temp;
      }
      //std::array<int, 4> search_box_st = {{etaBin_st - 1, etaBin_st + 1, phiBin_st - 1, phiBin_st + 1}}; // 3x3
      //std::array<int, 4> search_box_st = {{etaBin_st - delta_eta, etaBin_st + delta_eta, phiBin_st - delta_phi, phiBin_st + delta_phi}}; // f(uncertainties)

      for (int eta_i = search_box_st[0]; eta_i <= search_box_st[1]; ++eta_i) {
        for (int phi_i = search_box_st[2]; phi_i <= search_box_st[3]; ++phi_i) {
          auto simtracksters_in_box = simtracksterProp_tile_bw[simtracksterProp_tile_bw.globalBin(eta_i, phi_i)];
          simtracksters_near[nth].insert(
              std::end(simtracksters_near[nth]), std::begin(simtracksters_in_box), std::end(simtracksters_in_box));
        }
      }  // search in box ends, STS
    }

    ++nth;
  }  // track loop ends

  // Find uniquely linked TS + fill histos etc. for validation

  std::vector<unsigned> linked_unique;
  std::vector<unsigned> linked_unique_st;
  int NtrackPcol = trackPcol.size();

  std::cout << std::endl << "search results (track : tracksters): " << std::endl;

  for (int i = 0; i < NtrackPcol; ++i) {
    std::cout << std::endl << "track " << i << " :";
    for (auto j : tracksters_near[i]) {
      if (tracksters_near[i].empty())
        continue;
      std::cout << " " << j << " ";
      bool already_in = (std::find(linked_unique.begin(), linked_unique.end(), j) != linked_unique.end());
      if (!already_in) {
        linked_unique.push_back(j);
        const auto& t = trackstersLink[j];
        totEnergyLinked_TS += t.raw_energy();

        // tracksters that are linked
        h_tksters_linked_eta->Fill(t.barycenter().Eta());
        h_eff_phi->Fill(t.barycenter().Phi());  // divided later by h_tksters_tot_phi
        h_eff_pT->Fill(t.raw_pt());
        h_eff_E->Fill(t.raw_energy());
      }
    }  // tracksters_near loop ends

    for (auto k : simtracksters_near[i]) {
      if (simtracksters_near[i].empty())
        continue;
      bool already_in_st = (std::find(linked_unique_st.begin(), linked_unique_st.end(), k) != linked_unique_st.end());
      if (!already_in_st) {
        linked_unique_st.push_back(k);
        const auto& st = simTracksters[k];
        totEnergyLinked_STS += st.raw_energy();

        // simtracksters that are linked
        h_stksters_linked_eta->Fill(st.barycenter().Eta());
        h_st_eff_phi->Fill(st.barycenter().Phi());  // divided later
        h_st_eff_pT->Fill(st.raw_pt());
        h_st_eff_E->Fill(st.raw_energy());
      }
    }  // simtracksters_near loop ends

  }  // trackPcol loop ends

  totTrackstersLinked += linked_unique.size();
  totsimTrackstersLinked += linked_unique_st.size();

  if (CP_zlpusE + CP_zminusE > 0) {
    h_EfracLinked_TS->Fill(totEnergyLinked_TS / (CP_zlpusE + CP_zminusE));
    h_EfracLinked_STS->Fill(totEnergyLinked_STS / (CP_zlpusE + CP_zminusE));
  }

  std::cout << std::endl << "TSs linked:  ";
  for (int i : linked_unique) {
    std::cout << i << " ";
  }
  std::cout << std::endl << "STSs linked:  ";
  for (int i : linked_unique_st) {
    std::cout << i << " ";
  }
  std::cout << std::endl;

  std::vector<int> closestTracks(tracksterPcol.size(),
                                 -1);  // index of closest track for every trackster, arranged as in tracksterPcol
  nth = 0;

  // Calculating TS/STS-trk separations + finding closest tracks to tracksters

  for (auto const i : tracksterPcol) {
    double minSep = 600.;  // prevent trks from other side from being assigned as closest
    for (auto const j : trackPcol) {
      double sep = pow((i.first - j.first.first).Mag2(), 0.5);
      if (sep < minSep) {
        minSep = sep;
        closestTracks[nth] = j.second;
      }
    }
    if (minSep < 600) {
      h_distAll->Fill(minSep);
    }
    if (unconverted) {
      h_distCaloP->Fill(minSep);
    }
    if (std::find(linked_unique.begin(), linked_unique.end(), i.second) ==
        linked_unique.end()) {  // propagated but not linked
      h_distTS_notLinked->Fill(minSep);
    }

    ++nth;
  }  // tracksterPcol loop end

  for (auto const i : simtracksterPcol) {
    double minSep = 600.;  // prevent trks from other side from being assigned as closest
    for (auto const j : trackPcol) {
      double sep = pow((i.first - j.first.first).Mag2(), 0.5);
      if (sep < minSep) {
        minSep = sep;
      }
    }
    if (minSep < 600) {
      h_st_distAll->Fill(minSep);
    }

    if (std::find(linked_unique_st.begin(), linked_unique_st.end(), i.second) ==
        linked_unique_st.end()) {  // propagated but not linked
      h_distSTS_notLinked->Fill(minSep);
    }

  }  // simtracksterPcol loop end

  std::cout << "closest tracks to tracksters: ";
  for (int i : closestTracks) {
    std::cout << i << " ";
  }
  std::cout << std::endl;

  // e_tk dot e_PCA for S/TS
  for (auto trk_p : trackP_mom_fw) {  //fwd
    auto trkV = trk_p.Unit();
    for (auto v : ts_eigenv_fw) {  //TS
      double cos_theta = trkV.Dot(v);
      h_theta_tk_TS->Fill(abs(cos_theta));
      if (abs(cos_theta) < 0.8)
        badPCAs.push_back(NEvent);
    }
    for (auto v : sts_eigenv_fw) {  //STS
      double cos_theta = trkV.Dot(v);
      h_theta_tk_STS->Fill(abs(cos_theta));
      if (abs(cos_theta) < 0.8)
        badPCAs.push_back(NEvent);
    }
  }
  for (auto trk_p : trackP_mom_bw) {  //bwd
    auto trkV = trk_p.Unit();
    for (auto v : ts_eigenv_bw) {  //TS
      double cos_theta = trkV.Dot(v);
      h_theta_tk_TS->Fill(abs(cos_theta));
      if (abs(cos_theta) < 0.8)
        badPCAs.push_back(NEvent);
    }
    for (auto v : sts_eigenv_bw) {  //STS
      double cos_theta = trkV.Dot(v);
      h_theta_tk_STS->Fill(abs(cos_theta));
      if (abs(cos_theta) < 0.8)
        badPCAs.push_back(NEvent);
    }
  }

  std::cout << "Linking done" << std::endl;

}  // linkTracksters

void TiclDebugger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  ++NEvent;
  std::cout << std::endl << std::endl << "-----EVENT " << NEvent << "-----" << std::endl;

  static const char* particle_kind[] = {"gam", "e", "mu", "pi0", "h", "h0", "?", "!"};
  using namespace edm;
  using std::begin;
  using std::end;
  using std::iota;
  using std::sort;

  edm::Handle<std::vector<ticl::Trackster>> trackstersMergeH;
  iEvent.getByToken(trackstersMergeToken_, trackstersMergeH);
  auto const& tracksters = *trackstersMergeH.product();

  edm::Handle<std::vector<ticl::Trackster>> simTS_h;
  iEvent.getByToken(simTSToken_, simTS_h);
  auto const& simTracksters = *simTS_h.product();

  // local copies for the modified PCA
  std::vector<ticl::Trackster> tracksters_ = tracksters;
  auto& tracksters_mutable = tracksters_;

  std::vector<ticl::Trackster> simTracksters_ = simTracksters;
  auto& simTracksters_mutable = simTracksters_;

  edm::Handle<std::vector<CaloParticle>> caloParticlesH;
  iEvent.getByToken(caloParticlesToken_, caloParticlesH);
  auto const& caloParticles = *caloParticlesH.product();

  edm::Handle<std::vector<reco::CaloCluster>> layerClustersH;
  iEvent.getByToken(layerClustersToken_, layerClustersH);
  auto const& layerClusters = *layerClustersH.product();

  edm::Handle<edm::ValueMap<std::pair<float, float>>> clustersTimeH;
  iEvent.getByToken(clustersTime_token_, clustersTimeH);
  const auto& layerClustersTimes = *clustersTimeH;

  std::cout << "N tracksters " << tracksters_.size() << std::endl;
  std::cout << "N simTracksters " << simTracksters_.size() << std::endl;
  // PCA with cleaning
  assignPCAtoTracksters_mod(tracksters_mutable,
                            layerClusters,
                            layerClustersTimes,
                            rhtools_,
                            rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z());

  assignPCAtoTracksters_mod(simTracksters_mutable,
                            layerClusters,
                            layerClustersTimes,
                            rhtools_,
                            rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z());

  // Linking
  linkTracksters(iEvent, iSetup, tracksters_mutable, simTracksters_mutable, caloParticlesH);


  // Trackster diagonistics, by layer
  // TS
  int nonZero_trackster = 0;
  int nonZero_EM = 0;
  int nonZero_HAD = 0;
  std::array<double, 100> layer_Efrac_TS_EM = {};
  std::array<double, 100> layer_Efrac_TS_HAD = {};
  double stdDev_by_layer_EM_[100] = {};
  double stdDev_by_layer_HAD_[100] = {};

  for (auto& t : tracksters) {
    unsigned N = t.vertices().size();
    if (N == 0)
      continue;
    ++nonZero_trackster;
    bool EM = (t.ticlIteration() == TracksterIterIndex::EM);
    bool HAD = (t.ticlIteration() == TracksterIterIndex::HAD);

    Vector bary = t.barycenter();

    std::vector<double> LCenergies;
    for (unsigned i = 0; i < N; ++i) {
      LCenergies.push_back(layerClusters[t.vertices(i)].energy());
    }
    auto maxE_vert = std::distance(LCenergies.begin(), std::max_element(LCenergies.begin(), LCenergies.end()));
    auto maxE_layer = rhtools_.getLayerWithOffset(layerClusters[t.vertices(maxE_vert)].hitsAndFractions()[0].first);
    if (EM) {
      h_layerOfLargeE_TS_EM->Fill(maxE_layer);
      h_barycenters_EM->Fill(bary.Z(), pow(bary.Perp2(), 0.5) );
      ++nonZero_EM;
    }
    else if (HAD) {
      h_layerOfLargeE_TS_HAD->Fill(maxE_layer);
      h_barycenters_HAD->Fill(bary.Z(), pow(bary.Perp2(), 0.5) );
      ++nonZero_HAD;
    }

    auto vertices_by_layer = sortByLayer(t, layerClusters, rhtools_);
    double maxElayer_energy = 0;
    for (auto v : vertices_by_layer[maxE_layer])
      maxElayer_energy += layerClusters[t.vertices(v)].energy();

    auto vertT_by_layer = translateTrackster(t, layerClusters, rhtools_, maxE_layer);

    // energy weighted standard dev of LC by layer
    for (int i = 0; i < 100; ++i) {
      auto vertices_in_layer = vertT_by_layer[i];
      if (vertices_in_layer.size() <= 1)
        continue;

      auto bary = barycenterInLayer(t, layerClusters, vertices_in_layer);

      double numr = 0.;
      double denr = 0.;
      for (auto v : vertices_in_layer) {
        auto thisLC = layerClusters[t.vertices(v)];
        Vector coords(thisLC.x(), thisLC.y(), thisLC.z());
        double wt = thisLC.energy();
        double diff_sq = (coords - bary).Mag2();
        numr += wt * diff_sq;
        denr += wt;
      }
      double M = (vertices_in_layer.size() - 1) / (double)vertices_in_layer.size();
      double stdDev = pow(numr / (M * denr), 0.5);

      if (EM)
      stdDev_by_layer_EM_[i] += stdDev;
      else if (HAD)
      stdDev_by_layer_HAD_[i] += stdDev;
    }

    // E frac by layer
    for (int i = 0; i < 100; ++i) {
      auto vertices_in_layer = vertT_by_layer[i];
      //if (vertices_in_layer.empty()) continue;
      double e_layer = 0.;
      for (auto v : vertices_in_layer)
        e_layer += layerClusters[t.vertices(v)].energy();
      //e_layer /= maxEnergy;
      e_layer /= maxElayer_energy;
      if (EM)
      layer_Efrac_TS_EM[i] += e_layer;
      if (HAD)
      layer_Efrac_TS_HAD[i] += e_layer;
    }
  } // Tracksters

  if (nonZero_EM > 0)
    for (int i = 0; i < 100; ++i) {
      layer_Efrac_TS_EM[i] /= nonZero_EM;    
      stdDev_by_layer_EM[i] += (stdDev_by_layer_EM_[i] / nonZero_EM);
    }
  if (nonZero_HAD > 0)
    for (int i = 0; i < 100; ++i) {
      layer_Efrac_TS_HAD[i] /= nonZero_HAD;
      stdDev_by_layer_HAD[i] += (stdDev_by_layer_HAD_[i] / nonZero_HAD);
    }

  // STS
  nonZero_trackster = 0;
  std::array<double, 100> layer_Efrac_STS_ = {};
  for (auto& t : simTracksters) {
    unsigned N = t.vertices().size();
    if (N == 0)
      continue;
    ++nonZero_trackster;

    auto vertices_by_layer = sortByLayer(t, layerClusters, rhtools_);

    std::vector<double> LCenergies;
    for (unsigned i = 0; i < N; ++i) {
      LCenergies.push_back(layerClusters[t.vertices(i)].energy());
    }

    auto maxE_vert = std::distance(LCenergies.begin(), std::max_element(LCenergies.begin(), LCenergies.end()));
    auto maxE_layer = rhtools_.getLayerWithOffset(layerClusters[t.vertices(maxE_vert)].hitsAndFractions()[0].first);
    h_layerOfLargeE_STS->Fill(maxE_layer);
    double maxElayer_energy = 0;
    for (auto v : vertices_by_layer[maxE_layer])
      maxElayer_energy += layerClusters[t.vertices(v)].energy();

    if ((int)maxE_layer <= 7)
      h_lowEmaxLayer_STS_E->Fill(t.raw_energy());

    auto vertT_by_layer = translateTrackster(t, layerClusters, rhtools_, maxE_layer);

    // E frac by layer
    for (unsigned i = 0; i < 100; ++i) {
      auto vertices_in_layer = vertT_by_layer[i];
      //if (vertices_in_layer.empty()) continue;
      double e_layer = 0.;
      for (auto v : vertices_in_layer)
        e_layer += layerClusters[t.vertices(v)].energy();
      //e_layer /= maxEnergy;
      e_layer /= maxElayer_energy;
      layer_Efrac_STS_[i] += e_layer;
    }
  } // simTracksters

  if (nonZero_trackster > 0)
    for (int i = 0; i < 100; ++i) {
      layer_Efrac_STS_[i] /= nonZero_trackster;
    }

  for (int i = 0; i < 100; ++i) {
    layer_Efrac_mean_TS_EM[i] += layer_Efrac_TS_EM[i];
    layer_Efrac_mean_TS_HAD[i] += layer_Efrac_TS_HAD[i];
    layer_Efrac_mean_STS[i] += layer_Efrac_STS_[i];
  }

  std::vector<int> sorted_tracksters_idx(tracksters.size());
  iota(begin(sorted_tracksters_idx), end(sorted_tracksters_idx), 0);
  sort(begin(sorted_tracksters_idx), end(sorted_tracksters_idx), [&tracksters](int i, int j) {
    return tracksters[i].raw_energy() > tracksters[j].raw_energy();
  });

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

  // distances
  h_distAll = fs->make<TH1D>("allMinSep", "Trackster-Track closest distance, all", 100, 0., 100.);
  h_distAll->GetXaxis()->SetTitle("cm");
  h_distTS_notLinked = fs->make<TH1D>("ts_notLinked", "not linked TS-Tk", 50, 0., 10.);
  h_distTS_notLinked->GetXaxis()->SetTitle("cm");
  h_st_distAll = fs->make<TH1D>("st_allMinSep", "simTrackster-Track closest distance, all", 100, 0., 100.);
  h_st_distAll->GetXaxis()->SetTitle("cm");
  h_distSTS_notLinked = fs->make<TH1D>("sts_notLinked", "not linked STS-Tk", 100, 0., 40.);
  h_distSTS_notLinked->GetXaxis()->SetTitle("cm");
  h_distCaloP = fs->make<TH1D>("CaloPMinSep", "Trackster-Track closest distance, only unconverted", 100, 0., 40.);
  h_distCaloP->GetXaxis()->SetTitle("cm");

  h_NTracks = fs->make<TH1D>("NTracks", "Number of Tracks", 40, 0., 10.);
  //h_NTracks->GetXaxis()->SetTitle("#");
  h_NCaloParticles = fs->make<TH1D>("NCaloParticles", "Number of Caloparticles", 40, 0., 10.);

  // efficiencies
  h_tksters_linked_eta = fs->make<TH1D>("trackster_linked_eta", "Trackster linked - eta", 40, -4.0, 4.0);
  h_tksters_tot_eta = fs->make<TH1D>("trackster_eta", "Trackster - eta", 40, -4.0, 4.0);
  h_eff_eta = fs->make<TH1D>("eff_eta", "Efficiency - eta", 40, -4.0, 4.0);
  h_eff_eta->Sumw2();
  h_eff_eta->GetXaxis()->SetTitle("#eta of baryc");

  h_tksters_tot_phi = fs->make<TH1D>("trackster_phi", "Trackster - phi", 20, -1 * M_PI, M_PI);
  h_eff_phi = fs->make<TH1D>("eff_phi", "Efficiency - phi", 20, -1 * M_PI, M_PI);
  h_eff_phi->Sumw2();
  h_eff_phi->GetXaxis()->SetTitle("#phi of baryc");

  h_tksters_tot_pT = fs->make<TH1D>("trackster_pT", "Trackster - pT", 50, 0, 50);
  h_eff_pT = fs->make<TH1D>("eff_pT", "Efficiency - pT", 50, 0, 50);
  h_eff_pT->Sumw2();
  h_eff_pT->GetXaxis()->SetTitle("pT [GeV]");

  h_tksters_tot_E = fs->make<TH1D>("trackster_E", "Trackster - energy", 30, 0, 300);
  h_eff_E = fs->make<TH1D>("eff_E", "Efficiency - energy", 30, 0, 300);
  h_eff_E->Sumw2();
  h_eff_E->GetXaxis()->SetTitle("Energy [GeV]");

  // efficiencies for STS linking
  h_stksters_linked_eta = fs->make<TH1D>("simtrackster_linked_eta", "simTrackster linked - eta", 40, -4.0, 4.0);
  h_stksters_tot_eta = fs->make<TH1D>("simtrackster_eta", "simTrackster - eta", 40, -4.0, 4.0);
  h_st_eff_eta = fs->make<TH1D>("st_eff_eta", "Efficiency - eta", 40, -4.0, 4.0);
  h_st_eff_eta->Sumw2();
  h_st_eff_eta->GetXaxis()->SetTitle("#eta of baryc");

  h_stksters_tot_phi = fs->make<TH1D>("simtrackster_phi", "simTrackster - phi", 20, -1 * M_PI, M_PI);
  h_st_eff_phi = fs->make<TH1D>("st_eff_phi", "Efficiency - phi", 20, -1 * M_PI, M_PI);
  h_st_eff_phi->Sumw2();
  h_st_eff_phi->GetXaxis()->SetTitle("#phi of baryc");

  h_stksters_tot_pT = fs->make<TH1D>("simtrackster_pT", "simTrackster - pT", 50, 0, 50);
  h_st_eff_pT = fs->make<TH1D>("st_eff_pT", "Efficiency - pT", 50, 0, 50);
  h_st_eff_pT->Sumw2();
  h_st_eff_pT->GetXaxis()->SetTitle("pT [GeV]");

  h_stksters_tot_E = fs->make<TH1D>("simtrackster_E", "simTrackster - energy", 30, 0, 300);
  h_st_eff_E = fs->make<TH1D>("st_eff_E", "Efficiency - energy", 30, 0, 300);
  h_st_eff_E->Sumw2();
  h_st_eff_E->GetXaxis()->SetTitle("Energy [GeV]");

  // e_tk dot e_PCA
  h_theta_tk_STS = fs->make<TH1D>("tk_STS", "#hat{e}_{tk} #bullet #hat{e}_{PCA} - STS", 60, .7, 1.0);
  h_theta_tk_TS = fs->make<TH1D>("tk_TS", "#hat{e}_{tk} #bullet #hat{e}_{PCA} - TS", 60, .7, 1.0);

  // E fraction linked
  h_EfracLinked_TS = fs->make<TH1D>("Efrac_TS", "E frac linked - TS", 20, 0.0, 1.0);
  h_EfracLinked_TS->GetXaxis()->SetTitle("E_{TS}^{link}/E_{CP}^{sel}");
  h_EfracLinked_STS = fs->make<TH1D>("Efrac_STS", "E frac linked - STS", 20, 0.0, 1.0);
  h_EfracLinked_STS->GetXaxis()->SetTitle("E_{STS}^{link}/E_{CP}^{sel}");

  // E frac in layer
  h_EfracLayer_TS_EM = fs->make<TH1D>("EfracLayer_TS_EM", "mean energy in layer - TS", 100, 0.5, 100.5);
  h_EfracLayer_TS_EM->GetXaxis()->SetTitle("layer number");
  h_EfracLayer_TS_EM->GetYaxis()->SetTitle("frac E in layer wrt max E layer");
  h_EfracLayer_TS_HAD = fs->make<TH1D>("EfracLayer_TS_HAD", "mean energy in layer - TS", 100, 0.5, 100.5);
  h_EfracLayer_TS_HAD->GetXaxis()->SetTitle("layer number");
  h_EfracLayer_TS_HAD->GetYaxis()->SetTitle("frac E in layer wrt max E layer");

  h_EfracLayer_STS = fs->make<TH1D>("EfracLayer_STS", "mean energy in layer - STS", 100, 0.5, 100.5);
  h_EfracLayer_STS->GetXaxis()->SetTitle("layer number");
  h_EfracLayer_STS->GetYaxis()->SetTitle("frac E in layer wrt max E layer");

  // layer of max E LC of shower
  h_layerOfLargeE_TS_EM = fs->make<TH1D>("maxE_layer_TS_sel_EM", "max E LC layer", 50, 0.5, 50.5);
  h_layerOfLargeE_TS_EM->GetXaxis()->SetTitle("layer number");
  h_layerOfLargeE_TS_HAD = fs->make<TH1D>("maxE_layer_TS_sel_HAD", "max E LC layer", 50, 0.5, 50.5);
  h_layerOfLargeE_TS_HAD->GetXaxis()->SetTitle("layer number");
  h_layerOfLargeE_STS = fs->make<TH1D>("maxE_layer_STS_onlyCP", "max E LC layer", 50, 0.5, 50.5);
  h_layerOfLargeE_STS->GetXaxis()->SetTitle("layer number");

  h_lowEmaxLayer_STS_E = fs->make<TH1D>("lowEmaxLayer", "STS with max E LC layer <= 7", 50, 0, 50);
  h_lowEmaxLayer_STS_E->GetXaxis()->SetTitle("Energy [GeV]");

  h_stdDevLayer_TS_EM = fs->make<TH1D>("stdDevLayer_TS_EM", "", 100, 0.5, 100.5);
  h_stdDevLayer_TS_EM->GetYaxis()->SetTitle("E weighted std dev of LCs/layer");
  h_stdDevLayer_TS_EM->GetXaxis()->SetTitle("layer number (max E layer -> 50)");
  h_stdDevLayer_TS_HAD = fs->make<TH1D>("stdDevLayer_TS_HAD", "", 100, 0.5, 100.5);
  h_stdDevLayer_TS_HAD->GetYaxis()->SetTitle("E weighted std dev of LCs/layer");
  h_stdDevLayer_TS_HAD->GetXaxis()->SetTitle("layer number (max E layer -> 50)");

  // barycenters
  // (x,y) in histo -> (Z,R) in det 
  h_barycenters_EM = fs->make<TH2D>("barycenters_EM","Barycenters EM",45,300,525,50,25,275);
  h_barycenters_EM->GetYaxis()->SetTitle("R [cm]");
  h_barycenters_EM->GetXaxis()->SetTitle("Z [cm]");
  h_barycenters_HAD = fs->make<TH2D>("barycenters_HAD","Barycenters HAD",45,300,525,50,25,275);
  h_barycenters_HAD->GetYaxis()->SetTitle("R [cm]");
  h_barycenters_HAD->GetXaxis()->SetTitle("Z [cm]");

  hs_EfracLayer = fs->make<THStack>("s1", "frac E in layer wrt max E layer");
  hs_stdDevLayer_TS = fs->make<THStack>("s2", "E weighted std dev of LCs/layer");
}

void TiclDebugger::endJob() {
  std::cout << std::endl;
  std::cout << std::endl << "Total number of Tracksters = " << totTracksters << std::endl;
  std::cout << "Total number of Tracksters propagated = " << totTrackstersPropagated << std::endl;
  std::cout << "Total number of simTracksters = " << totsimTracksters << std::endl;
  std::cout << "Total number of simTracksters propagated = " << totsimTrackstersPropagated << std::endl;
  std::cout << "TS linking - glob efficiency = " << totTrackstersLinked / (double)totTrackstersPropagated << std::endl;
  std::cout << "STS linking - glob efficiency = " << totsimTrackstersLinked / (double)totsimTrackstersPropagated
            << std::endl;
  std::cout << "Number of simTrackster handles not valid = " << NsimTS_notValid << std::endl;
  std::cout << "sim TS from CP = " << totsimTSfromCP << std::endl;
  std::cout << "Number of tsos not valid = " << Ntsos_notvalid << std::endl;
  std::cout << "tracks not satisfying cutTk = " << noTracks << std::endl;
  std::cout << "Total no. of SimTracks crossed boundary = " << NSimTracksCrossedBoundary << std::endl;
  std::cout << "Bad events (PCA) : ";
  for (auto i : badPCAs)
    std::cout << i << " ";
  std::cout << std::endl;

  h_eff_eta->Divide(h_tksters_linked_eta, h_tksters_tot_eta);
  h_eff_phi->Divide(h_tksters_tot_phi);
  h_eff_pT->Divide(h_tksters_tot_pT);
  h_eff_E->Divide(h_tksters_tot_E);

  h_st_eff_eta->Divide(h_stksters_linked_eta, h_stksters_tot_eta);
  h_st_eff_phi->Divide(h_stksters_tot_phi);
  h_st_eff_pT->Divide(h_stksters_tot_pT);
  h_st_eff_E->Divide(h_stksters_tot_E);

  for (int i = 0; i < 100; ++i) {
    layer_Efrac_mean_TS_EM[i] /= NEvent;
    layer_Efrac_mean_TS_HAD[i] /= NEvent;
    layer_Efrac_mean_STS[i] /= NEvent;
    stdDev_by_layer_EM[i] /= NEvent;
    stdDev_by_layer_HAD[i] /= NEvent;
  }

  for (int i = 1; i <= 100; ++i) {
    h_EfracLayer_TS_EM->SetBinContent(i, layer_Efrac_mean_TS_EM[i - 1]);
    h_EfracLayer_TS_HAD->SetBinContent(i, layer_Efrac_mean_TS_HAD[i - 1]);
    h_EfracLayer_STS->SetBinContent(i, layer_Efrac_mean_STS[i - 1]);
    h_stdDevLayer_TS_EM->SetBinContent(i, stdDev_by_layer_EM[i - 1]);
    h_stdDevLayer_TS_HAD->SetBinContent(i, stdDev_by_layer_HAD[i - 1]);
  }

  h_EfracLayer_TS_EM->SetLineColor(kRed);
  h_EfracLayer_TS_EM->SetLineWidth(2);
  hs_EfracLayer->Add(h_EfracLayer_TS_EM);
  h_EfracLayer_TS_HAD->SetLineColor(kBlue);
  h_EfracLayer_TS_HAD->SetLineWidth(2);
  hs_EfracLayer->Add(h_EfracLayer_TS_HAD);
  h_EfracLayer_STS->SetLineColor(kOrange);
  h_EfracLayer_STS->SetLineWidth(2);
  hs_EfracLayer->Add(h_EfracLayer_STS);

  if (hs_stdDevLayer_TS) {
    h_stdDevLayer_TS_EM->SetLineColor(kRed);
    h_stdDevLayer_TS_EM->SetLineWidth(2);
    hs_stdDevLayer_TS->Add(h_stdDevLayer_TS_EM);
    h_stdDevLayer_TS_HAD->SetLineColor(kBlue);
    h_stdDevLayer_TS_HAD->SetLineWidth(2);
    hs_stdDevLayer_TS->Add(h_stdDevLayer_TS_HAD);
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TiclDebugger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersMerge", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("label_simTst", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("layer_clustersTime", edm::InputTag("hgcalLayerClusters", "timeLayerCluster"));
  desc.add<std::string>("detector", "HGCAL");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  desc.add<edm::InputTag>("linkedTracksters", edm::InputTag("ticlTrackstersMerge", "linkedTrackster"));
  desc.add<edm::InputTag>("linkedSimTracksters", edm::InputTag("ticlTrackstersMerge", "linkedSimTrackster"));
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  descriptions.add("ticlDebugger", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(TiclDebugger);
