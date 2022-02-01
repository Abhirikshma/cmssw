#include <cmath>
#include <array>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByPCAGeometric.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

using namespace ticl;

LinkingAlgoByPCAGeometric::LinkingAlgoByPCAGeometric(const edm::ParameterSet &conf) : LinkingAlgoBase(conf) {}

LinkingAlgoByPCAGeometric::~LinkingAlgoByPCAGeometric() {}

void LinkingAlgoByPCAGeometric::initialize(const HGCalDDDConstants *hgcons,
                                           const hgcal::RecHitTools rhtools,
                                           const edm::ESHandle<MagneticField> bfieldH,
                                           const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;

  buildFirstLayers();

  rhtools_ = rhtools;

  //bFieldProd = &es.getData(bfieldToken_);
  bfield_ = bfieldH;
  propagator_ = propH;
}

void LinkingAlgoByPCAGeometric::buildFirstLayers() {
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

void LinkingAlgoByPCAGeometric::linkTracksters(const std::vector<reco::Track> &tracks,
                                               const StringCutObjectSelector<reco::Track> cutTk,
                                               const std::vector<CaloParticle> &caloParticles,
                                               const std::vector<Trackster> &tracksters,
                                               const std::vector<Trackster> &simTracksters,
                                               std::vector<SuperTrackster> &resultTracksters,
                                               std::vector<SuperTrackster> &resultSimTracksters) {
  const double simTracksterEnergyCut = 0.;  // S/TS < this frac of CP energy not selected
  const double tracksterEnergyCut = 0.;
  const double delta_search = 0.02;  // search box delta in eta-phi

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  buildFirstLayers();

  // propagated point collections
  std::vector<std::pair<std::pair<Vector, Vector>, int>>
      trackPcol;  // propagated track points, errors and index of track in collection
  std::vector<std::pair<Vector, int>> tracksterPcol;     // same but for Tracksters, errors not included here yet
  std::vector<std::pair<Vector, int>> simtracksterPcol;  // simTracksters

  // tiles
  TICLLayerTile tracksterProp_tile_fw;
  TICLLayerTile tracksterProp_tile_bw;
  TICLLayerTile simtracksterProp_tile_fw;
  TICLLayerTile simtracksterProp_tile_bw;
  TICLLayerTile trackProp_tile_fw;
  TICLLayerTile trackProp_tile_bw;

  bool singleSC_zplus = false;
  bool singleSC_zminus = false;

  double CP_zlpusE = 0.;  // Energy cuts on S/TS as frac of CP energy
  double CP_zminusE = 0.;

  for (auto const &cp : caloParticles) {
    if (cp.g4Tracks()[0].crossedBoundary()) {
      if (cp.g4Tracks()[0].getPositionAtBoundary().Z() > 0) {
        singleSC_zplus = true;
        CP_zlpusE = cp.energy();
      } else if (cp.g4Tracks()[0].getPositionAtBoundary().Z() < 0) {
        singleSC_zminus = true;
        CP_zminusE = cp.energy();
      }
    }
  }  // CP

  // Propagate tracks to the HGCal front
  for (unsigned i = 0; i < tracks.size(); ++i) {
    const auto tk = tracks[i];
    if (!cutTk((tk))) {
      continue;
    }

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
    int iSide = int(tk.eta() > 0);
    TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());

    if (tsos.isValid()) {
      Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
      Vector tkP_momentum(tsos.globalMomentum().x(), tsos.globalMomentum().y(), tsos.globalMomentum().z());
      double x_err = pow(tsos.localError().positionError().xx(), 0.5);
      double y_err = pow(tsos.localError().positionError().yy(), 0.5);
      Vector track_err(x_err, y_err, 0.0);

      trackPcol.push_back(std::make_pair(std::make_pair(trackP, track_err), i));

      if (trackP.Eta() > 0 && singleSC_zplus)
        trackProp_tile_fw.fill(trackP.Eta(), trackP.Phi(), i);

      else if (trackP.Eta() < 0 && singleSC_zminus)
        trackProp_tile_bw.fill(trackP.Eta(), trackP.Phi(), i);
    }
  }

  // Propagate simTracksters
  std::vector<unsigned> selectedSTS_idx;

  for (unsigned i = 0; i < simTracksters.size(); ++i) {
    const auto &st = simTracksters[i];

    Vector baryc = st.barycenter();
    Vector directnv = st.eigenvectors(0);

    if (abs(directnv.Z()) < 0.00001)
      continue;

    bool zplus_reject = true;
    bool zminus_reject = true;

    if (baryc.Z() > 0 && singleSC_zplus) {
      if (st.raw_energy() > simTracksterEnergyCut * CP_zlpusE)
        zplus_reject = false;
    }

    else if (baryc.Z() < 0 && singleSC_zminus) {
      if (st.raw_energy() > simTracksterEnergyCut * CP_zminusE)
        zminus_reject = false;
    }

    if (!(zplus_reject && zminus_reject)) {
      float zVal = hgcons_->waferZ(1, true);
      zVal *= (baryc.Z() > 0) ? 1 : -1;

      double par = (zVal - baryc.Z()) / directnv.Z();

      double xOnSurface = par * directnv.X() + baryc.X();
      double yOnSurface = par * directnv.Y() + baryc.Y();

      Vector simtracksterP(xOnSurface, yOnSurface, zVal);
      selectedSTS_idx.push_back(i);
      simtracksterPcol.push_back(std::make_pair(simtracksterP, i));

      if (simtracksterP.Eta() > 0)
        simtracksterProp_tile_fw.fill(simtracksterP.Eta(), simtracksterP.Phi(), i);

      else if (simtracksterP.Eta() < 0)
        simtracksterProp_tile_bw.fill(simtracksterP.Eta(), simtracksterP.Phi(), i);
    }

  }  // STS

  // Propagate tracksters
  std::vector<unsigned> selectedTS_idx;
  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto &t = tracksters[i];

    Vector baryc = t.barycenter();
    Vector directnv = t.eigenvectors(0);
    Vector trackster_err(pow(t.sigmas()[0], 0.5), pow(t.sigmas()[1], 0.5), 0.0);

    bool zplus_reject = true;
    bool zminus_reject = true;

    if (baryc.Z() > 0 && singleSC_zplus) {
      if (t.raw_energy() > tracksterEnergyCut * CP_zlpusE)
        zplus_reject = false;
    }

    else if (baryc.Z() < 0 && singleSC_zminus) {
      if (t.raw_energy() > tracksterEnergyCut * CP_zminusE)
        zminus_reject = false;
    }

    if (!(zplus_reject && zminus_reject)) {
      float zVal = hgcons_->waferZ(1, true);
      zVal *= (baryc.Z() > 0) ? 1 : -1;

      double par = (zVal - baryc.Z()) / directnv.Z();

      double xOnSurface = par * directnv.X() + baryc.X();
      double yOnSurface = par * directnv.Y() + baryc.Y();

      Vector tracksterP(xOnSurface, yOnSurface, zVal);
      selectedTS_idx.push_back(i);
      tracksterPcol.push_back(std::make_pair(tracksterP, i));

      if (tracksterP.Eta() > 0)
        tracksterProp_tile_fw.fill(tracksterP.Eta(), tracksterP.Phi(), i);

      else if (tracksterP.Eta() < 0)
        tracksterProp_tile_bw.fill(tracksterP.Eta(), tracksterP.Phi(), i);
    }
  }  // TS

  // Search box over tiles + preliminary linking
  std::vector<unsigned> tracksters_near[trackPcol.size()];  // i-th: indices of tracksters 'linked' to track i
  std::vector<unsigned> simtracksters_near[trackPcol.size()];
  int nth = 0;
  for (auto i : trackPcol) {
    auto trackP = i.first.first;
    //auto track_err = i.first.second;

    double tk_eta = trackP.Eta();
    double tk_phi = trackP.Phi();

    double delta = delta_search;

    if (tk_eta > 0) {
      double eta_min = ((tk_eta - delta) > 0) ? (tk_eta - delta) : 0;

      std::array<int, 4> search_box =
          tracksterProp_tile_fw.searchBoxEtaPhi(eta_min, tk_eta + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto tracksters_in_box = tracksterProp_tile_fw[tracksterProp_tile_fw.globalBin(eta_i, phi_i)];
          tracksters_near[nth].insert(
              std::end(tracksters_near[nth]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // TS

      std::array<int, 4> search_box_st =
          simtracksterProp_tile_fw.searchBoxEtaPhi(eta_min, tk_eta + delta, tk_phi - delta, tk_phi + delta);
      if (search_box_st[2] > search_box_st[3]) {
        double temp = search_box_st[3];
        search_box_st[3] = search_box_st[2];
        search_box_st[2] = temp;
      }

      for (int eta_i = search_box_st[0]; eta_i <= search_box_st[1]; ++eta_i) {
        for (int phi_i = search_box_st[2]; phi_i <= search_box_st[3]; ++phi_i) {
          auto simtracksters_in_box = simtracksterProp_tile_fw[simtracksterProp_tile_fw.globalBin(eta_i, phi_i)];
          simtracksters_near[nth].insert(
              std::end(simtracksters_near[nth]), std::begin(simtracksters_in_box), std::end(simtracksters_in_box));
        }
      }  // STS
    }    // forward

    if (tk_eta < 0) {
      double eta_min = ((abs(tk_eta) - delta) > 0) ? (abs(tk_eta) - delta) : 0;

      std::array<int, 4> search_box =
          tracksterProp_tile_bw.searchBoxEtaPhi(eta_min, abs(tk_eta) + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto tracksters_in_box = tracksterProp_tile_bw[tracksterProp_tile_bw.globalBin(eta_i, phi_i)];
          tracksters_near[nth].insert(
              std::end(tracksters_near[nth]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // TS

      std::array<int, 4> search_box_st =
          simtracksterProp_tile_bw.searchBoxEtaPhi(eta_min, abs(tk_eta) + delta, tk_phi - delta, tk_phi + delta);
      if (search_box_st[2] > search_box_st[3]) {
        double temp = search_box_st[3];
        search_box_st[3] = search_box_st[2];
        search_box_st[2] = temp;
      }

      for (int eta_i = search_box_st[0]; eta_i <= search_box_st[1]; ++eta_i) {
        for (int phi_i = search_box_st[2]; phi_i <= search_box_st[3]; ++phi_i) {
          auto simtracksters_in_box = simtracksterProp_tile_bw[simtracksterProp_tile_bw.globalBin(eta_i, phi_i)];
          simtracksters_near[nth].insert(
              std::end(simtracksters_near[nth]), std::begin(simtracksters_in_box), std::end(simtracksters_in_box));
        }
      }  // STS
    }    // backward

    ++nth;
  }  // propagated tracks

  for (unsigned i = 0; i < trackPcol.size(); ++i) {
    auto const trackstersLinked = tracksters_near[i];
    const unsigned track = i;
    SuperTrackster resultTS(trackstersLinked, track);
    resultTracksters.push_back(resultTS);

    auto const simTrackstersLinked = simtracksters_near[i];
    SuperTrackster resultSTS(simTrackstersLinked, track);
    resultSimTracksters.push_back(resultSTS);
  }

}  // linkTracksters

void LinkingAlgoByPCAGeometric::fillPSetDescription(edm::ParameterSetDescription &desc) {
  LinkingAlgoBase::fillPSetDescription(desc);
}
