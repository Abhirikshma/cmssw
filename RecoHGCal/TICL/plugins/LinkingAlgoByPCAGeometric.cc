#include <cmath>
#include <ostream>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByPCAGeometric.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

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

math::XYZVector LinkingAlgoByPCAGeometric::propagateTrackster(
    const Trackster &t, const unsigned idx, float zVal, std::array<TICLLayerTile, 2>& tracksterTiles) {
  // any energy or caloparticle based selection has to be handled outside
  Vector baryc = t.barycenter();
  Vector directnv = t.eigenvectors(0);

  assert(abs(directnv.Z()) > 0.00001);

  zVal *= (baryc.Z() > 0) ? 1 : -1;

  double par = (zVal - baryc.Z()) / directnv.Z();
  double xOnSurface = par * directnv.X() + baryc.X();
  double yOnSurface = par * directnv.Y() + baryc.Y();
  Vector tPoint(xOnSurface, yOnSurface, zVal);

  if (tPoint.Eta() > 0)
    tracksterTiles[1].fill(tPoint.Eta(), tPoint.Phi(), idx);

  else if (tPoint.Eta() < 0)
    tracksterTiles[0].fill(tPoint.Eta(), tPoint.Phi(), idx);

  return tPoint;
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
                                               const std::vector<Trackster> &tracksters,
                                               std::vector<SuperTrackster> &resultLinked) {
  // Selections based on CaloParticles or energy have to be implemented outside

  const double delta = 0.02;  // search box delta in eta-phi

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  buildFirstLayers();

  // propagated point collections
  std::vector<std::pair<std::pair<Vector, Vector>, int>>
      trackPcol;  // propagated track points, errors and index of track in collection
  std::vector<std::pair<Vector, int>> tracksterPcol;     // same but for Tracksters, errors not included here yet

  // tiles
  std::array<TICLLayerTile, 2> tracksterPropTiles = {}; // layer 0 is bw, 1 is fw
  std::array<TICLLayerTile, 2> trackPropTiles = {};


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

      if (trackP.Eta() > 0)
        trackPropTiles[1].fill(trackP.Eta(), trackP.Phi(), i);

      else if (trackP.Eta() < 0)
        trackPropTiles[0].fill(trackP.Eta(), trackP.Phi(), i);
    }
  }  // Tracks

  // Propagate tracksters
  std::vector<unsigned> selectedTS_idx;
  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto &t = tracksters[i];
    Vector directnv = t.eigenvectors(0);

    if (abs(directnv.Z()) < 0.00001)
      continue;

    float zVal = hgcons_->waferZ(1, true);
    Vector tsP = propagateTrackster(t, i, zVal, tracksterPropTiles);
    tracksterPcol.push_back(std::make_pair(tsP, i));
    selectedTS_idx.push_back(i);
  }  // TS

  // Search box over trackster tiles + preliminary linking
  std::vector<unsigned> tracksters_near[trackPcol.size()];  // i-th element: vector of indices of tracksters 'linked' to track i
  int nth = 0;
  for (auto i : trackPcol) {
    auto trackP = i.first.first;
    //auto track_err = i.first.second;

    double tk_eta = trackP.Eta();
    double tk_phi = trackP.Phi();

    if (tk_eta > 0) {
      double eta_min = std::max(tk_eta - delta, 0.);

      const TICLLayerTile &tile = tracksterPropTiles[1];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, tk_eta + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          tracksters_near[nth].insert(
              std::end(tracksters_near[nth]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // TS
    }    // forward

    if (tk_eta < 0) {
      double eta_min = std::max(abs(tk_eta) - delta, 0.);

      const TICLLayerTile &tile = tracksterPropTiles[0];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, abs(tk_eta) + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          tracksters_near[nth].insert(
              std::end(tracksters_near[nth]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // TS
    }    // backward

    ++nth;
  }  // propagated tracks

  for (unsigned i = 0; i < trackPcol.size(); ++i) {
    auto const trackstersLinked = tracksters_near[i];
    const unsigned track = i;
    SuperTrackster resultTS(trackstersLinked, track);
    resultLinked.push_back(resultTS);
  }

}  // linkTracksters

void LinkingAlgoByPCAGeometric::fillPSetDescription(edm::ParameterSetDescription &desc) {
  LinkingAlgoBase::fillPSetDescription(desc);
}
