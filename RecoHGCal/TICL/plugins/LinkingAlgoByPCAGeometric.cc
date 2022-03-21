#include <cmath>
#include <ostream>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByPCAGeometric.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

using namespace ticl;

LinkingAlgoByPCAGeometric::LinkingAlgoByPCAGeometric(const edm::ParameterSet &conf) : LinkingAlgoBase(conf) {}

LinkingAlgoByPCAGeometric::~LinkingAlgoByPCAGeometric() {}

void LinkingAlgoByPCAGeometric::initialize(const HGCalDDDConstants *hgcons,
                                           const hgcal::RecHitTools rhtools,
                                           const edm::ESHandle<MagneticField> bfieldH,
                                           const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  buildLayers();

  //bFieldProd = &es.getData(bfieldToken_);
  bfield_ = bfieldH;
  propagator_ = propH;
}

math::XYZVector LinkingAlgoByPCAGeometric::propagateTrackster(
    const Trackster &t, const unsigned idx, float zVal, std::array<TICLLayerTile, 2>& tracksterTiles) {
  // any energy or caloparticle based selection has to be handled outside
  // need to only provide the positive Z co-ordinate of the surface to propagate to
  // the correct sign is calculated inside according to the barycenter of trackster 
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

void LinkingAlgoByPCAGeometric::buildLayers() {
  float zVal = hgcons_->waferZ(1, true);
  std::pair<double, double> rMinMax = hgcons_->rangeR(zVal, true);

  float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
  std::pair<double, double> rMinMax_interface = hgcons_->rangeR(zVal_interface, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());

    zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
    interfaceDisk_[iSide] = 
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second, zSide - 0.5, zSide + 0.5))
                                      .get());
  }
}

void LinkingAlgoByPCAGeometric::linkTracksters(const std::vector<reco::Track> &tracks,
                                               const edm::Handle<std::vector<reco::Track>> tkH,
                                               const StringCutObjectSelector<reco::Track> cutTk,
                                               const edm::OrphanHandle<std::vector<Trackster>> tsH,
                                               std::vector<TICLCandidate> &resultLinked) {
  // Selections based on CaloParticles or energy have to be implemented outside

  // search box deltas in eta-phi
  const double delta = 0.02;  // track -> trackster
  const double del_ts = 0.02;  // trackster CE-E -> CE-H
  const double del_tsHad = 0.02; // CE-H -> CE-H 

  const auto tracksters = *tsH;
  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  // propagated point collections
  // elements in the propagated points collecions are used
  // to look for potential linkages in the appropriate tiles
  std::vector<std::pair<std::pair<Vector, Vector>, unsigned>>
      trackPColl;  // propagated track points, errors and index of track in collection
  std::vector<std::pair<Vector, unsigned>> tkPropIntColl; // tracks propagated to lastLayerEE
  std::vector<std::pair<Vector, unsigned>> tsPropIntColl; // Tracksters in CE-E, propagated to lastLayerEE
  std::vector<std::pair<Vector, unsigned>> tsHadPropIntColl; // Tracksters in CE-H, propagated to lastLayerEE

  // tiles, layer 0 is bw, 1 is fw
  std::array<TICLLayerTile, 2> tracksterPropTiles = {}; // all Tracksters, propagated to layer 1
  std::array<TICLLayerTile, 2> tsPropIntTiles = {}; // all Tracksters, propagated to lastLayerEE
  std::array<TICLLayerTile, 2> tsHadPropIntTiles = {}; // Tracksters in CE-H, propagated to lastLayerEE
  std::array<TICLLayerTile, 2> dummyTiles = {}; // dummy tile to pass to the propagate function


  // Propagate tracks
  for (unsigned i = 0; i < tracks.size(); ++i) {
    const auto tk = tracks[i];
    if (!cutTk((tk))) {
      continue;
    }

    FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
    int iSide = int(tk.eta() > 0);

    // to the HGCal front
    TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    if (tsos.isValid()) {
      Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
      Vector tkP_momentum(tsos.globalMomentum().x(), tsos.globalMomentum().y(), tsos.globalMomentum().z());
      double x_err = pow(tsos.localError().positionError().xx(), 0.5);
      double y_err = pow(tsos.localError().positionError().yy(), 0.5);
      Vector track_err(x_err, y_err, 0.0);

      trackPColl.push_back(std::make_pair(std::make_pair(trackP, track_err), i));
    }
    // to lastLayerEE
    tsos = prop.propagate(fts, interfaceDisk_[iSide]->surface());
    if (tsos.isValid()) {
      Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());

      tkPropIntColl.push_back(std::make_pair(trackP, i));
    }
  }  // Tracks

  // Propagate tracksters
  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto &t = tracksters[i];
    Vector directnv = t.eigenvectors(0);

    if (abs(directnv.Z()) < 0.00001)
      continue;

    // to HGCal front
    float zVal = hgcons_->waferZ(1, true);
    Vector tsP = propagateTrackster(t, i, zVal, tracksterPropTiles);

    // to lastLayerEE
    zVal = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
    tsP = propagateTrackster(t, i, zVal, tsPropIntTiles);
  }  // TS

  // step 3: linking tracks -> all tracksters, at layer 1
  std::vector<unsigned> tracksters_near[tracks.size()] = {}; // i-th element: vector of indices of tracksters 'linked' to track i

  for (auto i : trackPColl) {
    auto trackP = i.first.first;
    const unsigned tkId = i.second;
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
          tracksters_near[tkId].insert(
              std::end(tracksters_near[tkId]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
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
          tracksters_near[tkId].insert(
              std::end(tracksters_near[tkId]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // TS
    }    // backward

  }  // propagated tracks

  // step 4: linking tracks -> all tracksters, at lastLayerEE
  std::vector<unsigned> tsNearTkAtInt[tracks.size()] = {};

  for (auto i : tkPropIntColl) {
    auto trackP = i.first;
    const auto tkId = i.second;

    double tk_eta = trackP.Eta();
    double tk_phi = trackP.Phi();

    if (tk_eta > 0) {
      double eta_min = std::max(tk_eta - delta, 0.);

      const TICLLayerTile &tile = tsPropIntTiles[1];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, tk_eta + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          tsNearTkAtInt[tkId].insert(
              std::end(tsNearTkAtInt[tkId]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // TS
    }    // forward

    if (tk_eta < 0) {
      double eta_min = std::max(abs(tk_eta) - delta, 0.);

      const TICLLayerTile &tile = tsPropIntTiles[0];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, abs(tk_eta) + delta, tk_phi - delta, tk_phi + delta);
      if (search_box[2] > search_box[3]) {
        double temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }

      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          tsNearTkAtInt[tkId].insert(
              std::end(tsNearTkAtInt[tkId]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      }  // TS
    }    // backward
  }

  // Trackster - trackster linking
  // tracksters after the interface are propagated back and filled in a tile
  // and those before are propagated to the interface and used
  // to search for potential linkages in the tile of propagated tracksters
  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto &t = tracksters[i];
    Vector directnv = t.eigenvectors(0);
    if (abs(directnv.Z()) < 0.00001)
      continue;

    float zVal = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z(); // interface of CE-E & CE-H

    Vector tsP = propagateTrackster(t, i, zVal, dummyTiles);
    
    if (abs(t.barycenter().Z()) < abs(zVal)) // trackster in CE-E
      tsPropIntColl.push_back(std::make_pair(tsP, i));
    else { // in CE-H
      tsHadPropIntTiles[(t.barycenter().Z() > 0) ? 1:0].fill(tsP.Eta(), tsP.Phi(), i);
      tsHadPropIntColl.push_back(std::make_pair(tsP, i));
    }
  }
  // step 2: linking Tracksters CE-E -> CE-H
  std::vector<unsigned> tsNearAtInt[tracksters.size()] = {};
  int tsMask2[tracksters.size()] = {0};

  for (auto i : tsPropIntColl) {
    auto ts_eta = i.first.Eta();
    auto ts_phi = i.first.Phi();
    const unsigned tsId = i.second;

    if (ts_eta > 0) {
      double eta_min = std::max(ts_eta - del_ts, 0.);
      const TICLLayerTile &tile = tsHadPropIntTiles[1];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, ts_eta + del_ts, ts_phi - del_ts, ts_phi + del_ts);
      if (search_box[2] > search_box[3]) {
        int temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }
      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          const auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          for (const unsigned t_i : tracksters_in_box) {
            if (!tsMask2[t_i]) {
              tsNearAtInt[tsId].push_back(t_i);
              tsMask2[t_i] = 1;
            }
          }
          //tsNearAtInt[tsId].insert(std::end(tsNearAtInt[tsId]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      } // TS
    } // forward
    if (ts_eta < 0) {
      double eta_min = std::max(abs(ts_eta) - del_ts, 0.);
      const TICLLayerTile &tile = tsHadPropIntTiles[0];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, abs(ts_eta) + del_ts, ts_phi - del_ts, ts_phi + del_ts);
      if (search_box[2] > search_box[3]) {
        int temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }
      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          const auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          for (const unsigned t_i : tracksters_in_box) {
            if (!tsMask2[t_i]) {
              tsNearAtInt[tsId].push_back(t_i);
              tsMask2[t_i] = 1;
            }
          }
          //tsNearAtInt[tsId].insert(std::end(tsNearAtInt[tsId]), std::begin(tracksters_in_box), std::end(tracksters_in_box));
        }
      } // TS
    } // backward

  } // tsPropIntColl

  // step 1: Linking Tracksters CE-H -> CE-H
  std::vector<unsigned> tsHadNearAtInt[tracksters.size()] = {};
  int tsMask1[tracksters.size()] = {0};

  for (auto i : tsHadPropIntColl) {
    double ts_eta = i.first.Eta();
    double ts_phi = i.first.Phi();
    const unsigned tsId = i.second;

    if (ts_eta > 0) {
      double eta_min = std::max(ts_eta - del_tsHad, 0.);
      const TICLLayerTile &tile = tsHadPropIntTiles[1];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, ts_eta + del_tsHad, ts_phi - del_tsHad, ts_phi + del_tsHad);
      if (search_box[2] > search_box[3]) {
        int temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }
      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          const auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          for (const unsigned t_i : tracksters_in_box) {
            if (!tsMask1[t_i]) {
              tsHadNearAtInt[tsId].push_back(t_i);
              tsMask1[t_i] = 1;
            }
          }
        }
      } // TS
    } // forward
    if (ts_eta < 0) {
      double eta_min = std::max(abs(ts_eta) - del_tsHad, 0.);
      const TICLLayerTile &tile = tsHadPropIntTiles[0];
      std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, abs(ts_eta) + del_tsHad, ts_phi - del_tsHad, ts_phi + del_tsHad);
      if (search_box[2] > search_box[3]) {
        int temp = search_box[3];
        search_box[3] = search_box[2];
        search_box[2] = temp;
      }
      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          const auto &tracksters_in_box = tile[tile.globalBin(eta_i, phi_i)];
          for (const unsigned t_i : tracksters_in_box) {
            if (!tsMask1[t_i]) {
              tsHadNearAtInt[tsId].push_back(t_i);
              tsMask1[t_i] = 1;
            }
          }
        }
      } // TS
    } // backward
  } // tsHadPropIntColl

  // make final collections

  int chargedMask[tracksters.size()] = {0};
  for (unsigned i = 0; i < tracks.size(); ++i) {
    if (tracksters_near[i].empty() && tsNearTkAtInt[i].empty()) continue; // nothing linked to track

    TICLCandidate chargedCandidate;
    for (const unsigned ts3_idx : tracksters_near[i]) { // tk -> ts
      if (!chargedMask[ts3_idx]) {
        chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts3_idx));
        chargedMask[ts3_idx] = 1;
      }
      for (const unsigned ts2_idx : tsNearAtInt[ts3_idx]) { // ts_EM -> ts_HAD
        if (!chargedMask[ts2_idx]) {
          chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts2_idx));
          chargedMask[ts2_idx] = 1;
        }
        for (const unsigned ts1_idx : tsHadNearAtInt[ts2_idx]) { // ts_HAD -> ts_HAD
          if (!chargedMask[ts1_idx]) {
            chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts1_idx));
            chargedMask[ts1_idx] = 1;
          }
        }
      }
      for (const unsigned ts1_idx : tsHadNearAtInt[ts3_idx]) { // ts_HAD -> ts_HAD
        if (!chargedMask[ts1_idx]) {
          chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts1_idx));
          chargedMask[ts1_idx] = 1;
        }
      }
    }

    for (const unsigned ts4_idx : tsNearTkAtInt[i]) { // do the same for ts -> ts links at the interface
      if (!chargedMask[ts4_idx]) {
        chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts4_idx));
        chargedMask[ts4_idx] = 1;
      }
      for (const unsigned ts2_idx : tsNearAtInt[ts4_idx]) {
        if (!chargedMask[ts2_idx]) {
          chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts2_idx));
          chargedMask[ts2_idx] = 1;
        }
        for (const unsigned ts1_idx : tsHadNearAtInt[ts2_idx]) {
          if (!chargedMask[ts1_idx]) {
            chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts1_idx));
            chargedMask[ts1_idx] = 1;
          }
        }
      }
      for (const unsigned ts1_idx : tsHadNearAtInt[ts4_idx]) {
        if (!chargedMask[ts1_idx]) {
          chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts1_idx));
          chargedMask[ts1_idx] = 1;
        }
      }
    }

    // do not create a candidate if no tracksters were added to candidate
    // can happen if all the tracksters linked to that track were already masked
    if (chargedCandidate.tracksters().size() > 0) {
      chargedCandidate.setTrackPtr(edm::Ptr<reco::Track>(tkH, i));
      resultLinked.push_back(chargedCandidate);
    }
  }

  std::vector<TICLCandidate> neutralCandidates;
  int neutralMask[tracksters.size()] = {0};
  for (unsigned i = 0; i < tracksters.size(); ++i) {
    if (chargedMask[i]) continue;

    TICLCandidate neutralCandidate;
    if (tsNearAtInt[i].empty() && tsHadNearAtInt[i].empty()) { // nothing linked to this ts
      if (!neutralMask[i]) {
        TICLCandidate neutralNoLinks;
        neutralNoLinks.addTrackster(edm::Ptr<Trackster>(tsH, i));
        neutralMask[i] = 1;
        neutralCandidates.push_back(neutralNoLinks);
      }
    }
    else { // at least one other trackster linked to this
    if (!neutralMask[i]) {
      neutralCandidate.addTrackster(edm::Ptr<Trackster>(tsH, i));
      neutralMask[i] = 1;
    }
    for (const unsigned ts2_idx : tsNearAtInt[i]) {
      if (chargedMask[ts2_idx]) continue;
      if (!neutralMask[ts2_idx]) {
        neutralCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts2_idx));
        neutralMask[ts2_idx] = 1;
      }
      for (const unsigned ts1_idx : tsHadNearAtInt[ts2_idx]) {
        if (chargedMask[ts1_idx]) continue;
        if (!neutralMask[ts1_idx]) {
          neutralCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts1_idx));
          neutralMask[ts1_idx] = 1;
        }
      }
    }
    for (const unsigned ts1_idx : tsHadNearAtInt[i]) {
      if (chargedMask[ts1_idx]) continue;
      if (!neutralMask[ts1_idx]) {
        neutralCandidate.addTrackster(edm::Ptr<Trackster>(tsH, ts1_idx));
        neutralMask[ts1_idx] = 1;
      }
    }
    neutralCandidates.push_back(neutralCandidate);
    }
  }
  
resultLinked.insert(std::end(resultLinked), std::begin(neutralCandidates), std::end(neutralCandidates));

}  // linkTracksters

void LinkingAlgoByPCAGeometric::fillPSetDescription(edm::ParameterSetDescription &desc) {
  LinkingAlgoBase::fillPSetDescription(desc);
}
