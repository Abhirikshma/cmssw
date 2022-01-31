#ifndef RECOHGCAL_TICL_TRACKSTERSPCA_MOD_H
#define RECOHGCAL_TICL_TRACKSTERSPCA_MOD_H

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <iostream>
#include <array>
#include <vector>
#include <set>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace ticl;

typedef math::XYZVector Vector;

unsigned getLayerFromLC(const reco::CaloCluster LC, hgcal::RecHitTools rhtools) {
  std::vector<std::pair<DetId, float>> thisclusterHits = LC.hitsAndFractions();
  auto layer = rhtools.getLayerWithOffset(thisclusterHits[0].first);
  return layer;
}

Vector barycenterInLayer(const Trackster &ts,
                         const std::vector<reco::CaloCluster> &layerClusters,
                         std::vector<unsigned> vertices_in_layer) {
  double wt_x = 0.;
  double wt_y = 0.;
  double wt_z = 0.;
  double wt = 0.;
  for (auto v : vertices_in_layer) { // in a layer
    auto thisLC = layerClusters[ts.vertices(v)];
    double thisE = thisLC.energy();
    wt_x += thisE*thisLC.x();
    wt_y += thisE*thisLC.y();
    wt_z += thisE*thisLC.z();
    wt += thisE;
  }
  Vector result(wt_x/wt, wt_y/wt, wt_z/wt);
  return result;
}

std::vector<std::vector<unsigned>> sortByLayer(const Trackster &ts,
                                               const std::vector<reco::CaloCluster> &layerClusters,
                                               hgcal::RecHitTools rhtools) {
  size_t N = ts.vertices().size();

  std::vector<std::vector<unsigned>> result;
  result.resize(rhtools.lastLayer()+1);

  for (unsigned i = 0; i < N; ++i) {
    auto thisLC = layerClusters[ts.vertices(i)];
    auto layer = getLayerFromLC(thisLC, rhtools);
    result[layer].push_back(i);
  }
  return result;
}

std::array<std::vector<unsigned>, 100> translateTrackster(const Trackster &ts,
                                        const std::vector<reco::CaloCluster> &layerClusters,
                                        hgcal::RecHitTools rhtools,
                                        unsigned maxE_layer,
                                        int shift_to = 50) {
  std::array<std::vector<unsigned>, 100> result; 
  int shift = shift_to - (int)maxE_layer;

  auto vert_by_layer = sortByLayer(ts, layerClusters, rhtools);

  for (unsigned i = 1; i <= rhtools.lastLayer(); ++i) {
    auto vert_in_layer = vert_by_layer[i];
    if (vert_in_layer.empty()) continue;
    for (auto v : vert_in_layer) result[i+shift-1].push_back(v);
  }
  return result;
}

std::array<double, 100> assignPCAtoTracksters_mod(std::vector<Trackster> &tracksters,
                                 const std::vector<reco::CaloCluster> &layerClusters,
                                 const edm::ValueMap<std::pair<float, float>> &layerClustersTime,
                                 hgcal::RecHitTools rhtools,
                                 double z_limit_em,
                                 bool energyWeight=true) {
  LogDebug("TrackstersPCA_Eigen") << "------- Eigen -------" << std::endl;
  int nonZero_trackster = 0;
  std::array<double, 100> layer_Efrac = {};
  for (auto &trackster : tracksters) {
    Eigen::Vector3d point;
    point << 0., 0., 0.;
    Eigen::Vector3d barycenter;
    barycenter << 0., 0., 0.;
    Eigen::Vector3d filtered_barycenter;
    filtered_barycenter << 0., 0., 0.;

    auto fillPoint = [&](const reco::CaloCluster &c, const float weight = 1.f) {
      point[0] = weight * c.x();
      point[1] = weight * c.y();
      point[2] = weight * c.z();
    };

    // Initialize this trackster with default, dummy values
    trackster.setRawEnergy(0.f);
    trackster.setRawEmEnergy(0.f);
    trackster.setRawPt(0.f);
    trackster.setRawEmPt(0.f);

    size_t N = trackster.vertices().size();
    if (N == 0) continue;
    float weight = 1.f / N;
    float weights2_sum = 0.f;
    Eigen::Vector3d sigmas;
    sigmas << 0., 0., 0.;
    Eigen::Vector3d sigmasEigen;
    sigmasEigen << 0., 0., 0.;
    Eigen::Matrix3d covM = Eigen::Matrix3d::Zero();

    std::vector<float> times;
    std::vector<float> timeErrors;
    std::set<uint32_t> usedLC;


    for (size_t i = 0; i < N; ++i) {
      auto fraction = 1.f / trackster.vertex_multiplicity(i);
      trackster.addToRawEnergy(layerClusters[trackster.vertices(i)].energy() * fraction);
      if (std::abs(layerClusters[trackster.vertices(i)].z()) <= z_limit_em)
        trackster.addToRawEmEnergy(layerClusters[trackster.vertices(i)].energy() * fraction);

      // Compute the weighted barycenter.
      if (energyWeight)
        weight = layerClusters[trackster.vertices(i)].energy() * fraction;
      fillPoint(layerClusters[trackster.vertices(i)], weight);
      for (size_t j = 0; j < 3; ++j)
        barycenter[j] += point[j];

      // Add timing from layerClusters not already used
      if ((usedLC.insert(trackster.vertices(i))).second) {
        float timeE = layerClustersTime.get(trackster.vertices(i)).second;
        if (timeE > 0.f) {
          times.push_back(layerClustersTime.get(trackster.vertices(i)).first);
          timeErrors.push_back(1. / pow(timeE, 2));
        }
      }
    }
    if (energyWeight && trackster.raw_energy())
      barycenter /= trackster.raw_energy();

    hgcalsimclustertime::ComputeClusterTime timeEstimator;
    std::pair<float, float> timeTrackster = timeEstimator.fixSizeHighestDensity(times, timeErrors);

    // LayerCluster filtering for the modified PCA
    auto vertices_by_layer = sortByLayer(trackster, layerClusters, rhtools);
    std::vector<double> layerClusterenergies; 

    for (unsigned i = 0; i <= rhtools.lastLayer(); ++i) {
      auto vertices_in_layer = vertices_by_layer[i];
      for (auto v : vertices_in_layer) {
        layerClusterenergies.push_back(layerClusters[trackster.vertices(v)].energy());
      }
    }

    // max E LC
    auto result_glob = std::max_element(layerClusterenergies.begin(), layerClusterenergies.end());
    auto maxE_vertex = std::distance(layerClusterenergies.begin(), result_glob);
    auto maxE_layer = getLayerFromLC(layerClusters[trackster.vertices(maxE_vertex)], rhtools);
    //double maxEnergy = layerClusters[trackster.vertices(maxE_vertex)].energy();
    double maxElayer_energy = 0; // energy in layer of max E LC
    for (auto vrt : vertices_by_layer[maxE_layer]) maxElayer_energy += layerClusters[trackster.vertices(vrt)].energy();

    std::vector<unsigned> filtered_idx; // higest energy vertices in a layer
    double filtered_energy = 0;

    for (unsigned i = 1; i <= rhtools.lastLayer(); ++i) {
      auto vertices_in_layer = vertices_by_layer[i];
      if (vertices_in_layer.empty()) continue;

      std::vector<double> energies_in_layer;
      for (auto vrt : vertices_in_layer) energies_in_layer.push_back(layerClusters[trackster.vertices(vrt)].energy());
    
      auto result_inlayer = std::max_element(energies_in_layer.begin(),energies_in_layer.end());
      unsigned maxE_id = std::distance(energies_in_layer.begin(),result_inlayer);
      
      if ((int)i >= (int)maxE_layer - 20 && (int)i <= (int)maxE_layer + 5) {
        auto filtered_vert = vertices_in_layer[maxE_id];
        filtered_idx.push_back(filtered_vert);
        
        const auto maxE_LC = layerClusters[trackster.vertices(filtered_vert)]; 
        fillPoint(maxE_LC, maxE_LC.energy()*(1.f/trackster.vertex_multiplicity(filtered_vert)));
        for (size_t j = 0; j < 3; ++j) filtered_barycenter[j] += point[j];
        filtered_energy += maxE_LC.energy();
      }
      
    }

    filtered_barycenter /= filtered_energy;

    // Translate trackster to have layer of max E LC at 50
    auto vertT_by_layer = translateTrackster(trackster, layerClusters, rhtools, maxE_layer);
    
    // E frac by layer for translated trackster
    for (unsigned i = 0; i < 100; ++i) {
      auto vertices_in_layer = vertT_by_layer[i];
      //if (vertices_in_layer.empty()) continue;
      double e_layer = 0.;
      for (auto v : vertices_in_layer) e_layer += layerClusters[trackster.vertices(v)].energy();
      //e_layer /= maxEnergy;
      e_layer /= maxElayer_energy;
      layer_Efrac[i] += e_layer;
    }

    std::cout << "trackster energy = " << trackster.raw_energy() << std::endl;

    // Compute the Covariance Matrix and the sum of the squared weights, used
    // to compute the correct normalization.
    // The barycenter has to be known.
    for (size_t i : filtered_idx) {
      fillPoint(layerClusters[trackster.vertices(i)]);
      if (energyWeight && trackster.raw_energy())
        weight = (layerClusters[trackster.vertices(i)].energy() / trackster.vertex_multiplicity(i)) / trackster.raw_energy();
        //weight = (layerClusters[trackster.vertices(i)].energy() / trackster.vertex_multiplicity(i)) / filtered_energy;
      weights2_sum += weight * weight;
      for (size_t x = 0; x < 3; ++x)
        for (size_t y = 0; y <= x; ++y) {
          covM(x, y) += weight * (point[x] - barycenter[x]) * (point[y] - barycenter[y]);
          covM(y, x) = covM(x, y);
        }
    }
 
    covM *= 1. / (1. - weights2_sum);

    // Perform the actual decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>::RealVectorType eigenvalues_fromEigen;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>::EigenvectorsType eigenvectors_fromEigen;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covM);
    if (eigensolver.info() != Eigen::Success) {
      eigenvalues_fromEigen = eigenvalues_fromEigen.Zero();
      eigenvectors_fromEigen = eigenvectors_fromEigen.Zero();
    } else {
      eigenvalues_fromEigen = eigensolver.eigenvalues();
      eigenvectors_fromEigen = eigensolver.eigenvectors();
    }

    // Compute the spread in the both spaces.
    for (size_t i = 0; i < N; ++i) {
      fillPoint(layerClusters[trackster.vertices(i)]);
      sigmas += weight * (point - barycenter).cwiseAbs2();
      Eigen::Vector3d point_transformed = eigenvectors_fromEigen * (point - barycenter);
      if (energyWeight && trackster.raw_energy())
        weight =
            (layerClusters[trackster.vertices(i)].energy() / trackster.vertex_multiplicity(i)) / trackster.raw_energy();
      sigmasEigen += weight * (point_transformed.cwiseAbs2());
    }
    sigmas /= (1. - weights2_sum);
    sigmasEigen /= (1. - weights2_sum);

    // Add trackster attributes
    trackster.setBarycenter(ticl::Trackster::Vector(barycenter));
    trackster.setTimeAndError(timeTrackster.first, timeTrackster.second);
    trackster.fillPCAVariables(
        eigenvalues_fromEigen, eigenvectors_fromEigen, sigmas, sigmasEigen, 3, ticl::Trackster::PCAOrdering::ascending);

    LogDebug("TrackstersPCA") << "Use energy weighting: " << energyWeight << std::endl;
    LogDebug("TrackstersPCA") << "\nTrackster characteristics: " << std::endl;
    LogDebug("TrackstersPCA") << "Size: " << N << std::endl;
    LogDebug("TrackstersPCA") << "Energy: " << trackster.raw_energy() << std::endl;
    LogDebug("TrackstersPCA") << "raw_pt: " << trackster.raw_pt() << std::endl;
    LogDebug("TrackstersPCA") << "Means:          " << barycenter[0] << ", " << barycenter[1] << ", " << barycenter[2]
                              << std::endl;
    LogDebug("TrackstersPCA") << "Time:          " << trackster.time() << " +/- " << trackster.timeError() << std::endl;
    LogDebug("TrackstersPCA") << "EigenValues from Eigen/Tr(cov): " << eigenvalues_fromEigen[2] / covM.trace() << ", "
                              << eigenvalues_fromEigen[1] / covM.trace() << ", "
                              << eigenvalues_fromEigen[0] / covM.trace() << std::endl;
    LogDebug("TrackstersPCA") << "EigenValues from Eigen:         " << eigenvalues_fromEigen[2] << ", "
                              << eigenvalues_fromEigen[1] << ", " << eigenvalues_fromEigen[0] << std::endl;
    LogDebug("TrackstersPCA") << "EigenVector 3 from Eigen: " << eigenvectors_fromEigen(0, 2) << ", "
                              << eigenvectors_fromEigen(1, 2) << ", " << eigenvectors_fromEigen(2, 2) << std::endl;
    LogDebug("TrackstersPCA") << "EigenVector 2 from Eigen: " << eigenvectors_fromEigen(0, 1) << ", "
                              << eigenvectors_fromEigen(1, 1) << ", " << eigenvectors_fromEigen(2, 1) << std::endl;
    LogDebug("TrackstersPCA") << "EigenVector 1 from Eigen: " << eigenvectors_fromEigen(0, 0) << ", "
                              << eigenvectors_fromEigen(1, 0) << ", " << eigenvectors_fromEigen(2, 0) << std::endl;
    LogDebug("TrackstersPCA") << "Original sigmas:          " << sigmas[0] << ", " << sigmas[1] << ", " << sigmas[2]
                              << std::endl;
    LogDebug("TrackstersPCA") << "SigmasEigen in PCA space: " << sigmasEigen[2] << ", " << sigmasEigen[1] << ", "
                              << sigmasEigen[0] << std::endl;
    LogDebug("TrackstersPCA") << "covM:     \n" << covM << std::endl;
    ++nonZero_trackster;
  }

  // mean E frac per layer
  if (nonZero_trackster > 0)
  for (int i = 0; i < 100; ++i) {
    layer_Efrac[i] /= nonZero_trackster;
  }
  return layer_Efrac;
}
#endif