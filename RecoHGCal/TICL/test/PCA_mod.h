#ifndef RECOHGCAL_TICL_TRACKSTERSPCA_MOD_H
#define RECOHGCAL_TICL_TRACKSTERSPCA_MOD_H

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <iostream>
#include <vector>
#include <set>

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace ticl;

void assignPCAtoTracksters_mod(std::vector<Trackster> &tracksters,
                                 const std::vector<reco::CaloCluster> &layerClusters,
                                 const edm::ValueMap<std::pair<float, float>> &layerClustersTime,
                                 hgcal::RecHitTools rhtools,
                                 double z_limit_em,
                                 bool energyWeight=true) {
  LogDebug("TrackstersPCA_Eigen") << "------- Eigen -------" << std::endl;

  for (auto &trackster : tracksters) {
    Eigen::Vector3d point;
    point << 0., 0., 0.;
    Eigen::Vector3d barycenter;
    barycenter << 0., 0., 0.;

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

    std::vector<reco::CaloCluster> layerClusters_by_layer[rhtools.lastLayer()+1]; //layers can be from 0 -> lastLayer()
    std::vector<unsigned> vertices_by_layer[rhtools.lastLayer()+1]; 
    std::vector<double> layerClusterenergies; 

    std::cout << std::endl;

    for (size_t i = 0; i < N; ++i) {
      std::cout << "LC " << i << "  layer : ";
      const reco::CaloCluster thislayerCluster = layerClusters[trackster.vertices(i)];
      std::vector<std::pair<DetId, float>> thisclusterHits = thislayerCluster.hitsAndFractions();
      auto layer = rhtools.getLayerWithOffset(thisclusterHits[0].first);

      //std::cout << rhtools.getLayer(thisclusterHits[0].first) << "  w/ offset : " << layer << "  Z : " << (rhtools.getPosition(thisclusterHits[0].first)).z();

      layerClusters_by_layer[layer].push_back(thislayerCluster);
      vertices_by_layer[layer].push_back(i);
      layerClusterenergies.push_back(thislayerCluster.energy());
      std::cout << std::endl;
    }

    auto result_glob = std::max_element(layerClusterenergies.begin(), layerClusterenergies.end());
    auto maxE_vertex = std::distance(layerClusterenergies.begin(), result_glob);
    auto maxE_layerCluster = layerClusters[trackster.vertices(maxE_vertex)];
    auto maxE_layer = rhtools.getLayerWithOffset((maxE_layerCluster.hitsAndFractions())[0].first);

    std::cout << "Max E layer = " << maxE_layer << std::endl;

    std::vector<unsigned> filtered_idx; // higest energy vertices in a layer

    for (unsigned i = 1; i <= rhtools.lastLayer(); ++i) {
      auto layerClusters_in_layer = layerClusters_by_layer[i];
      auto vertices_in_layer = vertices_by_layer[i];

      if (vertices_in_layer.empty()) continue;

      std::cout << " Layer " << i << " : ";
      for (auto i : vertices_in_layer) std::cout << "  " << i;

      std::vector<double> energies_in_layer;
      for (auto j : layerClusters_in_layer) energies_in_layer.push_back(j.energy());
    
      auto result_inlayer = std::max_element(energies_in_layer.begin(),energies_in_layer.end());
      unsigned maxE_idx = std::distance(energies_in_layer.begin(),result_inlayer);

      std::cout << "  - maxE = " << energies_in_layer[maxE_idx] << std::endl;
      
      if ((int)i >= (int)maxE_layer - 10 && (int)i <= (int)maxE_layer + 15)
      filtered_idx.push_back(vertices_in_layer[maxE_idx]);
    }

    std::cout << std::endl;

    // Compute the Covariance Matrix and the sum of the squared weights, used
    // to compute the correct normalization.
    // The barycenter has to be known.
    std::cout << "# before filtering : " << N << std::endl;
    std::cout << "# filtered idx: " << filtered_idx.size() << std::endl;
    for (size_t i : filtered_idx) {
      fillPoint(layerClusters[trackster.vertices(i)]);
      if (energyWeight && trackster.raw_energy())
        weight =
            (layerClusters[trackster.vertices(i)].energy() / trackster.vertex_multiplicity(i)) / trackster.raw_energy();
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
  }
}
#endif