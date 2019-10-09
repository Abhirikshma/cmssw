
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalCLUEAlgo.h"

// Geometry
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"

#include "DataFormats/DetId/interface/DetId.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "tbb/task_arena.h"
#include "tbb/tbb.h"
#include <limits>
#include <fstream>
#include <iostream>
#include <ctime>
#include <unordered_map>
#include <map>
#include <google/dense_hash_map>

using namespace hgcal_clustering;
//using std::dense_hash_map;

void HGCalCLUEAlgo::populate(const HGCRecHitCollection &hits) {}
std::unordered_map<uint32_t, int> HGCalCLUEAlgo::populate(const HGCRecHitCollection &hits, std::unordered_map<uint32_t, int> u) {
  // loop over all hits and create the Hexel structure, skip energies below ecut

  if (dependSensor_) {
    // for each layer and wafer calculate the thresholds (sigmaNoise and energy)
    // once
    computeThreshold();
  }
  
  for (unsigned int i = 0; i < hits.size(); ++i) {
    const HGCRecHit &hgrh = hits[i];
    DetId detid = hgrh.detid();
    unsigned int layerOnSide = (rhtools_.getLayerWithOffset(detid) - 1);

    // set sigmaNoise default value 1 to use kappa value directly in case of
    // sensor-independent thresholds
    float sigmaNoise = 1.f;
    if (dependSensor_) {
      int thickness_index = rhtools_.getSiThickIndex(detid);
      if (thickness_index == -1) thickness_index = 3;
      double storedThreshold = thresholds_[layerOnSide][thickness_index];
      sigmaNoise = v_sigmaNoise_[layerOnSide][thickness_index];

      if (hgrh.energy() < storedThreshold)
        continue;  // this sets the ZS threshold at ecut times the sigma noise
                   // for the sensor
    }
    if (!dependSensor_ && hgrh.energy() < ecut_) continue;
    const GlobalPoint position(rhtools_.getPosition(detid));

    int offset = ((rhtools_.zside(detid) + 1) >> 1)*maxlayer;
    int layer = layerOnSide + offset;
    cells_[layer].detid.emplace_back(detid);
    cells_[layer].x.emplace_back(position.x());
    cells_[layer].y.emplace_back(position.y());
    cells_[layer].weight.emplace_back(hgrh.energy());
    cells_[layer].sigmaNoise.emplace_back(sigmaNoise);

    /*if (u.count(detid) == 0){
    u[detid] = (cells_[layer].detid.size()-1);
    }*/
    u.emplace(detid.rawId(), cells_[layer].detid.size()-1);
    
  }
  return u;
}

/*std::unordered_map<DetId, int> makeMap()
{
  std::unordered_map<DetId, int> u = {};
  for (int i = 0; i < 106; i++)
    {
      auto& cellsOnLayer = cells_[i];
      unsigned int numberOfCells = cellsOnLayer.detid.size();
      for (unsigned int n = 0; n < numberOfCells; n++)
	{
	  u[cellsOnLayer.detid[n]] = n;
	}
    }
  std::cout << "unordered-map made\n";
  return u;
  }*/


void HGCalCLUEAlgo::prepareDataStructures(unsigned int l)
{
  auto cellsSize = cells_[l].detid.size();
  cells_[l].rho.resize(cellsSize,0);
  cells_[l].delta.resize(cellsSize,9999999);
  cells_[l].nearestHigher.resize(cellsSize,-1);
  cells_[l].clusterIndex.resize(cellsSize,-1);
  cells_[l].followers.resize(cellsSize);
  cells_[l].isSeed.resize(cellsSize,false);
  
}

// Create a vector of Hexels associated to one cluster from a collection of
// HGCalRecHits - this can be used directly to make the final cluster list -
// this method can be invoked multiple times for the same event with different
// input (reset should be called between events)
void HGCalCLUEAlgo::makeClusters(){}
void HGCalCLUEAlgo::makeClusters(const HGCalTopology& topoEE, const HGCalTopology& topoSi, const HGCalTopology& topoScin) {
  //std::ofstream out;
  //out.open("testCLUE_modrho.csv", std::ofstream::out | std::ofstream::trunc);
  // assign all hits in each layer to a cluster core

  //std::clock_t start;
  //double duration;
  tbb::this_task_arena::isolate([&] {
    tbb::parallel_for(size_t(0), size_t(2 * maxlayer + 2), [&](size_t i) {
	//HGCalLayerTiles lt;
	//lt.fill(cells_[i].x,cells_[i].y);

      prepareDataStructures(i);
      //std::cout << "Data structures prepared for layer "<<i << std::endl;

      /*float delta_c;  // maximum search distance (critical distance) for local
                  // density calculation
      if (i%maxlayer < lastLayerEE)
	{
	  delta_c = vecDeltas_[0];
	}
      else if (i%maxlayer < lastLayerFH)
	{
	  delta_c = vecDeltas_[1];
	}
      else
	{
	  delta_c = vecDeltas_[2];
	  }*/

      /*
      //Using std::map
      std::map<DetId, int> m;
      auto& cellsOnLayer = cells_[i];
      unsigned int numberOfCells = cellsOnLayer.detid.size();
      
      for (unsigned int i = 0; i < numberOfCells; i++)
	{
	  m[cellsOnLayer.detid[i]] = i;
	  }*/

      /*//Using std::unordered_map
      std::unordered_map<DetId, int> u_unordered = {};
      auto& cellsOnLayer = cells_[i];
      unsigned int numberOfCells = cellsOnLayer.detid.size();

      for (unsigned int i = 0; i < numberOfCells; i++)
	{
	  u_unordered[cellsOnLayer.detid[i]] = i;
	  }*/

     
      //Using google::dense_hash_map (fastest)
      struct eqstr
      {
	bool operator() (const uint32_t* s1, const uint32_t* s2) const
	{
	  return (s1 == s2);
	}
      };

      google::dense_hash_map<uint32_t, int> u_google_dense;
      u_google_dense.set_empty_key(-1);

      auto& cellsOnLayer = cells_[i];
      unsigned int numberOfCells = cellsOnLayer.detid.size();

      for (unsigned int i = 0; i < numberOfCells; i++)
	{
	  u_google_dense[cellsOnLayer.detid[i].rawId()] = i;
	}
      
      std::cout << "using google hash" << std::endl;
      
      
      localDensityModFindSeeds(i, topoEE, topoSi, topoScin, u_google_dense);
      numberOfClustersPerLayer_[i] = clusterModified(i, topoEE, topoSi, topoScin, u_google_dense);
      
      //std::cout << "CLUE" << std::endl;
      //calculateLocalDensity(lt, i, delta_c);
      //calculateDistanceToHigher(lt, i, delta_c);
      //numberOfClustersPerLayer_[i] = findAndAssignClusters(i,delta_c);  
      
      /*auto& cellsOnLayer = cells_[i];
      unsigned int numberOfCells = cellsOnLayer.detid.size();
      for (unsigned int m = 0; m < numberOfCells; m++)
	{
	  out << i <<","<< m <<","<<cellsOnLayer.x[m]<<","<<cellsOnLayer.y[m]<<","<<cellsOnLayer.rho[m]<<","<<cellsOnLayer.delta[m]<<","<<cellsOnLayer.nearestHigher[m]<<","<<cellsOnLayer.clusterIndex[m]<<std::endl; 
	  }*/


    });
  });
  //Now that we have the density per point we can store it
  //for(unsigned int i=0; i< 2 * maxlayer + 2; ++i) { setDensity(i);}
  //duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
  //std::cout <<" Time for all layers = " << duration << "\n";
}

std::vector<reco::BasicCluster> HGCalCLUEAlgo::getClusters(bool) {

  std::vector<int> offsets(numberOfClustersPerLayer_.size(),0);

  int maxClustersOnLayer = numberOfClustersPerLayer_[0];

  for(unsigned layerId = 1; layerId<offsets.size(); ++layerId)
  {
    offsets[layerId] = offsets[layerId-1] +numberOfClustersPerLayer_[layerId-1];

    maxClustersOnLayer = std::max(maxClustersOnLayer, numberOfClustersPerLayer_[layerId]);
  }


  auto totalNumberOfClusters = offsets.back()+numberOfClustersPerLayer_.back();
  clusters_v_.resize(totalNumberOfClusters);
  std::vector<std::vector<int> > cellsIdInCluster;
  cellsIdInCluster.reserve(maxClustersOnLayer);

  for(unsigned int layerId = 0; layerId < 2 * maxlayer + 2; ++layerId)
  {
    cellsIdInCluster.resize(numberOfClustersPerLayer_[layerId]);
    auto& cellsOnLayer = cells_[layerId];
    unsigned int numberOfCells = cellsOnLayer.detid.size();
    auto firstClusterIdx = offsets[layerId];
    
    for (unsigned int i = 0; i < numberOfCells; ++i )
    {   
      auto clusterIndex = cellsOnLayer.clusterIndex[i];
      if(clusterIndex != -1)
        cellsIdInCluster[clusterIndex].push_back(i);
    }
    

    std::vector<std::pair<DetId, float>> thisCluster;
    
    for(auto& cl: cellsIdInCluster)
    {
      auto position = calculatePosition(cl, layerId);
      float energy = 0.f;
      int seedDetId = -1;
      
      for(auto cellIdx : cl)
      {
        energy+= cellsOnLayer.weight[cellIdx];
        thisCluster.emplace_back(cellsOnLayer.detid[cellIdx],1.f);
        if(cellsOnLayer.isSeed[cellIdx])
        {
          seedDetId = cellsOnLayer.detid[cellIdx];
        }
      }
      auto globalClusterIndex = cellsOnLayer.clusterIndex[cl[0]] +  firstClusterIdx;
      
      clusters_v_[globalClusterIndex]=reco::BasicCluster(energy, position, reco::CaloID::DET_HGCAL_ENDCAP, thisCluster, algoId_);
      clusters_v_[globalClusterIndex].setSeed(seedDetId);
      thisCluster.clear();
    }

    cellsIdInCluster.clear();

  }

  return clusters_v_;

}



math::XYZPoint HGCalCLUEAlgo::calculatePosition(const std::vector<int> &v, const unsigned int layerId) const {
  
  float total_weight = 0.f;
  float x = 0.f;
  float y = 0.f;

  unsigned int maxEnergyIndex = 0;
  float maxEnergyValue = 0.f;
  
  auto& cellsOnLayer = cells_[layerId];


  // loop over hits in cluster candidate
  // determining the maximum energy hit
  for (auto i : v) {
    total_weight += cellsOnLayer.weight[i];
    if (cellsOnLayer.weight[i] > maxEnergyValue) {
      maxEnergyValue = cellsOnLayer.weight[i];
      maxEnergyIndex = i;
    }
  }

  // Si cell or Scintillator. Used to set approach and parameters
  auto thick = rhtools_.getSiThickIndex(cellsOnLayer.detid[maxEnergyIndex]);
  bool isSiliconCell = (thick != -1);


// TODO: this is recomputing everything twice and overwriting the position with log weighting position
  if(isSiliconCell)
  {
    float total_weight_log = 0.f;
    float x_log = 0.f;
    float y_log = 0.f;
    for (auto i : v) {
      float rhEnergy = cellsOnLayer.weight[i];
      float Wi = std::max(thresholdW0_[thick] + std::log(rhEnergy / total_weight), 0.);
      x_log += cellsOnLayer.x[i] * Wi;
      y_log += cellsOnLayer.y[i] * Wi;
      total_weight_log += Wi;
    }

    total_weight = total_weight_log;
    x = x_log;
    y = y_log;
  }
  else
  {
    for (auto i : v) {

      float rhEnergy = cellsOnLayer.weight[i];

      x += cellsOnLayer.x[i] * rhEnergy;
      y += cellsOnLayer.y[i] * rhEnergy;
      
    }
  }
  if (total_weight != 0.) {
    float inv_tot_weight = 1.f / total_weight;
    return math::XYZPoint(x * inv_tot_weight, y * inv_tot_weight, rhtools_.getPosition(cellsOnLayer.detid[maxEnergyIndex]).z());
  }
  else return math::XYZPoint(0.f, 0.f, 0.f);
}



void HGCalCLUEAlgo::calculateLocalDensity(const HGCalLayerTiles& lt, const unsigned int layerId, float delta_c)  
{

  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();

  for(unsigned int i = 0; i < numberOfCells; i++) 
  {
    std::array<int,4> search_box = lt.searchBox(cellsOnLayer.x[i] - delta_c, cellsOnLayer.x[i] + delta_c, cellsOnLayer.y[i] - delta_c, cellsOnLayer.y[i] + delta_c);
    
    for(int xBin = search_box[0]; xBin < search_box[1]+1; ++xBin) {
      for(int yBin = search_box[2]; yBin < search_box[3]+1; ++yBin) {
        
        int binId = lt.getGlobalBinByBin(xBin,yBin);
        size_t binSize = lt[binId].size();
        
        for (unsigned int j = 0; j < binSize; j++) {
          unsigned int otherId = lt[binId][j];
          if(distance(i, otherId, layerId) < delta_c) {
            cellsOnLayer.rho[i] += (i == otherId ? 1.f : 0.5f) * cellsOnLayer.weight[otherId];
          }
        }
      }
    }    
  }

}

//rho = sum of energy in a hexagonal grid around a cell
void HGCalCLUEAlgo::localDensityModFindSeeds(const unsigned int layerId, const HGCalTopology& topoEE, const HGCalTopology& topoSi, const HGCalTopology& topoScin, google::dense_hash_map<uint32_t, int> u)
{
  auto& cellsOnLayer = cells_[layerId];
  std::vector<DetId> neighborCells;
  float sumWeights = 0.f;
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  //unsigned int nClustersOnLayer = 0;
  float rho_c;

  std::cout << "Calculating local density and finding seeds!" << std::endl;
  //std::cout << "Total #cells =" << numberOfCells << std::endl;


  for(unsigned int i = 0; i < numberOfCells; i++)
    {
      DetId detid = cellsOnLayer.detid[i];
      uint32_t rawId;

      switch (detid.det()) {
      case DetId::HGCalEE:
	neighborCells = topoEE.neighbors(cellsOnLayer.detid[i]);
	break;
      case DetId::HGCalHSi:
	neighborCells = topoSi.neighbors(cellsOnLayer.detid[i]);
	break;
      case DetId::HGCalHSc:
	neighborCells = topoScin.neighbors(cellsOnLayer.detid[i]);
	break;
      default:
	std::cout<< "DetId does not belong to any! \n";
       
      }

      rho_c = kappa_ * cellsOnLayer.sigmaNoise[i];
      //neighborCells = topo.neighbors(cellsOnLayer.detid[i]);
      sumWeights = 0.f;

      /*for(auto j : neighborCells) 
	{
	  for(unsigned int k = 0; k < numberOfCells; k++) 
	    {
	      if (cellsOnLayer.detid[k] == j) 
		{
		  sumWeights += 0.5*cellsOnLayer.weight[k];
		  break;
		}
	    }
	    }*/
      
      for (auto j : neighborCells)
	{
	  rawId = j.rawId();
	  sumWeights += 0.5*cellsOnLayer.weight[u[rawId]];
	}
      cellsOnLayer.rho[i] = sumWeights + cellsOnLayer.weight[i];
      if (cellsOnLayer.rho[i] > 0.8*rho_c)
	{
	  cellsOnLayer.isSeed[i] = true;
	  //cellsOnLayer.clusterIndex[i] = nClustersOnLayer;
	  //nClustersOnLayer++;
	  //std::cout << "Cell " << i << " is seed!" << std::endl;
	}
      //else cellsOnLayer.clusterIndex[i] = -1;
      cellsOnLayer.followers[i].clear();
  
    }

  std::cout << "local density done!" << std::endl;
}

/*void HGCalCLUEAlgo::clusterModified(const unsigned int layerId)
{
  unsigned int nClustersOnLayer = 0;
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  float clusterEnergy = 0.f;
  std::vector<DetID> neighbors;
 

  for(i = 0; i < numberOfCells; i++)
    {
      if (cellsOnLayer.isSeed[i] == true) 
	{
	  //find out if it has atleast one neighbour that has non zero energy otherwise mark as outlier(?)
	  //Assign clusterIndex to neighbours
	  clusterEnergy = 0.f;
	  neighbors = HGCalTopology.neighbors (cellsOnLayer.detid[i]);// does it return the original cell?
	  for(auto j : neighbors) 
	    {
	      for(unsigned int k = 0; k < numberOfCells; k++)
		{
		  if (cellsOnLayer.detid[k] == j) //if cell is in the neighbourhood
		    {
		      clusterEnergy += cellsOnLayer.weight[k];
		      cellsOnLayer.clusterIndex[k] = cellsOnLayer.clusterIndex[i];
		      if (cellsOnLayer.isSeed[k] == true) //if seed in neighbourhood
			{
			  neighborsSecondary = HGCalTopology.neighbors (cellsOnLayer.detid[k]);
			  for(auto m : neighborsSecondary)
			    { 
			      for(unsigned int l = 0; l < numberOfCells; l++)
				{
				  if (cellsOnLayer.detid[l] == m)
				    {
				      clusterEnergy += cellsOnlayer.weight[l];
				      cellsOnLayer.clusterIndex[l] = cellsOnLayer.clusterIndex[i];
				    }   
				     
				}
			    }
			}
		    }
		}
	    }//one cluster done	    
	
	}
    }
    }*/

int HGCalCLUEAlgo::clusterModified(const unsigned int layerId, const HGCalTopology& topoEE, const HGCalTopology& topoSi, const HGCalTopology& topoScin,google::dense_hash_map<uint32_t, int> u)
{
  unsigned int nClustersOnLayer = 0;
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  std::vector<DetId> neighborCells;
  std::vector<int> localStack;
  int highestNeighbor = 0;
  float highestRho = 0.f;

  std::cout << "Clustering starts" << std::endl;


  for(unsigned int i = 0; i < numberOfCells; i++)
    {

      DetId detid = cellsOnLayer.detid[i]; 
      uint32_t rawId;
      //int u_coord = rhtools_.getCell(detid).first;
      //int v_coord = rhtools_.getCell(detid).second;

      switch (detid.det()) {
      case DetId::HGCalEE:
	neighborCells = topoEE.neighbors(cellsOnLayer.detid[i]);
	break;
      case DetId::HGCalHSi:
	neighborCells = topoSi.neighbors(cellsOnLayer.detid[i]);
	break;
      case DetId::HGCalHSc:
	neighborCells = topoScin.neighbors(cellsOnLayer.detid[i]);
	break;
      default:
	std::cout<<"DetId does not belong to any! \n";
       
      }

      highestNeighbor = 0;
      highestRho = 0.f;
      float thisRho = 0.f;
      //bool isOutlier = (cellsOnLayer.rho[i] == cellsOnLayer.weight[i]);
      
      //if (!isOutlier){
      //neighborCells = topo.neighbors(cellsOnLayer.detid[i]);
      
      /*for(auto j : neighborCells)
	{
	  for (unsigned int k =0; k < numberOfCells; k++)
	    {
	      if (cellsOnLayer.detid[k] == j)
		{
		  if (cellsOnLayer.rho[k] > highestRho)
		    {
		      highestNeighbor = k;
		      highestRho = cellsOnLayer.rho[k];
		    }
		  break;
		}
	    }
	    }*/

      for (auto j : neighborCells)
	{
	  rawId = j.rawId();
	  thisRho = cellsOnLayer.rho[u[rawId]]; 
	  if (thisRho > highestRho)
	    {
	      highestNeighbor = u[rawId];
	      highestRho = thisRho;
	    }
	}
      if (highestRho > cellsOnLayer.rho[i])
	{
	  if (cellsOnLayer.isSeed[i])
	    {
	      //cellsOnLayer.isSeed[highestNeighbor] = true;
	      //cellsOnLayer.clusterIndex[highestNeighbor] = cellsOnLayer.clusterIndex[i];
	      cellsOnLayer.isSeed[i] = false;
	      //cellsOnLayer.clusterIndex[i] = -1;
	    }
	  cellsOnLayer.followers[highestNeighbor].push_back(i);
	  cellsOnLayer.nearestHigher[i] = highestNeighbor;
	}
      else cellsOnLayer.nearestHigher[i] = -1;
      // }
      
    }

  for (unsigned int i=0; i < numberOfCells; i++)
    {
      cellsOnLayer.clusterIndex[i] = -1;
      if (cellsOnLayer.isSeed[i])
	{
	  cellsOnLayer.clusterIndex[i] = nClustersOnLayer;
	  nClustersOnLayer++;
	  localStack.push_back(i);
	}
    }


  std::cout << "# elements in localStack = " <<localStack.size() << "\n";

  while (!localStack.empty()) 
    {
      int endStack = localStack.back();
      //std::cout << "\n endStack = " << endStack << "\n";
      auto& thisSeed = cellsOnLayer.followers[endStack];
      //std::cout << "#followers = " << thisSeed.size() << "\n";
      //std::cout << "followers = "; 
      localStack.pop_back();
      if(thisSeed.empty()) continue;
	  

      // loop over followers
      for( int l : thisSeed)
	{
	  // pass id to a follower
	  cellsOnLayer.clusterIndex[l] = cellsOnLayer.clusterIndex[endStack];
	  // push this follower to localStack
	  //int z = 0;
	  //int isThere = 0;
	  /*for (int m : localStack)
	    {
	      if (l==m) 
		{
		  //const int itr = z;
		  localStack.erase(localStack.begin() + z);
		  isThere = 1;
		  break;
		}
	      z++;
	    }
	    if (!isThere) localStack.push_back(l);*/
	  //if (localStack.back() != l) localStack.push_back(l);
	  //else localStack.pop_back();
	  localStack.push_back(l);
	}
      
    }
  std::cout << "local stack empty!" << std::endl;

    
  return nClustersOnLayer;
}


void HGCalCLUEAlgo::calculateDistanceToHigher(const HGCalLayerTiles& lt, const unsigned int layerId, float delta_c) {


  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();

  for(unsigned int i = 0; i < numberOfCells; i++) {
    // initialize delta and nearest higher for i
    float maxDelta = std::numeric_limits<float>::max();
    float i_delta = maxDelta;
    int i_nearestHigher = -1;

    // get search box for ith hit
    // guarantee to cover a range "outlierDeltaFactor_*delta_c"
    auto range = outlierDeltaFactor_*delta_c;
    std::array<int,4> search_box = lt.searchBox(cellsOnLayer.x[i]  - range, cellsOnLayer.x[i] + range, cellsOnLayer.y[i] - range, cellsOnLayer.y[i] + range);
    
    // loop over all bins in the search box
    for(int xBin = search_box[0]; xBin < search_box[1]+1; ++xBin) {
      for(int yBin = search_box[2]; yBin < search_box[3]+1; ++yBin) {
        
        // get the id of this bin
        size_t binId = lt.getGlobalBinByBin(xBin,yBin);
        // get the size of this bin
        size_t binSize = lt[binId].size();

        // loop over all hits in this bin
        for (unsigned int j = 0; j < binSize; j++) {
          unsigned int otherId = lt[binId][j];

          float dist = distance(i, otherId, layerId);
          bool foundHigher = cellsOnLayer.rho[otherId] > cellsOnLayer.rho[i];


          // if dist == i_delta, then last comer being the nearest higher
          if(foundHigher && dist <= i_delta) {

            // update i_delta
            i_delta = dist;
            // update i_nearestHigher
            i_nearestHigher = otherId;
            
          }
        }
      }
    }

    bool foundNearestHigherInSearchBox = (i_delta != maxDelta);
    if (foundNearestHigherInSearchBox){
      cellsOnLayer.delta[i] = i_delta;
      cellsOnLayer.nearestHigher[i] = i_nearestHigher;
    } else {
      // otherwise delta is guaranteed to be larger outlierDeltaFactor_*delta_c
      // we can safely maximize delta to be maxDelta
      cellsOnLayer.delta[i] = maxDelta;
      cellsOnLayer.nearestHigher[i] = -1;
    }
  }
}


int HGCalCLUEAlgo::findAndAssignClusters(const unsigned int layerId, float delta_c ) {
   
  // this is called once per layer and endcap...
  // so when filling the cluster temporary vector of Hexels we resize each time
  // by the number  of clusters found. This is always equal to the number of
  // cluster centers...
  unsigned int nClustersOnLayer = 0;
  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  std::vector<int> localStack;
  // find cluster seeds and outlier  
  for(unsigned int i = 0; i < numberOfCells; i++) {
    float rho_c = kappa_ * cellsOnLayer.sigmaNoise[i];
    // initialize clusterIndex
    cellsOnLayer.clusterIndex[i] = -1;
    bool isSeed = (cellsOnLayer.delta[i] > delta_c) && (cellsOnLayer.rho[i] >= rho_c);
    bool isOutlier = (cellsOnLayer.delta[i] > outlierDeltaFactor_*delta_c) && (cellsOnLayer.rho[i] < rho_c);
    if (isSeed) 
    {
      cellsOnLayer.clusterIndex[i] = nClustersOnLayer;
      cellsOnLayer.isSeed[i] = true;
      nClustersOnLayer++;
      localStack.push_back(i);
    
    } else if (!isOutlier) {
      cellsOnLayer.followers[cellsOnLayer.nearestHigher[i]].push_back(i);   
    } 
  }
  // need to pass clusterIndex to their followers
  while (!localStack.empty()) {
    int endStack = localStack.back();
    auto& thisSeed = cellsOnLayer.followers[endStack];
    localStack.pop_back();

    // loop over followers
    for( int j : thisSeed){
      // pass id to a follower
      cellsOnLayer.clusterIndex[j] = cellsOnLayer.clusterIndex[endStack];
      // push this follower to localStack
      localStack.push_back(j);
    }
    
  }
  return nClustersOnLayer;
}

void HGCalCLUEAlgo::computeThreshold() {
  // To support the TDR geometry and also the post-TDR one (v9 onwards), we
  // need to change the logic of the vectors containing signal to noise and
  // thresholds. The first 3 indices will keep on addressing the different
  // thicknesses of the Silicon detectors, while the last one, number 3 (the
  // fourth) will address the Scintillators. This change will support both
  // geometries at the same time.

  if (initialized_) return;  // only need to calculate thresholds once

  initialized_ = true;

  std::vector<double> dummy;
  const unsigned maxNumberOfThickIndices = 3;
  dummy.resize(maxNumberOfThickIndices + 1, 0);  // +1 to accomodate for the Scintillators
  thresholds_.resize(maxlayer, dummy);
  v_sigmaNoise_.resize(maxlayer, dummy);

  for (unsigned ilayer = 1; ilayer <= maxlayer; ++ilayer) {
    for (unsigned ithick = 0; ithick < maxNumberOfThickIndices; ++ithick) {
      float sigmaNoise = 0.001f * fcPerEle_ * nonAgedNoises_[ithick] * dEdXweights_[ilayer] /
                         (fcPerMip_[ithick] * thicknessCorrection_[ithick]);
      thresholds_[ilayer - 1][ithick] = sigmaNoise * ecut_;
      v_sigmaNoise_[ilayer - 1][ithick] = sigmaNoise;
    }
    float scintillators_sigmaNoise = 0.001f * noiseMip_ * dEdXweights_[ilayer];
    thresholds_[ilayer - 1][maxNumberOfThickIndices] = ecut_ * scintillators_sigmaNoise;
    v_sigmaNoise_[ilayer - 1][maxNumberOfThickIndices] = scintillators_sigmaNoise;
  }
}

void HGCalCLUEAlgo::setDensity(const unsigned int layerId){

  auto& cellsOnLayer = cells_[layerId];
  unsigned int numberOfCells = cellsOnLayer.detid.size();
  for (unsigned int i = 0; i< numberOfCells; ++i) density_[ cellsOnLayer.detid[i] ] =   cellsOnLayer.rho[i] ;
  
}

Density HGCalCLUEAlgo::getDensity() {
  return density_;
}
