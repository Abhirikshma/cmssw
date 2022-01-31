#ifndef RecoHGCal_TICL_LinkingAlgoByPCAGeometric_H__
#define RecoHGCal_TICL_LinkingAlgoByPCAGeometric_H__

#include <memory>
#include <vector>
#include <string>
#include "RecoHGCal/TICL/plugins/LinkingAlgoBase.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

//#include "RecoHGCal/TICL/test/PCA_mod.h"


namespace ticl {
  class LinkingAlgoByPCAGeometric final : public LinkingAlgoBase {
  public:
		LinkingAlgoByPCAGeometric(const edm::ParameterSet& conf, edm::ConsumesCollector& sumes);
		~LinkingAlgoByPCAGeometric() override;

    void initialize(const edm::EventSetup& es) override;

		void linkTracksters(const edm::Event &, 
                        const edm::EventSetup &,
                        std::vector<Trackster> &,
                        std::vector<SuperTrackster> &) override;

		static void fillPSetDescription(edm::ParameterSetDescription& desc);

	private:
		typedef Trackster::IterationIndex TracksterIterIndex;
  	typedef math::XYZVector Vector;

		void buildFirstLayers();
		//bool isCP(const ticl::Trackster& t, edm::Handle<std::vector<CaloParticle>>);

		const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
    edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdcToken_;
    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldToken_;
    const std::string detector_;
    const std::string propName_;
    const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorToken_;
		//const edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersMergeToken_;
		const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTSToken_;
		const edm::EDGetTokenT<reco::TrackCollection> trackColToken_;
		const edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
		const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClustersToken_;

    std::once_flag initializeGeometry_;

    const HGCalDDDConstants* hgcons_;
    const StringCutObjectSelector<reco::Track> cutTk_;

    std::unique_ptr<GeomDet> firstDisk_[2];

    MagneticField const* bFieldProd_;
    hgcal::RecHitTools rhtools_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;
    
  };
}
#endif
