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

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

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
                        const std::vector<reco::Track> &,
                        const StringCutObjectSelector<reco::Track>,
                        const std::vector<CaloParticle> &,
                        const std::vector<Trackster> &,
                        const std::vector<Trackster> &,
                        std::vector<SuperTrackster> &,
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
		
    std::once_flag initializeGeometry_;

    const HGCalDDDConstants* hgcons_;

    std::unique_ptr<GeomDet> firstDisk_[2];

    MagneticField const* bFieldProd_;
    hgcal::RecHitTools rhtools_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;
    
  };
}
#endif
