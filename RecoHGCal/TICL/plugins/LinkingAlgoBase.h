#ifndef RecoHGCal_TICL_LinkingAlgoBase_H__
#define RecoHGCal_TICL_LinkingAlgoBase_H__

#include <memory>
#include <vector>
#include "DataFormats/HGCalReco/interface/Common.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/SuperTrackster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace ticl {
  class LinkingAlgoBase {
  public:
    LinkingAlgoBase(const edm::ParameterSet& conf, edm::ConsumesCollector& sumes) {}
    
    virtual ~LinkingAlgoBase(){};

    virtual void initialize(const edm::EventSetup& es) = 0;

    virtual void linkTracksters(const edm::Event& evt,
                                const edm::EventSetup& es,
                                const std::vector<reco::Track>& tracks,
                                const StringCutObjectSelector<reco::Track> cutTk,
                                const std::vector<CaloParticle>& caloParticles,
                                const std::vector<Trackster> &tracksters,
                                const std::vector<Trackster> &simTracksters,
                                std::vector<SuperTrackster>& resultTracksters,
                                std::vector<SuperTrackster>& resultSimTracksters) = 0;

    static void fillPSetDescription(edm::ParameterSetDescription& desc) {};

  };
}  // namespace ticl

#endif
