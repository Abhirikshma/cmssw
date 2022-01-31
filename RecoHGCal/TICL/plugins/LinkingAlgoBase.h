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
                                std::vector<Trackster> &tracksters,
                                std::vector<SuperTrackster>& result) = 0;

    static void fillPSetDescription(edm::ParameterSetDescription& desc) {};

  };
}  // namespace ticl

#endif
