
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <fstream>  

#include "TH3F.h"
#include "TTree.h"
#include "TGraph2D.h"

class SimHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimHitAnalyzer(const edm::ParameterSet&);
  ~SimHitAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  std::ofstream csv_out_;
  TH3F* hXY_;
  TTree* tree_;
  std::vector<float> tree_globalX, tree_globalY, tree_globalZ;
  std::vector<float> tree_localX, tree_localY, tree_localZ;
  std::vector<float> tree_energyLoss;
  std::vector<float> tree_time;
  std::vector<unsigned int> tree_detId;
  std::vector<unsigned int> tree_nDigis_perDetId;
  int tree_event;
  int tree_nHits;
  std::vector<float> tree_occupency;
  std::vector<unsigned int> tree_Digis_rows;
  std::vector<unsigned int> tree_Digis_columns;

  typedef edm::Ref<
    edm::DetSetVector<Phase2TrackerDigi>,
    Phase2TrackerDigi,
    edm::refhelper::FindForDetSetVector<Phase2TrackerDigi>
> Ref_Phase2TrackerDigi_;
  // ----------member data ---------------------------
   edm::EDGetTokenT<std::vector<PSimHit>> simHitToken_;
  //edm::EDGetTokenT<std::vector<SimVertex>> simHitToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  edm::EDGetTokenT<edm::DetSetVector<Phase2TrackerDigi>> digiToken_;
  //edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<edm::Ref<edm::DetSetVector<Phase2TrackerDigi>,Phase2TrackerDigi,edm::refhelper::FindForDetSetVector<Phase2TrackerDigi>>>> ttClusterIncToken_;
  /*
  edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> ttClusterIncToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> ttClusterAcceptToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> ttClusterRejectToken_;
  */

  //edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>> ttClusterIncToken_;
  //dm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>> ttClusterAcceptToken_;
  //edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>> ttClusterRejectToken_;
  //edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>> ttClusterIncToken_;
  //edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>> ttClusterAcceptToken_;
  //m::EDGetTokenT<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>> ttClusterRejectToken_;

  //edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<edm::Ref<edm::DetSetVector<Phase2TrackerDigi>,Phase2TrackerDigi,edm::refhelper::FindForDetSetVector<Phase2TrackerDigi>>>> ttClusterAcceptToken_;
  //edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<edm::Ref<edm::DetSetVector<Phase2TrackerDigi>,Phase2TrackerDigi,edm::refhelper::FindForDetSetVector<Phase2TrackerDigi>>>> ttClusterRejectToken_;



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

SimHitAnalyzer::SimHitAnalyzer(const edm::ParameterSet& iConfig)
    :simHitToken_(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simHitTag"))),
     tkGeomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
     topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
     digiToken_(consumes<edm::DetSetVector<Phase2TrackerDigi>>(iConfig.getParameter<edm::InputTag>("digiHitTag") ))
     /*
     ttClusterIncToken_ (consumes<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> (iConfig.getParameter<edm::InputTag>("clusterIncHitTag"))),
     ttClusterAcceptToken_ (consumes<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> (iConfig.getParameter<edm::InputTag>("clusterAcceptHitTag"))),
     ttClusterRejectToken_ (consumes<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> (iConfig.getParameter<edm::InputTag>("clusterRejectHitTag")))
     
     ttClusterIncToken_(consumes<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>>(
        iConfig.getParameter<edm::InputTag>("clusterIncHitTag"))),
    ttClusterAcceptToken_(consumes<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>>(
        iConfig.getParameter<edm::InputTag>("clusterAcceptHitTag"))),
    ttClusterRejectToken_(consumes<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>>(
        iConfig.getParameter<edm::InputTag>("clusterRejectHitTag")))
    */
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

SimHitAnalyzer::~SimHitAnalyzer() {
}


// ------------ method called for each event  ------------
void SimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //if( iEvent.id().event() != 629) return ;  
  std::vector<float> globalX, globalY, globalZ;
  using namespace edm;
  edm::LogVerbatim("SimHitAnalyzer") << "Analyzing event: " << iEvent.id();
  edm::Handle<std::vector<PSimHit>> simHitsHandle;
  iEvent.getByToken(simHitToken_, simHitsHandle);

  if (!simHitsHandle.isValid()) {
  edm::LogWarning("SimHitAnalyzer") 
      << "SimHit handle is not valid. Event: " << iEvent.id();
  } 
  /*
  else if (simHitsHandle->empty()) {
  edm::LogWarning("SimHitAnalyzer") 
      << "SimHit collection is valid but empty. Event: " << iEvent.id();
  } 
  else {
  edm::LogInfo("SimHitAnalyzer") << "Got " << simHitsHandle->size() << " hits";
  }*/

tree_globalX.clear();
tree_globalY.clear();
tree_globalZ.clear();
tree_localX.clear();
tree_localY.clear();
tree_localZ.clear();
tree_energyLoss.clear();
tree_time.clear();
tree_detId.clear();
tree_nDigis_perDetId.clear();
tree_occupency.clear();
tree_Digis_rows.clear();
tree_Digis_columns.clear();
  
const TrackerGeometry& tkGeom = iSetup.getData(tkGeomToken_);
const TrackerTopology& tTopo = iSetup.getData(topoToken_);


edm::Handle<std::vector<PSimHit>> hits;
iEvent.getByToken(simHitToken_, hits);
tree_event = iEvent.id().event(); 
tree_nHits = hits->size();
for (auto const & hit : *hits) {

     GlobalPoint globalPos = tkGeom.idToDet(hit.detUnitId())->surface().toGlobal(hit.localPosition());	 
     LocalPoint localPos = hit.localPosition(); 
    edm::LogInfo("HitPosition") << "x=" << globalPos.x() << " y=" << globalPos.y() << " z=" << globalPos.z();
    auto pos = hit.entryPoint();
    tree_globalX.push_back(globalPos.x());
    tree_globalY.push_back(globalPos.y());
    tree_globalZ.push_back(globalPos.z());
    tree_localX.push_back(localPos.x());
    tree_localY.push_back(localPos.y());
    tree_localZ.push_back(localPos.z());
    tree_energyLoss.push_back(hit.energyLoss());
    tree_time.push_back(hit.timeOfFlight());
    tree_detId.push_back(hit.detUnitId());


    hXY_->Fill(pos.x(),pos.y(), pos.z());
 }

   //tree_->Fill();

edm::Handle<edm::DetSetVector<Phase2TrackerDigi>> digis;
iEvent.getByToken(digiToken_, digis);

if (digis.isValid()) {
    // DIGIs are stored per Det Id : so you have a collection of Digis per DetIDs
    // for loop over detIds
    for (const auto& DSV : *digis) {
        uint32_t detid = DSV.id;
	tree_nDigis_perDetId.push_back(DSV.size());

        for (const auto& digi : DSV) {
	     //Loop of over digis per a DetId. 
            // Example: print or fill histos
            edm::LogInfo("Phase2TrackerDigi") << "DetId " << detid
                                              << " row=" << digi.row()
                                              << " col=" << digi.column();	    
	    tree_Digis_rows.push_back(digi.row());
	    tree_Digis_columns.push_back(digi.column());

        }
    }
}
/*
edm::Handle<edmNew::DetSetVector<TTCluster<Phase2TrackerDigi>>> ttClusters;
iEvent.getByToken(ttClusterIncToken_, ttClusters);


if (ttClusters.isValid()) {
    for (const auto& DSV : *ttClusters) {
        uint32_t detid = DSV.id();
        edm::LogInfo("TTCluster") << "DetId " << detid
                                  << " has " << DSV.size() << " clusters";
        for (const auto& cluster : DSV) {
            edm::LogInfo("TTCluster") << "  Cluster size = " << cluster.getHits().size();
        }
    }
}
*/

tree_->Fill();


/*
for (const auto& det : tkGeom.dets()) {
    const DetId id = det->geographicalId();
    edm::LogInfo("DetId") << id.rawId() ;
    std::cout << "test test test "<< std::endl;
    if (!det) continue;
    if (id.det() == DetId::Tracker) {
        GlobalPoint pos = det->position();

        edm::LogInfo("LayerPosition") 
            << ", Position (cm) = (" 
            << pos.x() << ", " 
            << pos.y() << ", " 
            << pos.z() << ")";

          
         std::ostringstream nameStream;
         edm::LogInfo("SubDet Id is, ") << id.subdetId() ; 
        if (id.subdetId() == StripSubdetector::TOB) {
            nameStream << "TOB, Layer " << tTopo.tobLayer(id)
                       << ", Rod " << tTopo.tobRod(id)
               	       << ", Module " << tTopo.tobModule(id)
                       << ", Stereo = " << tTopo.tobStereo(id);
	    if (tTopo.isLower(id))  edm::LogInfo("Sensor") << "This is the LOWER sensor (stereo=0)";
            else if (tTopo.isUpper(id)) edm::LogInfo("Sensor") << "This is the UPPER sensor (stereo=1)";


	    if(!tTopo.isLower(id) && !tTopo.isUpper(id)){
               // Save this DetId in the CSV
               csv_out_ << id.rawId() << "\n";
            }

        } else if (id.subdetId() == PixelSubdetector::PixelBarrel) {
            nameStream << "PixelBarrel, Layer " << tTopo.pxbLayer(id)
                       << ", Ladder " << tTopo.pxbLadder(id);
        } else if (id.subdetId() == StripSubdetector::TID) {
            nameStream << "TID, Side " << tTopo.tidSide(id)
                       << ", Wheel " << tTopo.tidWheel(id);
        } else {
            nameStream << "Unknown subdet = " << id.subdetId();
        }

        edm::LogInfo("DetDecoded") << "Decoded: " << nameStream.str();
        
    }
    
}
*/

}

// ------------ method called once each job just before starting event loop  ------------
void SimHitAnalyzer::beginJob() {
  // please remove this method if not needed
  edm::Service<TFileService> fs;
  csv_out_.open("detid_map.csv");
  csv_out_ << "DetId\n"; 

  //hXY_ = fs->make<TH2F>("hXY", "SimHit XY positions;X [cm];Y [cm]", 100, -6.5, 100, 200, -100, 100);
  hXY_ = fs->make<TH3F>("hXY", "SimHit XY positions", 100, -6.5, 6.5, 100, -50, 50, 100, -70, 70 );
  tree_ = fs->make<TTree>("SimHitsTree", "SimHits Tree");

  tree_->Branch("globalX", &tree_globalX);
  tree_->Branch("globalY", &tree_globalY);
  tree_->Branch("globalZ", &tree_globalZ);

  tree_->Branch("localX", &tree_localX);
  tree_->Branch("localY", &tree_localY);
  tree_->Branch("localZ", &tree_localZ);

  tree_->Branch("energyLoss", &tree_energyLoss);
  tree_->Branch("time", &tree_time);

  tree_->Branch("detId", &tree_detId);
  tree_->Branch("Event", &tree_event);
  tree_->Branch("nHits", &tree_nHits);
  tree_->Branch("nDigis_perDetId", &tree_nDigis_perDetId);
  tree_->Branch("occupency", &tree_occupency);
  tree_->Branch("Digis_columns", &tree_Digis_columns);
  tree_->Branch("Digis_rows", &tree_Digis_rows);
}

// ------------ method called once each job just after ending the event loop  ------------
void SimHitAnalyzer::endJob() {
  // please remove this method if not needed
  csv_out_.close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("simHitTag", edm::InputTag("g4SimHits", "TrackerHitsTIBLowTof"));
  desc.add<edm::InputTag>("digiHitTag", edm::InputTag("mix", "Tracker", "DIGI"));
  desc.add<edm::InputTag>("clusterIncHitTag", edm::InputTag("TTClustersFromPhase2TrackerDigis"   "ClusterInclusive"   "DIGI") ) ; 
  //desc.add<edm::InputTag>("clusterAcceptHitTag", edm::InputTag("TTClustersFromPhase2TrackerDigis"   "ClusterAccepted"   "DIGI") ) ; 
  //desc.add<edm::InputTag>("clusterRejectHitTag", edm::InputTag("TTClustersFromPhase2TrackerDigis"   "ClusterRejected"   "DIGI") ) ; 

  //desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimHitAnalyzer);
