
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
  TH3F* hXY_;
  TTree* tree_;
  float tree_globalX, tree_globalY, tree_globalZ;
  float tree_localX, tree_localY, tree_localZ;
  float tree_energyLoss, tree_time;
  uint32_t tree_detId;


  // ----------member data ---------------------------
   edm::EDGetTokenT<std::vector<PSimHit>> simHitToken_;
  //edm::EDGetTokenT<std::vector<SimVertex>> simHitToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

SimHitAnalyzer::SimHitAnalyzer(const edm::ParameterSet& iConfig)
    :simHitToken_(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simHitTag"))), tkGeomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),     topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>())
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
  else if (simHitsHandle->empty()) {
  edm::LogWarning("SimHitAnalyzer") 
      << "SimHit collection is valid but empty. Event: " << iEvent.id();
  } 
  else {
  edm::LogInfo("SimHitAnalyzer") << "Got " << simHitsHandle->size() << " hits";
  }

  
const TrackerGeometry& tkGeom = iSetup.getData(tkGeomToken_);
const TrackerTopology& tTopo = iSetup.getData(topoToken_);
edm::Handle<std::vector<PSimHit>> hits;
iEvent.getByToken(simHitToken_, hits);
   
for (auto const & hit : *hits) {

     GlobalPoint globalPos = tkGeom.idToDet(hit.detUnitId())->surface().toGlobal(hit.localPosition());	  
    //GlobalPoint globalPos = tkGeom.idToDet(hit.detUnitId())->surface().toGlobal(hit.entryPoint());	  
    edm::LogInfo("HitPosition") << "x=" << globalPos.x() << " y=" << globalPos.y() << " z=" << globalPos.z();
    auto pos = hit.entryPoint();
    //auto mom = hit.momentumAtEntry();
    //double tree_localX  = pos.x(), tree_localY  = pos.y(), tree_localZ  = pos.z();
    tree_globalX = globalPos.x();
    tree_globalY = globalPos.y();
    tree_globalZ = globalPos.z();
    LocalPoint localPos = hit.localPosition();
    tree_localX = localPos.x();
    tree_localY = localPos.y();
    tree_localZ = localPos.z();
    tree_energyLoss = hit.energyLoss();
    tree_time = hit.timeOfFlight();
    const DetId id = hit.detUnitId();
    tree_detId = id.rawId();
    tree_->Fill();
    //globalX.push_back(globalPos.x());
   // globalY.push_back(globalPos.y());
    //globalZ.push_back(globalPos.z());  
    //hXY_->Fill(globalPos.x(), globalPos.y(), globalPos.z());

    hXY_->Fill(pos.x(),pos.y(), pos.z());
    /*double px = mom.x(), py = mom.y(), pz = mom.z();
    edm::LogInfo("SimHitAnalyzer") << "DetUnitId = " << hit.detUnitId()
                                 << ", TOF = " << hit.timeOfFlight()
                                 << ", pAbs = " << hit.pabs()
    				 << ", pos = "  <<  pos
				 << ", mom = " << mom 
				 << ", x = " << x ;*/
    // record or print them
 }
   //TGraph* gXY = new TGraph();
   //for (size_t i = 0; i < globalX.size(); ++i){
   //gXY->SetPoint(i, globalX[i], globalY[i]);
   //}
   //gXY->SetTitle("SimHits (XY);x [cm];y [cm];z [cm]");
   //gXY->Draw("P0");  // or use "PC" to color points

 


/*
  edm::Handle<std::vector<SimVertex>> simVertexHandle;
  //iEvent.getByToken(simVertexToken_, simVertexHandle);
  iEvent.getByToken(simHitToken_, simVertexHandle);

if (simVertexHandle.isValid()) {
    for (const auto& vertex : *simVertexHandle) {
        edm::LogInfo("SimHitAnalyzer")
            << "Vertex ID: " << vertex.vertexId()
            << " Position (cm): " << vertex.position().x() << ", "
            << vertex.position().y() << ", "
            << vertex.position().z();
    }
} else {
    edm::LogWarning("SimHitAnalyzer") << "SimVertex collection not valid!";
}
*/


//const TrackerGeometry& tkGeom = iSetup.getData(tkGeomToken_);
//const TrackerTopology& tTopo = iSetup.getData(topoToken_);
/*
for (const auto& det : tkGeom.dets()) {
    const DetId id = det->geographicalId();
    edm::LogInfo("DetId") << id.rawId() ;
    if (!det) continue;
    if (id.det() == DetId::Tracker) {
        GlobalPoint pos = det->position();

        edm::LogInfo("LayerPosition") 
            << ", Position (cm) = (" 
            << pos.x() << ", " 
            << pos.y() << ", " 
            << pos.z() << ")";

        
         std::ostringstream nameStream;

        if (id.subdetId() == StripSubdetector::TOB) {
            nameStream << "TOB, Layer " << tTopo.tobLayer(id)
                       << ", Rod " << tTopo.tobRod(id)
               	       << ", Module " << tTopo.tobModule(id)
                       << ", Stereo = " << tTopo.tobStereo(id);
	    if (tTopo.isLower(id))  edm::LogInfo("Sensor") << "This is the LOWER sensor (stereo=0)";
            else if (tTopo.isUpper(id)) edm::LogInfo("Sensor") << "This is the UPPER sensor (stereo=1)";
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
  //hXY_ = fs->make<TH2F>("hXY", "SimHit XY positions;X [cm];Y [cm]", 100, -6.5, 100, 200, -100, 100);
  hXY_ = fs->make<TH3F>("hXY", "SimHit XY positions", 100, -6.5, 6.5, 100, -50, 50, 100, -70, 70 );
  tree_ = fs->make<TTree>("SimHitsTree", "SimHits Tree");

  tree_->Branch("globalX", &tree_globalX, "globalX/F");
  tree_->Branch("globalY", &tree_globalY, "globalY/F");
  tree_->Branch("globalZ", &tree_globalZ, "globalZ/F");

  tree_->Branch("localX", &tree_localX, "localX/F");
  tree_->Branch("localY", &tree_localY, "localY/F");
  tree_->Branch("localZ", &tree_localZ, "localZ/F");

  tree_->Branch("energyLoss", &tree_energyLoss, "energyLoss/F");
  tree_->Branch("time", &tree_time, "time/F");

  tree_->Branch("detId", &tree_detId, "detId/i");
}

// ------------ method called once each job just after ending the event loop  ------------
void SimHitAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("simHitTag", edm::InputTag("g4SimHits", "TrackerHitsTIBLowTof"));
  //desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimHitAnalyzer);
