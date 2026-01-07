import FWCore.ParameterSet.Config as cms
#from Configuration.Geometry.GeometryExtendedRun4D500Reco_cfi import *
from Alignment.CommonAlignmentProducer.GlobalPosition_Fake_cff import *
process = cms.Process("ANALYZE")

#process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring('file:/eos/home-h/hrejebsf/CRACK/step2.root'))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring('file:step2.root'))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
#eventsToProcess = cms.untracked.VEventRange('1:1:629')
process.load("Geometry.CMSCommonData.cmsExtendedGeometryRun4D500XML_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D500Reco_cff")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
#process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_dd4hep_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Alignment.CommonAlignmentProducer.GlobalPosition_Fake_cff")
#process.printGeo = cms.EDAnalyzer("PrintTrackerGeometry")
process.SimHitAnalyzer = cms.EDAnalyzer("SimHitAnalyzer",
    simHitTag = cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof","SIM"),  # Change this tag as needed
    digiHitTag   = cms.InputTag("mix", "Tracker", "DIGI"),
    clusterIncHitTag = cms.InputTag("TTClustersFromPhase2TrackerDigis"   "ClusterInclusive"   "DIGI"), 
    #clusterAcceptHitTag = cms.InputTag("TTClustersFromPhase2TrackerDigis"   "ClusterAccepted"   "DIGI"), 
    #clusterRejectHitTag = cms.InputTag("TTClustersFromPhase2TrackerDigis"   "ClusterRejected"   "DIGI") 
)
#process.XMLFromDD4hep.sourceType = cms.string("DDD")
#process.XMLFromDD4hep.dumpTree = cms.untracked.bool(True)
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO'),  # or 'DEBUG', 'WARNING'
        default   = cms.untracked.PSet(
            limit = cms.untracked.int32(-1)       # unlimited messages
        ),
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("simHits_test.root")
)
process.p = cms.Path(process.SimHitAnalyzer)

