import FWCore.ParameterSet.Config as cms
#from Configuration.Geometry.GeometryExtendedRun4D500Reco_cfi import *
from Alignment.CommonAlignmentProducer.GlobalPosition_Fake_cff import *
process = cms.Process("ANALYZE")

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring('file:step1_cosmic.root'))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(200000))
#eventsToProcess = cms.untracked.VEventRange('1:1:629')
process.load("Geometry.CMSCommonData.cmsExtendedGeometryRun4D500XML_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D500Reco_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Alignment.CommonAlignmentProducer.GlobalPosition_Fake_cff")
#process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.printGeo = cms.EDAnalyzer("PrintTrackerGeometry")
process.simHitAnalyzer = cms.EDAnalyzer("SimHitAnalyzer",
    simHitTag = cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof","SIM")  # Change this tag as needed
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
    fileName = cms.string("simHits.root")
)
process.p = cms.Path(process.simHitAnalyzer)

