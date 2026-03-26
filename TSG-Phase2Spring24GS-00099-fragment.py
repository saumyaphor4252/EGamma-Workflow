import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8ConcurrentGeneratorFilter",
	comEnergy = cms.double(14000.0),
	maxEventsToPrint = cms.untracked.int32(1),
	pythiaHepMCVerbosity = cms.untracked.bool(False),
	pythiaPylistVerbosity = cms.untracked.int32(1),
	PythiaParameters = cms.PSet(        
		pythia8CommonSettingsBlock,
		pythia8CP5SettingsBlock,
		processParameters = cms.vstring(
			'NewGaugeBoson:ffbar2gmZZprime = on',
			'Zprime:gmZmode = 3', # only pure Z' contribution
			'32:m0 = 6000', # mass of Z'
			'32:onMode = off', # switch off all of the Z' decay
			'32:onIfAny = 11', # switch on the Z'->tau-tau+
		),
		parameterSets = cms.vstring(
			'pythia8CommonSettings',
			'pythia8CP5Settings',
			'processParameters'
        )       
	)
)

ProductionFilterSequence = cms.Sequence(generator)
