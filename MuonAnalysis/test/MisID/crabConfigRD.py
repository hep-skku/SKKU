from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'prod_RD_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 200
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
config.Data.lumiMask = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T3_KR_KISTI'
#crab checkwrite --site=T3_KR_KISTI --lfn=/store/group/CAT/

import os
era = os.environ["ERA"]

pdName = "JetHT"
options = {
    "2016B":("Run2016B-23Sep2016-v3",[]),
    "2016C":("Run2016C-23Sep2016-v1",[]),
    "2016D":("Run2016D-23Sep2016-v1",[]),
    "2016E":("Run2016E-23Sep2016-v1",[]),
    "2016F":("Run2016F-23Sep2016-v1",[]),
    "2016G":("Run2016G-23Sep2016-v1",[]),
    "2016Hv2":("Run2016H-PromptReco-v2",['globalTag=80X_dataRun2_Prompt_v14']),
    "2016Hv3":("Run2016H-PromptReco-v3",['globalTag=80X_dataRun2_Prompt_v14'])
}

sdName, params = options[era]

config.General.requestName = "%s_%s" % (pdName, era)
config.General.workArea = "MuonMisID"
config.Data.inputDataset = '/%s/%s/AOD' % (pdName, sdName)
config.Data.outputDatasetTag = '%s_%s' % (pdName, sdName)
config.Data.outLFNDirBase = '/store/user/jhgoh/MuonMisID/20161125_1'
config.JobType.pyCfgParams = params

