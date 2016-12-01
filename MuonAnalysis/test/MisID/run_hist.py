#!/usr/bin/env python

import sys, os
from ROOT import *

gSystem.CompileMacro("MisID.C", "k")

proof = TProof.Open("")

modes = ["ks", "phi", "lamb"]
fileLists = [
    "JetHT_Run2016B-23Sep2016-v3.txt",
    "JetHT_Run2016C-23Sep2016-v1.txt",
    "JetHT_Run2016E-23Sep2016-v1.txt",
    "JetHT_Run2016F-23Sep2016-v1.txt",
    "JetHT_Run2016G-23Sep2016-v1.txt",
    "JetHT_Run2016H-PromptReco-v2.txt",
    "JetHT_Run2016H-PromptReco-v3.txt",
]

for mode in modes:
    chain = TChain("%s/tree" % mode)

    for fileList in fileLists:
        for f in open(fileList).readlines():
            f = f.strip()
            if not f.endswith(".root"): continue
            chain.Add(f)

    chain.SetProof()
    chain.Process("MisID.C+", "mode=%s" % mode)

