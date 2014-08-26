#!/usr/bin/env python

from ROOT import *

def analyze(fileName):
    f = TFile(fileName)
    #h = f.Get("SEvt/hnvertex")
    h = f.Get("SEvt/hMt")
    h.Sumw2()
    h.Scale(1./h.Integral())
    #if 'central_non' in fileName:
    #    h.Scale(1e-8) # Scale by lumi/nEvents

    #h = f.Get("SEvt/hnjets")
    result = []
    for i in range(1, h.GetNbinsX()+1):
        x  = h.GetBinCenter(i)
        xLo = h.GetBinLowEdge(i)
        xHi = xLo + h.GetBinWidth(i)
        y  = h.GetBinContent(i)
        ey = h.GetBinError(i)
        result.append((x, xLo, xHi, y, ey, ey))
    return result
