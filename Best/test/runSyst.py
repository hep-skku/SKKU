#!/usr/bin/env python
import sys, os

from ROOT import *
if os.path.exists("rootlogon.C"):
    gROOT.ProcessLine(".x rootlogon.C")

from array import array
from math import *

sys.path += "."
from analysis import analyze

histPath = "hist"

central = "hist_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root"

systs = {
    "bJES-":("hist_JESDn_MSDecays_central_non_FlavorPureBottom_pt35_nj4_nb2_lq0_nlj3.root",
             "hist_JESDn_MSDecays_central_non_FlavorPureCharm_pt35_nj4_nb2_lq0_nlj3.root",
             "hist_JESDn_MSDecays_central_non_FlavorPureGluon_pt35_nj4_nb2_lq0_nlj3.root",
             "hist_JESDn_MSDecays_central_non_FlavorPureQuark_pt35_nj4_nb2_lq0_nlj3.root",),
    "bJES+":("hist_JESUp_MSDecays_central_non_FlavorPureBottom_pt35_nj4_nb2_lq0_nlj3.root",
             "hist_JESUp_MSDecays_central_non_FlavorPureCharm_pt35_nj4_nb2_lq0_nlj3.root",
             "hist_JESUp_MSDecays_central_non_FlavorPureGluon_pt35_nj4_nb2_lq0_nlj3.root",
             "hist_JESUp_MSDecays_central_non_FlavorPureQuark_pt35_nj4_nb2_lq0_nlj3.root",),

    "JER-":("hist_JERDn_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root",),
    "JER+":("hist_JERUp_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root",),

    "JES-":("hist_JESDn_MSDecays_central_non_CorrelationGroupIntercalibration_pt35_nj4_nb2_lq0_nlj3.root",
            "hist_JESDn_MSDecays_central_non_CorrelationGroupMPFInSitu_pt35_nj4_nb2_lq0_nlj3.root",
            "hist_JESDn_MSDecays_central_non_CorrelationGroupUncorrelated_pt35_nj4_nb2_lq0_nlj3.root",),
    "JES+":("hist_JESUp_MSDecays_central_non_CorrelationGroupIntercalibration_pt35_nj4_nb2_lq0_nlj3.root",
            "hist_JESUp_MSDecays_central_non_CorrelationGroupMPFInSitu_pt35_nj4_nb2_lq0_nlj3.root",
            "hist_JESUp_MSDecays_central_non_CorrelationGroupUncorrelated_pt35_nj4_nb2_lq0_nlj3.root",),

    "matching-":("hist_MSDecays_matchingdown_non_v2_pt35_nj4_nb2_lq0_nlj3.root",),
    "matching+":("hist_MSDecays_matchingup_non_pt35_nj4_nb2_lq0_nlj3.root",),

    "MCgenerator":("hist_powheg_tauola_non_pt35_nj4_nb2_lq0_nlj3.root",),

    "Pileup-":("hist_PuDn_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root",),
    "Pileup+":("hist_PuUp_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root",),

    "PileupJES-":("hist_JESDn_MSDecays_central_non_PileUpPtBB_pt35_nj4_nb2_lq0_nlj3.root",
                  "hist_JESDn_MSDecays_central_non_PileUpPtEC1_pt35_nj4_nb2_lq0_nlj3.root",
                  "hist_JESDn_MSDecays_central_non_PileUpPtEC2_pt35_nj4_nb2_lq0_nlj3.root",),
    "PileupJES+":("hist_JESUp_MSDecays_central_non_PileUpPtBB_pt35_nj4_nb2_lq0_nlj3.root",
                  "hist_JESUp_MSDecays_central_non_PileUpPtEC1_pt35_nj4_nb2_lq0_nlj3.root",
                  "hist_JESUp_MSDecays_central_non_PileUpPtEC2_pt35_nj4_nb2_lq0_nlj3.root",),

    "PtWeight":("hist_Ptweight_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root",),

    "Qscale-":("hist_MSDecays_scaledown_non_pt35_nj4_nb2_lq0_nlj3.root",),
    "Qscale+":("hist_MSDecays_scaleup_non_pt35_nj4_nb2_lq0_nlj3.root",),

    "SFb-":("hist_SFbDn_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root",),
    "SFb+":("hist_SFbUp_MSDecays_central_non_pt35_nj4_nb2_lq0_nlj3.root",),

##### TODO: Tune is difficult
#    "Tune+":("hist_TuneP11mpiUp_non_pt35_nj4_nb2_lq0_nlj3.root"),
#    "Tune-":("hist_TuneP11TeV_non_pt35_nj4_nb2_lq0_nlj3.root",),
#            "hist_TuneP11_non_pt35_nj4_nb2_lq0_nlj3.root",
#            "hist_TuneP11noCR_non_pt35_nj4_nb2_lq0_nlj3.root",
}

combineIgnore = []
combineByScalarSum = ["bJES"]
combineByMaximum = []

# Perform analysis for each files and get results
# Results should be in the format of [(x, xDn, xUp, y, yErrDn, yErrUp), (), (), ...]
# So, if the result is a single measurement, one can set the value for example, [(0, 0, 0, 172.5, 0.2, 0.3)]
result0 = analyze(os.path.join(histPath, central))
errors = {}
for syst in systs:
    # Store combined systematic uncertainty for this uncertainty group
    errsUp = [0]*len(result0)
    errsDn = [0]*len(result0)

    if   syst in combineIgnore: combineBy = 'ignore'
    elif syst in combineByScalarSum: combineBy = 'scalarSum'
    elif syst in combineByMaximum: combineBy = 'maximum'
    else: combineBy = 'squareSum'

    # Dnop over all systematic variations in this group
    for fileName in systs[syst]:
        result = analyze(os.path.join(histPath, fileName))
        if len(result) != len(result0):
            print "!!! FATAL: length of result is not same !!!"
            os.exit(1)

        for iBin, (x, xDn, xUp, y, eyDn, eyUp) in enumerate(result):
            y0, eyDn0, eyUp0 = result0[iBin][3:]
            dy = y-y0
            # Choose stat error of this measurement if stat. error is large
            if   dy < 0 and abs(dy) < eyDn: dy = -abs(eyDn)
            elif dy > 0 and abs(dy) < eyUp: dy = eyUp
            elif dy == 0: # Conservative choice if difference is 0
                if abs(eyDn) < eyUp: dy = eyUp
                else: dy = -abs(eyDn)

            # Combine with previous one depending on combine mode
            if   combineBy == 'ignore': pass
            elif combineBy == 'scalarSum':
                if dy < 0: errsDn[iBin] += abs(dy)
                else: errsUp[iBin] += abs(dy)
            elif combineBy == 'maximum':
                if dy < 0: errsDn[iBin] = max(errsDn[iBin], abs(dy))
                else: errsUp[iBin] = max(errsUp[iBin], abs(dy))
            elif combineBy == 'squareSum':
                if dy < 0: errsDn[iBin] += dy**2
                else: errsUp[iBin] += dy**2
            else:
                print "!!! FATAL: combine method is wrong"
                os.exit(1)
    # Take square root if this is square-sum mode
    if combineBy == 'squareSum':
        errsUp = [sqrt(x) for x in errsUp]
        errsDn = [sqrt(x) for x in errsDn]

    errors[syst] = (errsUp, errsDn)

## Re-Combine total error
errsDn = [0]*len(result0)
errsUp = [0]*len(result0)
for syst in sorted(errors.keys()):
    if syst not in errors: continue
    if len(errors) == 0: break
    isAsymm = False
    if syst[-1] in ('+', '-') and \
       syst[:-1]+'+' in systs and syst[:-1]+'-' in systs:
        isAsymm = True
        errsAUp, errsADn = errors.pop(syst[:-1]+'+')
        errsBUp, errsBDn = errors.pop(syst[:-1]+'-')
    else:
        errsAUp, errsADn = errors.pop(syst)
        errsBUp, errsBDn = errsAUp, errsADn

    for iBin in range(len(errsAUp)):
        errUp = max(errsAUp[iBin], errsBUp[iBin])
        errDn = max(errsADn[iBin], errsBDn[iBin])
        # Use symmetric error if the error is only one side
        if not isAsymm:
            errUp = errDn = max(errUp, errDn)
        errsUp[iBin] += errUp**2
        errsDn[iBin] += errDn**2

errsUp = [sqrt(x) for x in errsUp]
errsDn = [sqrt(x) for x in errsDn]

## Produce summary plot
xbins = array('d', [x[1] for x in result0])
if len(xbins) != result0: xbins.append(result0[-1][1])
h = TH1F("h", "h", len(xbins)-1, xbins)
grpSyst  = TGraphAsymmErrors()
grpTotal = TGraphAsymmErrors()
for iBin in range(len(result0)):
    x, xDn, xUp, y, eyDn0, eyUp0 = result0[iBin]
    eyDn = errsDn[iBin]
    eyUp = errsUp[iBin]
    #print "%d: %f + %f - %f (stat.) +%f -%f (syst.)" % (iBin, y, eyUp0, eyDn0, eyUp, eyDn)
    h.SetBinContent(iBin+1, y)
    h.SetBinError(iBin+1, max(eyDn0, eyUp0))
    grpSyst.SetPoint(iBin, x, y)
    grpSyst.SetPointError(iBin, x-xDn, xUp-x, eyDn, eyUp)
    grpTotal.SetPoint(iBin, x, y)
    grpTotal.SetPointError(iBin, x-xDn, xUp-x, eyDn0+eyDn, eyUp0+eyUp)

def drawRatioPlot(h, grpTotal, grpSyst=None):
    c = TCanvas("c", "c", 500, 600)
    c.Divide(1,2)

    pad1 = c.cd(1)
    pad1.SetPad(0, 0.3, 1, 1)
    pad1.SetBottomMargin(0)
    pad1.SetTopMargin(0.12)
    pad1.SetLeftMargin(0.2)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(1.70)
    h.GetYaxis().SetLabelSize(0.05)
    h.Draw()
    if grpSyst != None:
        grpSyst.SetFillColor(kYellow+1)
        grpSyst.SetFillStyle(3001)
        grpSyst.Draw("pz2")
    grpTotal.Draw("pz")

    pad2 = c.cd(2)
    pad2.SetPad(0, 0, 1, 0.3)
    pad2.SetBottomMargin(0.3)
    pad2.SetTopMargin(0)
    pad2.SetLeftMargin(0.2)
    hRatio = h.Clone()
    hRatio.Reset()
    hRatio.SetName("hRatio_%s" % h.GetName())
    hRatio.GetYaxis().SetTitle("Data/MC")
    hRatio.GetYaxis().SetTitleSize(0.12)
    hRatio.GetYaxis().SetTitleOffset(0.75)
    hRatio.GetYaxis().SetLabelSize(0.12)
    hRatio.GetYaxis().SetNdivisions(505)
    hRatio.GetXaxis().SetTitleSize(0.15)
    hRatio.GetXaxis().SetTitleOffset(0.90)
    hRatio.GetXaxis().SetLabelSize(0.12)
    hRatio.SetMinimum(0.5)
    hRatio.SetMaximum(1.5)

    grpRatioTotal = TGraphAsymmErrors()
    if grpSyst != None: grpRatioSyst = TGraphAsymmErrors()
    for i in range(grpTotal.GetN()):
        y0 = h.GetBinContent(i+1)
        x = grpTotal.GetX()[i]
        y = grpTotal.GetY()[i]
        if y > 0:
            r = y0/y
            erDn = r*hypot(grpTotal.GetEYlow()[i]/y, h.GetBinError(i+1)/y0)
            erUp = r*hypot(grpTotal.GetEYhigh()[i]/y, h.GetBinError(i+1)/y0)
        else:
            r, erDn, erUp = 0, 0, 1e9
        exDn = grpTotal.GetEXlow()[i]
        exUp = grpTotal.GetEXhigh()[i]

        grpRatioTotal.SetPoint(i, x, r)
        grpRatioTotal.SetPointError(i, exDn, exUp, erDn, erUp)

    grpRatioSyst = TGraphAsymmErrors()
    grpRatioSyst.SetFillColor(kYellow+1)
    grpRatioSyst.SetFillStyle(3001)
    if grpSyst != None:
        for i in range(grpSyst.GetN()):
            y0 = h.GetBinContent(i+1)
            x = grpSyst.GetX()[i]
            y = grpSyst.GetY()[i]
            if y > 0:
                r = y0/y
                erDn = r*grpSyst.GetEYlow()[i]/y
                erUp = r*grpSyst.GetEYhigh()[i]/y
            else:
                r, erDn, erUp = 0, 0, 1e9
            exDn = grpSyst.GetEXlow()[i]
            exUp = grpSyst.GetEXhigh()[i]

            grpRatioSyst.SetPoint(i, x, r)
            grpRatioSyst.SetPointError(i, exDn, exUp, erDn, erUp)

    hRatio.Draw()
    grpRatioTotal.Draw("pz0")
    if grpRatioSyst != None: grpRatioSyst.Draw("pz02")

    c.Update()

    grpRatioTotal.SetEditable(False)
    grpRatioSyst.SetEditable(False)
    return [c, hRatio, grpRatioTotal, grpRatioSyst]

grpTotal.SetEditable(False)
grpSyst.SetEditable(False)
c = drawRatioPlot(h, grpTotal, grpSyst)
