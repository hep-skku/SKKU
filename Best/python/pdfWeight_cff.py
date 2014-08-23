import FWCore.ParameterSet.Config as cms

pdfWeight = cms.EDProducer("PDFWeightsProducer",
    pdfName = cms.string("cteq66.LHgrid"),
    altPdfNames = cms.vstring(["CT10.LHgrid","MSTW2008nlo68cl.LHgrid"]),#,"NNPDF23_nlo_as_0117.LHgrid"]),
)

