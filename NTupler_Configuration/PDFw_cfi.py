import FWCore.ParameterSet.Config as cms

PDFPSets = cms.PSet(

        PDFWcteq66_0 = cms.PSet(
               src = cms.InputTag("pdfWeights","cteq66"),
               method = cms.string("DoubleVVar"), #get a vector of double from event
               index = cms.uint32(0)
               ),
    
        # eigenvector sets
        PDFWcteq66_1 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(1)
            ),

        PDFWcteq66_2 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(2)    
            ),

        PDFWcteq66_3 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(3)
            ),

        PDFWcteq66_4 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(4)
            ),

        PDFWcteq66_5 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(5)
            ),

        PDFWcteq66_6 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(6)
            ),

        PDFWcteq66_7 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(7)
            ),

        PDFWcteq66_8 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(8)
            ),

        PDFWcteq66_9 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(9)
            ),

        PDFWcteq66_10 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(10)
            ),

        PDFWcteq66_11 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(11)
            ),

        PDFWcteq66_12 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(12)
            ),

        PDFWcteq66_13 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(13)
            ),

        PDFWcteq66_14 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(14)
            ),

        PDFWcteq66_15 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(15)
            ),

        PDFWcteq66_16 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(16)
            ),

        PDFWcteq66_17 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(17)
            ),

        PDFWcteq66_18 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(18)
            ),

        PDFWcteq66_19 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(19)
            ),

        PDFWcteq66_20 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(20)
            ),

        PDFWcteq66_21 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(21)
            ),

        PDFWcteq66_22 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(22)
            ),

        PDFWcteq66_23 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(23)
            ),

        PDFWcteq66_24 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(24)
            ),

        PDFWcteq66_25 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(25)
            ),

        PDFWcteq66_26 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(26)
            ),

        PDFWcteq66_27 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(27)
            ),

        PDFWcteq66_28 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(28)
            ),

        PDFWcteq66_29 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(29)
            ),

        PDFWcteq66_30 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(30)
            ),

        PDFWcteq66_31 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(31)
            ),

        PDFWcteq66_32 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(32)
            ),

        PDFWcteq66_33 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(33)
            ),

        PDFWcteq66_34 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(34)
            ),

        PDFWcteq66_35 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(35)
            ),

        PDFWcteq66_36 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(36)
            ),

        PDFWcteq66_37 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(37)
            ),

        PDFWcteq66_38 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(38)
            ),

        PDFWcteq66_39 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(39)
            ),

        PDFWcteq66_40 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(40)
            ),

        PDFWcteq66_41 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(41)
            ),

        PDFWcteq66_42 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(42)
            ),

        PDFWcteq66_43 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(43)
            ),

        PDFWcteq66_44 = cms.PSet(
            src = cms.InputTag("pdfWeights","cteq66"),
            method = cms.string("DoubleVVar"),
            index = cms.uint32(44)
            )
)
