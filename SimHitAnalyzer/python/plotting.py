import ROOT

# Set ROOT to batch mode (no GUI or window display)
#ROOT.gROOT.SetBatch(True)

# Open the ROOT file and navigate to the directory
f = ROOT.TFile.Open("simHits.root")
f.cd("simHitAnalyzer")

# Retrieve the histogram
h = f.Get("hXY")

# Create a canvas in batch mode
c = ROOT.TCanvas("c", "c", 800, 600)

# Draw and save with minimal overhead
h.Draw("COLZ")  # Required to render the histogram into canvas before saving
c.SaveAs("hXY.pdf")

