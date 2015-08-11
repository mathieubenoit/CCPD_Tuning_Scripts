						#
						#
						#
		#################################################################################
		# Call script: python -i threshold_CCPDv2.py DATAFILENAME ROOTFILENAME.root		#
		#																				#
		# 	Script to plot the threshold distribution from								#
		#	threshold scan using python functions 										#
		#################################################################################
						#
						#
						#
						#
# --------packages 

import sys
sys.path.insert(0, "/Users/misaelcaloz/Documents/GitHub/FunctionsRepository")
#sys.path.insert(0, "/Users/misaelcaloz/Documents/CERN master thesis/python scripts/functions")
#sys.path.append("/usr/local/lib/root")
#import ROOT
from ROOT import gROOT, TCanvas, TF1, TGraph, TLegend, TMath, TMultiGraph, TFile, TH2F, TH1F, TDirectory,TH2D, TH1D, TLatex, gPad, TH2I, TLine, TAxis, kTRUE, gStyle, TKey, THStack, kRed, kOrange, kBlue, kGreen
#from ROOT import *
import ROOT
from array import array
import numpy
import math
import pylab
from functions import *

# -------- Define sys.argv

FileDataThrScan = sys.argv[1]
rootFileName=sys.argv[2]

# -------- Create directories

outFile=TFile(rootFileName,"RECREATE")
GraphDir = outFile.mkdir("Graphs")
HighThrDir = outFile.mkdir("High thresh")
HighSigmaDir = outFile.mkdir("High sigma")
HighChi2Dir = outFile.mkdir("High Chi2")
# -------- Reset dictionaries and lists

TDAC_value_dic = {}
Threshold_value_dic = {}
Sigma_value_dic = {}
Chi2_value_dic = {}
Scurve_plot_dic = {}
Data_pointsX_dic = {}
Data_pointsY_dic = {}
DIC_AnalyseThresholdScan = {}

graphList = []

#-------- conversion factor injection -> electrons

injToElectrons=1660./.39	# WARNING: different depending the version of the chip: v2: 1660/0.39 v4:1660/0.25

# -------- Fill dictionaries with function AnalyseThresholdScan() in functions.pyc

DIC_AnalyseThresholdScan = AnalyseThresholdScan(FileDataThrScan)
TDAC_value_dic = DIC_AnalyseThresholdScan[0]
Threshold_value_dic = DIC_AnalyseThresholdScan[1]
Sigma_value_dic = DIC_AnalyseThresholdScan[2]
Chi2_value_dic = DIC_AnalyseThresholdScan[3]
Scurve_plot_dic = DIC_AnalyseThresholdScan[4]
Data_pointsX_dic = DIC_AnalyseThresholdScan[5]
Data_pointsY_dic = DIC_AnalyseThresholdScan[6]

# -------- Set histograms and fitting function

myfit = TF1("myfit", "[0]+0.5*TMath::Erf([1]*([2]+x))",-1.,3000.)
TDAC2D = TH2F("TDAC2D", "TDAC 2D plot;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
TDAC1D = TH1F("TDACDist", "TDAC distribution;Threshold [e];nb of pixels", 16,0, 16)
thresh2D = TH2F("thresh2D", "Thresh2D;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
thresh1D = TH1F("threshDist", "Threshold distribution;Threshold [e];nb of pixels", 300,-.1, 4000)
sigma2D = TH2F("sigma2D", "Sigma 2D plot;#Row;#Column;Sigma [e]",24,0,24, 36,12,48)
sigma1D = TH1F("sigmaDist", "Sigma distribution;Sigma [e];nb of pixels", 300,-.1, 500)
chi2_1D = TH1F("chi2Dist", "Chi2 distribution; ;nb of pixels", 300,0, 1)
chi2_2D = TH2F("chi2_2D", "Chi2 2D plot;#Row;#Column;chi2",24,0,24, 36,12,48)
correlation = TH2F("correlation", "Threshold vs TDAC;TDAC;Threshold [e]",16,0,16, 50,0,2000)
eyeDiagram = TH2D("eyeDiagram","S-curves eye diagram;injection [e];Probe",69*2/3,0,0.69*injToElectrons*2/3,128,0,1.5) 

#'$\ detected / simulated$'",69*2/3,0,0.69*injToElectrons*2/3,128,0,1.5)

# ------- Create Scurves for each pixel

outFile.cd("Graphs")

CanvScurves = TCanvas("S-curves_all")
CanvScurves.SetGrid()
CanvScurves.SetFillColor(0)
CanvScurves.cd()

cnt = 0
for r in range(24):
	for c in range(12,48):
		graph = Scurve_plot_dic["r"+str(r)+"_c"+str(c)+""]
		Xaxis = graph.GetXaxis()
		Xaxis.SetLimits(0,2000)
		graph.GetHistogram().SetMaximum(2)          
		graph.GetHistogram().SetMinimum(0)
		outFile.cd("Graphs")
		graph.Write()
		
		# select data
		
		if Threshold_value_dic["r"+str(r)+"_c"+str(c)+""] > 800:
			outFile.cd("High thresh")
			graph.Write()
		if Sigma_value_dic["r"+str(r)+"_c"+str(c)+""] > 100:
			outFile.cd("High sigma")
			graph.Write()			
		if Chi2_value_dic["r"+str(r)+"_c"+str(c)+""] > 0.01:
			outFile.cd("High Chi2")
			graph.Write()			
		
		
		if cnt ==0:	
			graph.Draw("AC")
		else:
			graph.Draw("C")		# method with Draw("same") also works
		
		cnt += 1
		graphList.append(graph)
gPad.Update()

# ------- Create histograms

for r in range(24):
	for c in range(12,48):
		TDAC_value = TDAC_value_dic["r"+str(r)+"_c"+str(c)+""]
		Threshold_value = Threshold_value_dic["r"+str(r)+"_c"+str(c)+""]
		Sigma_value = Sigma_value_dic["r"+str(r)+"_c"+str(c)+""]
		Chi2_value = Chi2_value_dic["r"+str(r)+"_c"+str(c)+""]
		Data_pointsX_list = Data_pointsX_dic["r"+str(r)+"_c"+str(c)+""]
		Data_pointsY_list = Data_pointsY_dic["r"+str(r)+"_c"+str(c)+""]
		
		# ---- Fill eyeDiagram 
		
		for i in range(len(Data_pointsX_list)):
			Data_pointsX_value = Data_pointsX_list[i]
			Data_pointsY_value = Data_pointsY_list[i]
			
			eyeDiagram.Fill(Data_pointsX_value,Data_pointsY_value)	
		
		# ----- Fill other histograms
		TDAC1D.Fill(TDAC_value)
		TDAC2D.Fill(r,c,TDAC_value)
		thresh1D.Fill(Threshold_value)
		thresh2D.Fill(r,c,Threshold_value)
		sigma1D.Fill(Sigma_value)
		sigma2D.Fill(r,c,Sigma_value)
		chi2_1D.Fill(Chi2_value)
		chi2_2D.Fill(r,c,Chi2_value)
		correlation.Fill(TDAC_value,Threshold_value)

outFile.cd()

# ------- Customize plots

# - eyeDiagram

eyeDiagram.SetAxisRange(0.,2000.,"X")
eyeDiagram.SetAxisRange(0.,1.4,"Y")

# - thresh 2D

thresh2D.SetTitleSize(0.025,"xyz")
thresh2D.SetTitleOffset(1.3,"z")

thresh2D.SetAxisRange(-1,1500.,"Z")
thresh2D.GetXaxis().SetNdivisions(32)
thresh2D.GetYaxis().SetNdivisions(64)
thresh2D.SetLabelSize(0.02,"X")
thresh2D.SetLabelSize(0.02,"Y")
thresh2D.GetZaxis().SetLabelSize(0.025)

x_axis_thresh2D = thresh2D.GetXaxis()
x_axis_thresh2D.CenterLabels(kTRUE)
y_axis_thresh2D = thresh2D.GetYaxis()
y_axis_thresh2D.CenterLabels(kTRUE)

# - thresh 1D

thresh1D.SetAxisRange(-1.,2000.,"X")
thresh1D.SetFillColor(38)
thresh1D.SetFillColorAlpha(38, 0.70)

# - Correlation

correlation.GetXaxis().SetNdivisions(16)
x_axis_correlation = correlation.GetXaxis()
x_axis_correlation.CenterLabels(kTRUE)

# - TDAC 2D

TDAC2D.GetZaxis().SetRangeUser(-1, 16)

# - sigma 1D

sigma1D.Rebin()

# ------- Create canvas

Canv4Plots = ROOT.TCanvas("Thresh_2D","Thresh_2D",800,600)
ROOT.SetOwnership(Canv4Plots,False) # TODO: DOESN'T WORK...

Canv4Plots.Divide(2,2)
Canv4Plots.cd(1)
gPad.SetGrid()
gStyle.SetOptStat("rmen")
thresh1D.Draw()
gPad.Update()


Canv4Plots.cd(2)
gStyle.SetOptStat("e")
gStyle.SetStatX(0.9)
gStyle.SetStatY(0.95)

thresh2D.GetZaxis().SetTitleOffset(1.3)
gPad.SetRightMargin(0.15)

thresh2D.Draw("colz")
TDAC2D.Draw("sametext")
TDAC2D.SetMarkerSize(0.9)
list_TLine=[]
for i in range(0,25):
	line_vert =  TLine(i,12,i,48)
	line_vert.Draw()
	list_TLine.append(line_vert)

for i in range(12,49):
	line_hor =  TLine(0,i,24,i)
	line_hor.Draw()
	list_TLine.append(line_hor)
	
gPad.Update()

Canv4Plots.cd(3)
gPad.SetGrid()
gStyle.SetOptStat("men")
gStyle.SetStatX(0.9)
gStyle.SetStatY(0.9)
correlation.Draw("colztext")

gPad.Update()

Canv4Plots.cd(4)
gPad.SetGrid()
gStyle.SetOptStat("n")
eyeDiagram.Draw("colz")

gPad.Update()

# ------- Create threshold 2D map 


Canv2DThresh = ROOT.TCanvas("2D_Thresh","Thresh_2D",800,600)
ROOT.SetOwnership(Canv4Plots,False) # TODO: DOESN'T WORK...

Canv2DThresh.cd()
gStyle.SetOptStat("e")
gStyle.SetStatX(0.9)
gStyle.SetStatY(0.95)

thresh2D.GetZaxis().SetTitleOffset(1.3)
gPad.SetRightMargin(0.15)

thresh2D.Draw("colz")
TDAC2D.Draw("sametext")
TDAC2D.SetMarkerSize(0.9)
list_TLine=[]
for i in range(0,25):
	line_vert =  TLine(i,12,i,48)
	line_vert.Draw()
	list_TLine.append(line_vert)

for i in range(12,49):
	line_hor =  TLine(0,i,24,i)
	line_hor.Draw()
	list_TLine.append(line_hor)
	
gPad.Update()

Canv2DThresh.Close()

# ------- Write plots 

outFile.cd()

thresh1D.Write()
thresh2D.Write()
TDAC1D.Write()
TDAC2D.Write()
correlation.Write()
sigma1D.Write()
sigma2D.Write()
chi2_1D.Write()
chi2_2D.Write()
eyeDiagram.Write()
CanvScurves.Write()
Canv4Plots.Write()
Canv2DThresh

#CanvScurves.Close()

#raw_input("Press Enter to continue...")


