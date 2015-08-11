						#
						#
						#
		#################################################################################
		# Call script: python -i threshold_CCPDv2_Tuning_usingFunctions.py				#
		#	OUTPUTFILE.root dataTDAC.txt data2TDAC.txt .........						#
		#																				#
		# 	Script to plot the threshold distribution from								#
		#	threshold scan																#
		#################################################################################
						#
						#
						#

#=======================================
# --------packages 

import sys
sys.path.insert(0, "/Users/misaelcaloz/Documents/GitHub/FunctionsRepository")
#sys.path.insert(0, "/Users/misaelcaloz/Documents/CERN master thesis/python scripts/functions")
sys.path.append("/usr/local/lib/root")
#import ROOT
from ROOT import gROOT, TCanvas, TF1, TGraph, TLegend, TMath, TMultiGraph, TFile, TH2F, TH1F, TDirectory,TH2D, TH1D, TLatex, gPad, TH2I, TLine, TAxis, kTRUE, gStyle, TKey, THStack, kRed, kOrange, kBlue, kGreen
from array import array
import numpy
import math
import pylab
import os
from decimal import Decimal
from functions import *

#------------ input data
rootFileName=sys.argv[1]

datafileTDAC2 = sys.argv[2]
datafileTDAC4 = sys.argv[3]
datafileTDAC6 = sys.argv[4]
datafileTDAC9 = sys.argv[5]
datafileTDAC11 = sys.argv[6]
#datafileTDAC14 = sys.argv[7]

dataTDAC2 = open(datafileTDAC2)
dataTDAC4 = open(datafileTDAC4)
dataTDAC6 = open(datafileTDAC6)
dataTDAC9 = open(datafileTDAC9)
dataTDAC11 = open(datafileTDAC11)
#dataTDAC14 = open(datafileTDAC14)

#dataTDAC_dic = {"dataTDAC4": 2,
#				"dataTDAC4": 4,
#                "dataTDAC6": 6,
#                "dataTDAC6": 10,
#                "dataTDAC6": 14,
#                "dataTDAC6": 14
#                }
			
FileTDAC_list = [sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]] #,sys.argv[7]]
TDAC_list = [2,4,6,10,14]

#-------- conversion factor injection -> electrons

injToElectrons=1660./.39			# WARNING: different depending the version of the chip: v2: 1660/0.39 v4: 1660/0.25

# ----------- Target threshold

target_threshold = 700.

#------------ create file.root

outFile=TFile(rootFileName,"RECREATE")

TDAC2Dir = outFile.mkdir("TDAC 2")
graphDirTDAC2 = TDAC2Dir.mkdir("Graphs")
FailFitDirTDAC2 = TDAC2Dir.mkdir("Fail_fit")
highThreshDirTDAC2 = TDAC2Dir.mkdir("High_thresh")
highSigmaDirTDAC2 = TDAC2Dir.mkdir("High_sigma")
PlotsTDAC2Dir = TDAC2Dir.mkdir("Plots")

TDAC4Dir = outFile.mkdir("TDAC 4")
graphDirTDAC4 = TDAC4Dir.mkdir("Graphs")
FailFitDirTDAC4 = TDAC4Dir.mkdir("Fail_fit")
highThreshDirTDAC4 = TDAC4Dir.mkdir("High_thresh")
highSigmaDirTDAC4 = TDAC4Dir.mkdir("High_sigma")
PlotsTDAC4Dir = TDAC4Dir.mkdir("Plots")

TDAC6Dir = outFile.mkdir("TDAC 6")
graphDirTDAC6 = TDAC6Dir.mkdir("Graphs")
FailFitDirTDAC6 = TDAC6Dir.mkdir("Fail_fit")
highThreshDirTDAC6 = TDAC6Dir.mkdir("High_thresh")
highSigmaDirTDAC6 = TDAC6Dir.mkdir("High_sigma")
PlotsTDAC6Dir = TDAC6Dir.mkdir("Plots")

TDAC9Dir = outFile.mkdir("TDAC 10")
graphDirTDAC9 = TDAC9Dir.mkdir("Graphs")
FailFitDirTDAC9 = TDAC9Dir.mkdir("Fail_fit")
highThreshDirTDAC9 = TDAC9Dir.mkdir("High_thresh")
highSigmaDirTDAC9 = TDAC9Dir.mkdir("High_sigma")
PlotsTDAC9Dir = TDAC9Dir.mkdir("Plots")

TDAC11Dir = outFile.mkdir("TDAC 14")
graphDirTDAC11 = TDAC11Dir.mkdir("Graphs")
FailFitDirTDAC11 = TDAC11Dir.mkdir("Fail_fit")
highThreshDirTDAC11 = TDAC11Dir.mkdir("High_thresh")
highSigmaDirTDAC11 = TDAC11Dir.mkdir("High_sigma")
PlotsTDAC11Dir = TDAC11Dir.mkdir("Plots")

#TDAC14Dir = outFile.mkdir("TDAC 14")
#graphDirTDAC14 = TDAC14Dir.mkdir("Graphs")
#FailFitDirTDAC14 = TDAC14Dir.mkdir("Fail_fit")
#highThreshDirTDAC14 = TDAC14Dir.mkdir("High_thresh")
#highSigmaDirTDAC14 = TDAC14Dir.mkdir("High_sigma")
#PlotsTDAC14Dir = TDAC14Dir.mkdir("Plots")

PlotsDir = outFile.mkdir("TDACsuggestion")
PixelsDir = PlotsDir.mkdir("Pixels")
FailThrDir = PlotsDir.mkdir("Failing thresh")
HighThrDic = PlotsDir.mkdir("High threshold")
GoodPixDic = PlotsDir.mkdir("Good pixels")

Pix1Dic = PlotsDir.mkdir("0 point")
Pix1Dic = PlotsDir.mkdir("1 point")
Pix2Dic = PlotsDir.mkdir("2 points")
Pix3Dic = PlotsDir.mkdir("3 points")
Pix4Dic = PlotsDir.mkdir("4 points")
Pix5Dic = PlotsDir.mkdir("5 points")
Pix6Dic = PlotsDir.mkdir("6 points")

PixWarning6points = Pix6Dic.mkdir("Warning")
PixGood6points = Pix6Dic.mkdir("Good pixels")


AllThresh_dic = {}
AllSigma_dic = {}
AllChi2_dic = {}


for n, TDACcnt in zip(range(2,7), TDAC_list):
	FileDataThrScan = sys.argv[n]
		
	# -------- Reset dictionaries and lists

	TDAC_value_dic = {}
	Threshold_value_dic = {}
	Sigma_value_dic = {}
	Chi2_value_dic = {}
	Scurve_plot_dic = {}
	Data_pointsX_dic = {}
	Data_pointsY_dic = {}
	graphList = []
	DIC_AnalyseThresholdScan = {}
	
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
	thresh2D = TH2F("thresh2D", "Threshold 2D plot;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
	thresh1D = TH1F("threshDist", "Threshold distribution;Threshold [e];nb of pixels", 300,-.1, 4000)
	sigma2D = TH2F("sigma2D", "Sigma 2D plot;#Row;#Column;Sigma [e]",24,0,24, 36,12,48)
	sigma1D = TH1F("sigmaDist", "Sigma distribution;Sigma [e];nb of pixels", 300,-.1, 500)
	chi2_1D = TH1F("chi2Dist", "Chi2 distribution; ;nb of pixels", 300,0, 1)
	chi2_2D = TH2F("chi2_2D", "Chi2 2D plot;#Row;#Column;chi2",24,0,24, 36,12,48)
	correlation = TH2F("correlation", "Threshold vs TDAC;TDAC;Threshold [e]",16,0,16, 50,0,2000)
	eyeDiagram = TH2D("eyeDiagram","S-curves eye diagram;injection [e];Probe",69*2/3,0,0.69*injToElectrons*2/3,128,0,1.5) 
		
	
	
	if n == 2:
		TDAC2Dir.cd("Graphs")
	
	elif n == 3:
		TDAC4Dir.cd("Graphs")
	
	elif n == 4:
		TDAC6Dir.cd("Graphs")
	
	elif n == 5:
		TDAC9Dir.cd("Graphs")
	
	elif n == 6:
		TDAC11Dir.cd("Graphs")

#	elif n == 7:
#		TDAC14Dir.cd("Graphs")
	
	else: 
		print 'ERROR'
		
		
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
			graph.GetHistogram().SetMaximum(1.2)          
			graph.GetHistogram().SetMinimum(0)
			graph.SetMarkerStyle(3)
			graph.Draw("AlP")
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

			AllThresh_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Threshold_value
			AllSigma_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Sigma_value
			AllChi2_dic ["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] = Chi2_value



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

	# ------- Create plots in canvas

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

	# ------- Write plots 

	if n == 2:
		TDAC2Dir.cd()
	
	elif n == 3:
		TDAC4Dir.cd()
	
	elif n == 4:
		TDAC6Dir.cd()
	
	elif n == 5:
		TDAC9Dir.cd()
	
	elif n == 6:
		TDAC11Dir.cd()

#	elif n == 7:
#		TDAC14Dir.cd()
	
	else: 
		print 'ERROR'
		
		
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

	CanvScurves.Close()
	Canv4Plots.Close()
	
	thresh1D.Delete()
	thresh2D.Delete()
	TDAC1D.Delete()
	TDAC2D.Delete()
	correlation.Delete()
	sigma1D.Delete()
	sigma2D.Delete()
	chi2_1D.Delete()
	chi2_2D.Delete()
	eyeDiagram.Delete()
	CanvScurves.Delete()
	Canv4Plots.Delete()
	
	
	
# =====================================================================	

# -----------        Analysis and suggestions of TDAC for every pixels

# =====================================================================



# Define function to set TDAC value and create a .txt file 


file_TDAC = open("suggested_TDACfileTuningInterpolation.txt", "w")  # "r" instead of "w" in order to take the original suggested_TDACfile as an input instead of create TDACfile with all 6

file_TDAC_list = []
lineNum_TDAC=0

for counter in range(1,1441): # old way: create a TDAC all 6 file and modify it
	if counter < 289:
		file_TDAC_list.append(15)
	elif 288 < counter < 1153:
		file_TDAC_list.append(6)
	elif counter > 1152:
		file_TDAC_list.append(15)

# Define function to set the TDAC

def setTDAC(CMOSrow,CMOScol,newTDAC_value):
	counter2=0
	for i in range(60) :
		for j in range(24):
			if j==CMOSrow and i==CMOScol:
				file_TDAC_list[counter2]=newTDAC_value
			counter2 += 1



# Initialize all dictionaries and lists just to be sure

AllThresh_dic = {}
AllSigma_dic = {}
AllChi2_dic = {}
AllData_pointsX_dic = {}
AllData_pointsY_dic = {}
AllSuggTDAC_dic = {}
AllExpectThresh_dic = {}
AllThr_vs_TDAC_graphs_dic = {}
AllThresh_simpl_dic = {}
DIC_AnalyseTDACThresholdScansV2 = {}

# Call function to suggest TDAC
	
DIC_AnalyseTDACThresholdScansV2 = AnalyseTDACThresholdScansV2(FileTDAC_list,TDAC_list,target_threshold)	

AllThresh_dic = DIC_AnalyseTDACThresholdScansV2[0]
AllSigma_dic = DIC_AnalyseTDACThresholdScansV2[1]
AllChi2_dic = DIC_AnalyseTDACThresholdScansV2[2]
AllData_pointsX_dic = DIC_AnalyseTDACThresholdScansV2[3]
AllData_pointsY_dic = DIC_AnalyseTDACThresholdScansV2[4]
AllSuggTDAC_dic = DIC_AnalyseTDACThresholdScansV2[5]
AllExpectThresh_dic = DIC_AnalyseTDACThresholdScansV2[6]
AllThr_vs_TDAC_graphs_dic = DIC_AnalyseTDACThresholdScansV2[7]
AllThresh_simpl_dic = DIC_AnalyseTDACThresholdScansV2[8]

# Plot graphics

PlotsDir.cd("Pixels")

Sugg_TDAC2D = TH2F("Sugg_TDAC2D", "Suggested TDAC 2D plot;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
Sugg_TDAC1D = TH1F("Sugg_TDACDist", "Suggested TDAC distribution;Threshold [e];nb of pixels", 16,0, 16)
Exp_thresh2D = TH2F("Exp_thresh2D", "Expected Threshold 2D plot;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)
Exp_thresh1D = TH1F("Exp_threshDist", "Expected Threshold distribution;Threshold [e];nb of pixels", 300,-.1, 2000)
points6Warning2D = TH2F("points6Warning", "Warning pixels location;#Row;#Column;Threshold [e]",24,0,24, 36,12,48)

for r in range(24):
	for c in range(12,48):
		PlotsDir.cd("Pixels")
		graph = AllThr_vs_TDAC_graphs_dic["r"+str(r)+"_c"+str(c)+""]
		Xaxis = graph.GetXaxis()
#		Xaxis.SetLimits(0,2000)
#		graph.GetHistogram().SetMaximum(1.2)          
#		graph.GetHistogram().SetMinimum(0)
		graph.SetMarkerStyle(3)
		graph.Draw("AlP")
		graph.Write()
		
		points6_list = []
		for TDACcnt in TDAC_list:
			if AllThresh_simpl_dic["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""] > 1:
#
				points6_list.append(AllThresh_simpl_dic["r"+str(r)+"_c"+str(c)+"_TDAC"+str(TDACcnt)+""])
#		
		outFile.cd("TDACsuggestion")
#		
#		if len(points6_list) == 6:
#			PlotsDir.cd("6 points")
#			graph.Write()
#			
#			if points6_list[2] > points6_list[3]:
#				Pix6Dic.cd("Warning")
#				graph.Write()
#				points6Warning2D.Fill(r,c)
#			else:
#				Pix6Dic.cd("Good pixels")
#				graph.Write()				
#				
			
		if len(points6_list) == 5:
			PlotsDir.cd("5 points")
			graph.Write()
		if len(points6_list) == 4:
			PlotsDir.cd("4 points")
			graph.Write()
		if len(points6_list) == 3:
			PlotsDir.cd("3 points")
			graph.Write()
		if len(points6_list) == 2:
			PlotsDir.cd("2 points")
			graph.Write()
		if len(points6_list) == 1:
			PlotsDir.cd("1 point")
			graph.Write()
		if len(points6_list) == 0:
			PlotsDir.cd("0 point")
			graph.Write()

		outFile.cd("TDACsuggestion")


outFile.cd("TDACsuggestion")

for r in range(24):
	for c in range(12,48):
		Exp_thresh_value = AllExpectThresh_dic["r"+str(r)+"_c"+str(c)+""]
		Sugg_TDAC_value = AllSuggTDAC_dic ["r"+str(r)+"_c"+str(c)+""]
		
		Exp_thresh1D.Fill(Exp_thresh_value)
		Sugg_TDAC1D.Fill(Sugg_TDAC_value)
		Exp_thresh2D.Fill(r,c,Exp_thresh_value)
		Sugg_TDAC2D.Fill(r,c,Sugg_TDAC_value)
		
		# Set TDAC value for each pixel
		
		setTDAC(r,c,Sugg_TDAC_value)
		
		
# --- Customize plots

# thresh 1D

Exp_thresh1D.SetAxisRange(-1.,2000.,"X")
Exp_thresh1D.SetFillColor(38)
Exp_thresh1D.SetFillColorAlpha(38, 0.70)

# threshd 2D

Exp_thresh2D.SetTitleSize(0.025,"xyz")
Exp_thresh2D.SetTitleOffset(1.3,"z")

Exp_thresh2D.SetAxisRange(-1,1500.,"Z")
Exp_thresh2D.GetXaxis().SetNdivisions(32)
Exp_thresh2D.GetYaxis().SetNdivisions(64)
Exp_thresh2D.SetLabelSize(0.02,"X")
Exp_thresh2D.SetLabelSize(0.02,"Y")
Exp_thresh2D.GetZaxis().SetLabelSize(0.025)

Exp_x_axis_thresh2D = Exp_thresh2D.GetXaxis()
Exp_x_axis_thresh2D.CenterLabels(kTRUE)
Exp_y_axis_thresh2D = Exp_thresh2D.GetYaxis()
Exp_y_axis_thresh2D.CenterLabels(kTRUE)
	
# TDAC 2D
	
Sugg_TDAC2D.GetZaxis().SetRangeUser(-1, 16)

# TDAC 1D

Sugg_TDAC1D.SetAxisRange(-1,16,"X")

# Warning 6 points
		
points6Warning2D.SetTitleSize(0.025,"xyz")
points6Warning2D.SetTitleOffset(1.3,"z")

points6Warning2D.GetXaxis().SetNdivisions(32)
points6Warning2D.GetYaxis().SetNdivisions(64)
points6Warning2D.SetLabelSize(0.02,"X")
points6Warning2D.SetLabelSize(0.02,"Y")
points6Warning2D.GetZaxis().SetLabelSize(0.025)

Exp_x_axis_thresh2D = points6Warning2D.GetXaxis()
Exp_x_axis_thresh2D.CenterLabels(kTRUE)
Exp_y_axis_thresh2D = points6Warning2D.GetYaxis()
Exp_y_axis_thresh2D.CenterLabels(kTRUE)



		
		
Exp_thresh1D.Write()
Exp_thresh2D.Write()
Sugg_TDAC1D.Write()
Sugg_TDAC2D.Write()
points6Warning2D.Write()
	
	
#	
#
##CanvGraphPix = TCanvas 
#
#GraphWarning = TGraph()
#GraphWarning.SetPoint(0,0,857)
#GraphWarning.SetPoint(1,1,925)
#GraphWarning.SetPoint(2,2,1054)
#GraphWarning.SetPoint(3,3,1128)
#GraphWarning.SetPoint(4,4,1214)
#GraphWarning.SetPoint(5,5,1294)
#GraphWarning.SetPoint(6,6,1433)
#GraphWarning.SetPoint(7,7,1500)
#GraphWarning.SetPoint(8,8,854)
#GraphWarning.SetPoint(9,9,930)
#GraphWarning.SetPoint(10,10,853)
#GraphWarning.SetPoint(11,11,1130)
#GraphWarning.SetPoint(12,12,1215)
#GraphWarning.SetPoint(13,13,1294)
#GraphWarning.SetPoint(14,14,1430)
#GraphWarning.SetPoint(15,15,1517)
#
#Xaxis = GraphWarning.GetXaxis()
#Xaxis.SetLimits(0.,16.)
#GraphWarning.GetHistogram().SetMaximum(2000.)
#GraphWarning.GetHistogram().SetMinimum(0.)
#
#GraphWarning.Draw("AC*")
#outFile.cd()
#GraphWarning.Write()
#
#GraphWarning2 = TGraph()
#GraphWarning2.SetPoint(0,0,849)
#GraphWarning2.SetPoint(1,1,925)
#GraphWarning2.SetPoint(2,2,1020)
#GraphWarning2.SetPoint(3,3,1094)
#GraphWarning2.SetPoint(4,4,1225)
#GraphWarning2.SetPoint(5,5,1300)
#GraphWarning2.SetPoint(6,6,1406)
#GraphWarning2.SetPoint(7,7,1500)
#GraphWarning2.SetPoint(8,8,850)
#GraphWarning2.SetPoint(9,9,925)
#GraphWarning2.SetPoint(10,10,1019)
#GraphWarning2.SetPoint(11,11,1097)
#GraphWarning2.SetPoint(12,12,1225)
#GraphWarning2.SetPoint(13,13,1305)
#GraphWarning2.SetPoint(14,14,1400)
#GraphWarning2.SetPoint(15,15,1489)
#
#Xaxis = GraphWarning2.GetXaxis()
#Xaxis.SetLimits(0.,16.)
#GraphWarning2.GetHistogram().SetMaximum(2000.)
#GraphWarning2.GetHistogram().SetMinimum(0.)
#
#GraphWarning2.Draw("AC*")
#outFile.cd()
#GraphWarning2.Write()
#
#
## -- Canvas warning pixels
#CanvWarningPixels = ROOT.TCanvas("WarningPixels","WarningPixels",400,400)
#CanvWarningPixels.cd()
#
#
#gStyle.SetOptStat("e")
#gStyle.SetStatX(0.9)
#gStyle.SetStatY(0.95)
#
#points6Warning2D.GetZaxis().SetTitleOffset(1.3)
#gPad.SetRightMargin(0.15)
#
#points6Warning2D.Draw("colz")
#list_TLine2=[]
#for i in range(0,25):
#	line_vert =  TLine(i,12,i,48)
#	line_vert.Draw()
#	list_TLine2.append(line_vert)
#
#for i in range(12,49):
#	line_hor =  TLine(0,i,24,i)
#	line_hor.Draw()
#	list_TLine2.append(line_hor)
#	
#gPad.Update()
#
#
#CanvWarningPixels.Write()
#
#
#
#
#
#
#












# Write TDACfile.txt (must be at the end of the code) 

cnt=0
for i in file_TDAC_list:
	
	if cnt <1439:
		file_TDAC.write(""+str(file_TDAC_list[cnt])+"\r\n")    # \r in order to cope windows
	else:                                                           #necessary to avoid the last blank line in the output file
		file_TDAC.write(str(file_TDAC_list[cnt]))              #check if this is necessary for the software ---> it is not TODO: remove	
	cnt += 1
file_TDAC.close()



#------------------------------ notes

#other
#execfile("~/Documents/CERN\ master\ thesis/python\ scripts/ThresholdScan/CCPDv2/t.pyt.py")
#os.system("root -l") # +str(sys.argv[2])+"")

#os.system("2+2")


#execfile("t.py")
