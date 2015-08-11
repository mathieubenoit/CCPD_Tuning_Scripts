# VERSION 4
# analyses noise occupancy from STcontrol and suggest plot histogtam of sugg GDAC
# SYS ARGUMENTS: 1 --> GDAC value list (for a fixed TDAC)

#call the script:

#noise_occ_gDAC_suggestion test.root [0,4,7,11,15] [[0.79,0.792,0.794,0.796,0.798,0.8,0.802,0.804,0.806,0.808,0.81,0.812,0.814,0.816,0.818,0.82,0.822,0.824,0.826,0.828,0.83],[0.806,0.808,0.81,0.812,0.814,0.816,0.818,0.82,0.822,0.824,0.826,0.828,0.83,0.832,0.834,0.836,0.838,0.84,0.842],[0.81,0.812,0.814,0.816,0.818,0.82,0.822,0.824,0.826,0.828,0.83,0.832,0.834,0.836,0.838,0.84,0.842,0.844,0.846,0.848,0.85],[0.826,0.828,0.83,0.832,0.834,0.836,0.838,0.84,0.842,0.844,0.846,0.848,0.85,0.852,0.854,0.856,0.858,0.86,0.862,0.864],[0.84,0.842,0.844,0.846,0.848,0.85,0.852,0.854,0.856,0.858,0.86,0.862,0.864,0.866,0.868,0.87,0.88]]

import sys
import ROOT
from ROOT import gROOT, TCanvas, TF1, TGraph, TLegend, TMath, TMultiGraph, TFile, TH2F, TH1F, TDirectory,TH2D, TH1D, TLatex, gPad, TH2I, TLine, TAxis, kTRUE, TKey, THStack, kRed, kOrange, kBlue, kGreen
from array import array
import math
import ast # to get list in sys argv

#-----inputs-----
rootFileName = sys.argv[1]
GDACList = ast.literal_eval(sys.argv[3]) # get list from sys arg, LIST COMPONENT MUST BE A LIST
TDACScanList = ast.literal_eval(sys.argv[2])
#-----dic-------
data_dic={}
noise_dic={}
NoiseEdgeDic = {}
SuggTDACDicFEFrame = {}
SuggTDACDicCMOSFrame = {}
#----TFile-----
outFile = TFile(rootFileName,"RECREATE")
for i in TDACScanList:
	TDACdir = outFile.mkdir("TDAC_"+str(i)+"")
	GraphNoiseEdgeDir = TDACdir.mkdir("GraphsNoiseEdge") # TODO: en rajouter
	TGraphNoiseEdgeDir = TDACdir.mkdir("TGraphsNoiseEdge")
GraphGDACSuggDir = outFile.mkdir("GraphsGDACSugg")
outFile.cd()
#------plots--------
#noiseocc_map_FE = TH2F("noiseocc_map_FE","noise occ map FE;#Col;#Row",6,0.5,6.5,16,146.5,162.5)
#noiseocc_map_hvcmos = TH2F("noiseocc_map_hvcmos","noise occ map hvcmos;#Row;#Col",24,-0.5,23.5,36,11.5,47.5)

#TargetGDAC = 0.8339


#------prerequisites:

# Create TDAC list compatible with Ivan's soft
file_TDAC_list = []
lineNum_TDAC=0
for counter in range(1,1441): # old way: create a TDAC all 6 file and modify it
	if counter < 289:
		file_TDAC_list.append(0)	# to be deleted
	elif 288 < counter < 1153:
		file_TDAC_list.append(0) # to be deleted
	elif counter > 1152:
		file_TDAC_list.append(0) # to be deleted

#--------Functions

#-----Definition changeframe FE --> HVCMOS map ------
cntROW = 0
HVCMOScol = 0
HVCMOSrow = 0
def changeframe(FEcol,FErow):
	global cntROW
	global HVCMOScol
	global HVCMOSrow
	cntROW = 0
	for n in range(8):
		if FErow == 147 + 2*n or FErow == 148 + 2*n : 
			HVCMOScol = 24 + cntROW	# VNC1: 24, VNC2: 25, VNC3: 26
		for m in range(6):
			if FEcol == 1 + m and FErow == 147 + 2*n :
				HVCMOSrow = 11 - 2*m
			elif FEcol == 1 + m and FErow == 148 + 2*n : 
				HVCMOSrow = 10 - 2*m
		cntROW += 3
	return (HVCMOScol,HVCMOSrow)

#-----Set TDAC value in a list which will be written in a txt file ------
def setTDAC(CMOSrow,CMOScol,newTDAC_value):
	cnt = 576
	for c in range(24,48):
		for r in range(12):
			if r == CMOSrow and c == CMOScol:
				file_TDAC_list[cnt] = newTDAC_value
				print "c = ", c , "r = ",r, " cnt = ", cnt, " TDAC = ", newTDAC_value
#			else:
#				print "row and column not found"
			cnt += 1
		cnt += 12








#-----import histogram from ST control and fill dic------
cnt = 0
for TDACScanVal in TDACScanList:
	for GDACVal in GDACList[cnt]:
		name = "NO_VNC1_TDAC"+str(TDACScanVal)+"_Thr"+str(GDACVal)+".root"	
		fileTDAC = ROOT.TFile(name)
		canvasTDAC = fileTDAC.Get("pixcan")
		k = TKey()
		k = fileTDAC.GetKey("pixcan")
		histoTDAC = k.ReadObj().FindObject("Occup_0_00_MA")

		for col in range(147,163):
			for row in range(1,7):
				noise=histoTDAC.GetBinContent(row,col)
				noise_dic["GDAC"+str(GDACVal)+"row"+str(row)+"_col"+str(col)+"_TDAC"+str(TDACScanVal)+""] = noise
	cnt += 1

#-----TGraph Noise vs GDAC------
cnt = 0
for TDACScanVal in TDACScanList:
	NoiseEdgeHist = TH1F("NoiseEdgeHist","NoiseEdgeHist;NoiseEdgeValue;#pixels",50,0.8,0.9)
	for col in range(147,163):
		for row in range(1,7):
			NoiseList = []
			for GDACVal in GDACList[cnt]:
				NoiseVal = noise_dic["GDAC"+str(GDACVal)+"row"+str(row)+"_col"+str(col)+"_TDAC"+str(TDACScanVal)+""]
				NoiseList.append(NoiseVal)
			
			y=array("d",NoiseList)
			x=array("d",GDACList[cnt])
			GDAC_vs_TDACTGraph = TGraph(len(x),x,y)			


			#-----Noise edge analysis------
			if all(i >= 10 for i in NoiseList) == True: # check if all noise point are above the selection (here 10)
				print 'TDAC ',TDACScanVal,'r',row,' c',col
				print len(NoiseList)-1, "CAUTION"
				NoiseEdgeX = len(NoiseList)-1 # TODO: this is really not a good way since we believe now that this pixel is not noisy 
				SuggGDAC = GDACList[cnt][NoiseEdgeX]					
			else:
				NoiseEdgeX = next(x[0] for x in enumerate(NoiseList) if x[1] < 10)
				SuggGDAC = GDACList[cnt][NoiseEdgeX]
				
			NoiseEdgeDic["row"+str(row)+"_col"+str(col)+"_TDAC"+str(TDACScanVal)+""] = SuggGDAC
			
			#-----TGraph options------
			GDAC_vs_TDACTGraph.SetNameTitle("row"+str(row)+"_col"+str(col)+"","row"+str(row)+"_col"+str(col)+"")		
			GDAC_vs_TDACTGraph.SetTitle("row"+str(row)+"_col"+str(col)+"_referenceGDAC = "+str(SuggGDAC)+"")
			XaxisTGraph = GDAC_vs_TDACTGraph.GetXaxis()
			XaxisTGraph.SetLimits(0.7,0.9)
#			GDAC_vs_TDACTGraph.GetHistogram().SetMaximum(0.9)          
#			GDAC_vs_TDACTGraph.GetHistogram().SetMinimum(0.7)				
			outFile.cd("TDAC_"+str(TDACScanVal)+"/TGraphsNoiseEdge")			
			GDAC_vs_TDACTGraph.Draw("ALP*")
			GDAC_vs_TDACTGraph.Write()
			
			#-----NoiseEdge distribution------
			NoiseEdgeHist.Fill(SuggGDAC)

	outFile.cd("TDAC_"+str(TDACScanVal)+"")
	NoiseEdgeHist.Write()
	cnt += 1
	
	#-----Target the GDAC----
	if TDACScanList.count(7) == 0:
		print 'Error: TDAC 7 FILE NOT FOUND'
		print 'Run not completed, script stopped'
		sys.exit() # exit script	
	elif TDACScanVal == 7:		
		TargetGDAC = NoiseEdgeHist.GetMean(1)




#-----hist noise vs GDAC------						
#cnt = 0
#for TDACScanVal in TDACScanList:
#	NoiseEdgeHist = TH1F("NoiseEdgeHist","NoiseEdgeHist;NoiseEdgeValue;#pixels",50,0.8,0.9)
#	for col in range(147,163):
#		for row in range(1,7):
#			DataPointsList = []
#			Noise_vs_GDACHist = TH1F("Noise_vs_GDACHist","Noise_vs_GDACHist;GDAC;Noise",100,0.7,0.9)
#
#			for GDACVal in GDACList[cnt]:
#				DataPoint = int(noise_dic["GDAC"+str(GDACVal)+"row"+str(row)+"_col"+str(col)+"_TDAC"+str(TDACScanVal)+""])
#	#			print DataPoint
#				DataPointsList.append(DataPoint)
#				Noise_vs_GDACHist.Fill(GDACVal,DataPoint)
#				
#			print DataPointsList
#
#			if DataPointsList.count(0) > 0:
#				print next(x[0] for x in enumerate(DataPointsList) if x[1] < 10)
#				NoiseEdgeX = next(x[0] for x in enumerate(DataPointsList) if x[1] < 10)
#				SuggGDAC = GDACList[cnt][NoiseEdgeX]
#			else:
#				print len(DataPointsList)-1, "CAUTION"
#				NoiseEdgeX = len(DataPointsList)-1 # TODO: this is really not a good way since we believe now that this pixel is not noisy 
#				SuggGDAC = GDACList[cnt][NoiseEdgeX]
#			print 'sugg GDAC = ', SuggGDAC
#			print ''
#			Noise_vs_GDACHist.SetNameTitle("row"+str(row)+"_col"+str(col)+"_SuggGDAC = "+str(SuggGDAC)+"","row"+str(row)+"_col"+str(col)+"_SuggGDAC = "+str(SuggGDAC)+"")
#			Noise_vs_GDACHist.SetTitle("row"+str(row)+"_col"+str(col)+"_SuggGDAC = "+str(SuggGDAC)+"")
#			outFile.cd("TDAC_"+str(TDACScanVal)+"/GraphsNoiseEdge")
#	#		Noise_vs_GDACHist.Draw()	
#			Noise_vs_GDACHist.Write()
#	#		Noise_vs_GDACHist.Delete()		
#					
#			NoiseEdgeDic["row"+str(row)+"_col"+str(col)+"_TDAC"+str(TDACScanVal)+""] = SuggGDAC
#			NoiseEdgeHist.Fill(SuggGDAC)
#	outFile.cd("TDAC_"+str(TDACScanVal)+"")
#	NoiseEdgeHist.Write()
#	cnt += 1

#-----GDAC vs TDAC plot------
ExpGDACDist = TH1F("ExpGDACDist","ExpGDACDist;ExpGDAC;#pixels",100,0.8,0.9)
for col in range(147,163):
	for row in range(1,7):
		GDACPointsList = []
		for TDACScanVal in TDACScanList:
			GDACVal = NoiseEdgeDic["row"+str(row)+"_col"+str(col)+"_TDAC"+str(TDACScanVal)+""]
			GDACPointsList.append(GDACVal)
			
#			print 'col ',col,' row ',row
#			print GDACPointsList
#			print ''
#			
			
			#GDAC_vs_TDACGraph = TGraph()

		y=array("d",GDACPointsList)
		x=array("d",TDACScanList)
		GDAC_vs_TDACGraph = TGraph(len(x),x,y)
		

		
		#-----TDACSugg------
		AbsoluteDiffList = []
		DiffList = []
		for i in GDACPointsList:
			AbsoluteDiffList.append(abs(TargetGDAC - i))
			DiffList.append(TargetGDAC - i)
		
		FirstDiffListPoint = min(AbsoluteDiffList)
		FirstPointIndex = AbsoluteDiffList.index(FirstDiffListPoint)
		FirstGDACPoint = GDACPointsList[FirstPointIndex]
		FirstTDACPoint = TDACScanList[FirstPointIndex]
		
		if DiffList[FirstPointIndex] > 0: # if target threshold is HIGHER than the closest point
			SecondPointIndex = FirstPointIndex + 1
			SecondGDACPoint = GDACPointsList[SecondPointIndex]
			SecondTDACPoint = TDACScanList[SecondPointIndex]
			
			p = (SecondGDACPoint - FirstGDACPoint) / (SecondTDACPoint - FirstTDACPoint)
			b = (FirstGDACPoint - p * FirstTDACPoint)
			SuggTDACFloat = (TargetGDAC - b)/p
			SuggTDACInt = int(round(SuggTDACFloat,0))
			ExpGDAC = p * SuggTDACInt + b
			
			RangeMin = FirstTDACPoint
			RangeMax = SecondTDACPoint
			MyFitInterpol=TF1("MyFitInterpol", "pol1",RangeMin,RangeMax)
			GDAC_vs_TDACGraph.Fit("MyFitInterpol","R")
			p_fit = MyFitInterpol.GetParameter(1)
			b_fit = MyFitInterpol.GetParameter(0)
			
			print p,'--->',p_fit
			print b,'--->',b_fit
			
		elif DiffList[FirstPointIndex] < 0: # if target threshold is LOWER than the closest point
			SecondPointIndex = FirstPointIndex - 1
			SecondGDACPoint = GDACPointsList[SecondPointIndex]
			SecondTDACPoint = TDACScanList[SecondPointIndex]

			p = (FirstGDACPoint - SecondGDACPoint) / (FirstTDACPoint - SecondTDACPoint)
			b = (SecondGDACPoint - p * SecondTDACPoint)
			SuggTDACFloat = (TargetGDAC - b)/p
			SuggTDACInt = int(round(SuggTDACFloat,0))
			ExpGDAC = p * SuggTDACInt + b

			RangeMin = FirstTDACPoint
			RangeMax = SecondTDACPoint
			MyFitInterpol=TF1("MyFitInterpol", "pol1",RangeMin,RangeMax)
			GDAC_vs_TDACGraph.Fit("MyFitInterpol","R")
			p_fit = MyFitInterpol.GetParameter(1)
			b_fit = MyFitInterpol.GetParameter(0)
			
			print p,'--->',p_fit
			print b,'--->',b_fit

		elif DiffList[FirstPointIndex] == 0:
			SuggTDACInt = int(TDACScanList[FirstPointIndex])
			ExpGDAC = FirstGDACPoint
		else:
			print 'Error: SecondGDACPoint not found'
			print 'Run not completed, script stopped'
			sys.exit() # exit script
		
		
		#---convert FE frame --> CMOS frame
		CMOSFrame = changeframe(row,col) # WARNING: row and col inverted...
		CMOScol = CMOSFrame[0]
		CMOSrow = CMOSFrame[1]
		
		SuggTDACDicFEFrame ["row"+str(row)+"_col"+str(col)+""] = SuggTDACInt
		SuggTDACDicCMOSFrame ["CMOSrow"+str(CMOSrow)+"_CMOScol"+str(CMOScol)+""] = SuggTDACInt
		
		#----set suggested TDAC value in list
		setTDAC(CMOSrow,CMOScol,SuggTDACInt)
		
		
		#-----TGraph options------
		GDAC_vs_TDACGraph.SetNameTitle("r"+str(row)+"_c"+str(col)+"_CMOSrow"+str(CMOSrow)+"_CMOScol"+str(CMOScol)+"","row"+str(row)+"_col"+str(col)+"")		
		GDAC_vs_TDACGraph.SetTitle("row"+str(row)+"_col"+str(col)+"_TargetGDAC = "+str(TargetGDAC)+"_FirstGDACPoint = "+str(FirstGDACPoint)+"_FirstTDACPoint = "+str(FirstTDACPoint)+"_SecondGDACPoint = "+str(SecondGDACPoint)+"_SecondTDACPoint = "+str(SecondTDACPoint)+"_SuggTDACFloat = "+str(SuggTDACFloat)+"_SuggTDACInt = "+str(SuggTDACInt)+"_ExpGDAC = "+str(ExpGDAC)+"")
		Xaxis = GDAC_vs_TDACGraph.GetXaxis()
		Xaxis.SetLimits(-1,16)
#		GDAC_vs_TDACGraph.GetHistogram().SetMaximum(0.9)          
#		GDAC_vs_TDACGraph.GetHistogram().SetMinimum(0.7)				
		GDAC_vs_TDACGraph.SetMaximum(0.865)          
		GDAC_vs_TDACGraph.SetMinimum(0.79)
		outFile.cd("GraphsGDACSugg")
		GDAC_vs_TDACGraph.Draw("AL*")
		GDAC_vs_TDACGraph.Write()
		
		#----ExpGDACDist
		ExpGDACDist.Fill(ExpGDAC)
		
outFile.cd()
ExpGDACDist.Write()	


# Write TDACfile.txt (must be at the end of the code) 
file_TDAC = open("SuggTDACFileNoiseTuning.txt", "w")
for i in file_TDAC_list:
	file_TDAC.write(""+str(i)+"\n")
file_TDAC.close()

#ROOT.gApplication.Run()

