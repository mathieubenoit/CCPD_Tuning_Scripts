# VERSION 4
# analyses noise occupancy from STcontrol and suggest plot histogtam of sugg GDAC
# SYS ARGUMENTS: 1 --> GDAC value list (for a fixed TDAC)

import sys
#sys.path.append("/usr/local/lib/root")
import ROOT
from ROOT import gROOT, TCanvas, TF1, TGraph, TLegend, TMath, TMultiGraph, TFile, TH2F, TH1F, TDirectory,TH2D, TH1D, TLatex, gPad, TH2I, TLine, TAxis, kTRUE, TKey, THStack, kRed, kOrange, kBlue, kGreen
from array import array
import math
import ast # to get list in sys argv

#-----inputs-----
rootFileName = sys.argv[1]
GDACList = ast.literal_eval(sys.argv[3]) # get list from sys arg
print GDACList[0]
TDACScanList = ast.literal_eval(sys.argv[2])
TDACScanVal = TDACScanList[0]
#-----dic-------
data_dic={}
noise_dic={}
NoiseEdgeDic = {}
#----TFile-----
outFile = TFile(rootFileName,"RECREATE")
GraphDir = outFile.mkdir("Pixels")
outFile.cd()
#------plots--------
#noiseocc_map_FE = TH2F("noiseocc_map_FE","noise occ map FE;#Col;#Row",6,0.5,6.5,16,146.5,162.5)
#noiseocc_map_hvcmos = TH2F("noiseocc_map_hvcmos","noise occ map hvcmos;#Row;#Col",24,-0.5,23.5,36,11.5,47.5)
NoiseEdgeHist = TH1F("NoiseEdgeHist","NoiseEdgeHist;NoiseEdgeValue;#pixels",50,0.8,0.9)


for GDACVal in GDACList:
	name = "NO_VNC1_TDAC"+str(TDACScanVal)+"_Thr"+str(GDACVal)+".root"	
	fileTDAC = ROOT.TFile(name)
	canvasTDAC = fileTDAC.Get("pixcan")
	k = TKey()
	k = fileTDAC.GetKey("pixcan")
	histoTDAC = k.ReadObj().FindObject("Occup_0_00_MA")

	for col in range(147,163):
		for row in range(1,7):
			noise=histoTDAC.GetBinContent(row,col)
			noise_dic["GDAC"+str(GDACVal)+"row"+str(row)+"_col"+str(col)+""] = noise

for col in range(147,163):
	for row in range(1,7):
		DataPointsList = []
		Noise_vs_GDACHist = TH1F("Noise_vs_GDACHist","Noise_vs_GDACHist;GDAC;Noise",100,0.7,0.9)

		for GDACVal in GDACList:
			DataPoint = int(noise_dic["GDAC"+str(GDACVal)+"row"+str(row)+"_col"+str(col)+""])
#			print DataPoint
			DataPointsList.append(DataPoint)
			Noise_vs_GDACHist.Fill(GDACVal,DataPoint)
		print DataPointsList

		if DataPointsList.count(0) > 0:
			print next(x[0] for x in enumerate(DataPointsList) if x[1] < 10)
			NoiseEdgeX = next(x[0] for x in enumerate(DataPointsList) if x[1] < 10)
			SuggGDAC = GDACList[NoiseEdgeX]
		else:
			print len(DataPointsList)-1, "CAUTION"
			NoiseEdgeX = len(DataPointsList)-1
			SuggGDAC = GDACList[NoiseEdgeX]
		print 'sugg GDAC = ', SuggGDAC
		print ''
		Noise_vs_GDACHist.SetNameTitle("row"+str(row)+"_col"+str(col)+"_SuggGDAC = "+str(SuggGDAC)+"","row"+str(row)+"_col"+str(col)+"_SuggGDAC = "+str(SuggGDAC)+"")
		Noise_vs_GDACHist.SetTitle("row"+str(row)+"_col"+str(col)+"_SuggGDAC = "+str(SuggGDAC)+"")
		outFile.cd("Pixels")
#		Noise_vs_GDACHist.Draw()	
		Noise_vs_GDACHist.Write()
#		Noise_vs_GDACHist.Delete()		
		
		
		NoiseEdgeDic["row"+str(row)+"_col"+str(col)+""] = SuggGDAC
		NoiseEdgeHist.Fill(SuggGDAC)
		

outFile.cd()
NoiseEdgeHist.Write()
c1=TCanvas("SuggGDACHist","SuggGDACHist", 400,300)
NoiseEdgeHist.Draw()
gPad.Update()





#
#
#
#
#
#
#
#
#
#
#print 'name ',name 
#
#
#
#
#fileTDAC = ROOT.TFile(name)
#canvasTDAC = fileTDAC.Get("pixcan")
#k = TKey()
#k = fileTDAC.GetKey("pixcan")
#histoTDAC = k.ReadObj().FindObject("Occup_0_00_MA")
#
#
#for col in range(147,163):
#    for row in range(1,7):
#        noise=histoTDAC.GetBinContent(row,col)
#        noiseocc_map_FE.Fill(row,col,noise)
#        noise_dic["row"+str(row)+"_col"+str(col)+""]=noise
##        data_dic["TDAC"+str(numbTDAC)+"_row"+str(row)+"_col"+str(col)+""]=noise
#
##print 'ok'
##noiseocc_map_FE.SetNameTitle("noiseocc_"+name+"","noiseocc_"+name+"")
##noiseocc_map_FE.Write()
##print 'ok'
#
#
##c1=TCanvas("NoiseOccMap_FE","NoiseOccMap_FE", 400,300)
##c1.cd()
##noiseocc_map_FE.Draw("colztext")
##c1.SetRightMargin(0.36)
##noiseocc_map_FE.GetZaxis().SetTitleOffset(1.3)
##list_TLine3=[]
#
##for i in range(0,13):
##    line_vert =  TLine(i+0.5,141.5,i+0.5,165.5)
##    line_vert.Draw()
##    list_TLine3.append(line_vert)
#
##for i in range(141,166):
##    line_hor =  TLine(0.5,i+0.5,12.5,i+0.5)
##    line_hor.Draw()
##    list_TLine3.append(line_hor)
#
##c1.Modified()
#
#
##def ChangeFrame(FErow,FEcol,value):
##    for vert in range(12):
##        for hor in range(6):
##            for i in range(4):
#
#
#
#
#for i in range(142,166):
#    for j in range(1,13):
#        for n in range(12):
#            if i == 142 + 2*n:
#                HVCMOScol = i-130+n
#            elif i == 143+2*n:
#                HVCMOScol = i-131+n
#            
#            for m in range(6):
#                if j == 12-2*m and i == 142+2*n:
#                    HVCMOSrow = 1+4*m
#                    print 'FErow | FEcol | HVcol | HVrow'
#                    print i,' |  ',j,'  |  ', HVCMOScol,' | ', HVCMOSrow
#                    print ''
#                elif j == 12-2*m and i == 143+2*n:
#                    HVCMOSrow = 0+4*m
#                    print 'FErow | FEcol | HVcol | HVrow'
#                    print i,' |  ',j,'  |  ', HVCMOScol,'  |  ', HVCMOSrow
#                    print ''
#                elif j == 11-2*m and i == 142+2*n:
#                    HVCMOSrow = 2+4*m
#                    print 'FErow | FEcol | HVcol | HVrow'
#                    print i,' |  ',j,'  |  ', HVCMOScol,'  |  ', HVCMOSrow
#                    print ''
#                elif j == 11-2*m and i == 143+2*n:
#                    HVCMOSrow = 3+4*m
#                    print 'FErow | FEcol | HVcol | HVrow'
#                    print i,' |  ',j,'  |  ', HVCMOScol,'  |  ', HVCMOSrow
#                    print ''
#
#
#def changeframe(FEcol,FErow):
#    for n in range(12):
#        if FErow == 142 + 2*n:
#            HVCMOScol = FErow-130+n
#        elif FErow == 143+2*n:
#            HVCMOScol = FErow-131+n
#
#        for m in range(6):
#            if FEcol == 12-2*m and FErow == 142+2*n:
#                HVCMOSrow = 1+4*m
#            elif FEcol == 12-2*m and FErow == 143+2*n:
#                HVCMOSrow = 0+4*m
#            elif FEcol == 11-2*m and FErow == 142+2*n:
#                HVCMOSrow = 2+4*m
#            elif FEcol == 11-2*m and FErow == 143+2*n:
#                HVCMOSrow = 3+4*m
#    return  (HVCMOScol,HVCMOSrow)
#
#
#print changeframe(5,155)
#print changeframe(5,155)[1]
#
#print changeframe(1,142)[1]
#
#for col in range(142,166):
#    for row in range(1,13):
#        colhvframe = changeframe(row,col)[0]
#        rowhvframe = changeframe(row,col)[1]
#        noise=histoTDAC.GetBinContent(row,col)
#        noiseocc_map_hvcmos.Fill(rowhvframe,colhvframe,noise)
#
#
#c2=TCanvas("NoiseOccMap_HVCMOS","NoiseOccMap_HVCMOS", 1600,800)
##c2.cd()
#c2.Divide(2,1)
#c2.cd(1)
#
#noiseocc_map_FE.Draw("colztext")
#noiseocc_map_FE.GetZaxis().SetTitleOffset(1.3)
#list_TLine2=[]
#
#for i in range(0,7):
#    line_vert =  TLine(i+0.5,146.5,i+0.5,162.5)
#    line_vert.Draw()
#    list_TLine2.append(line_vert)
#
#for i in range(146,163):
#    line_hor =  TLine(0.5,i+0.5,6.5,i+0.5)
#    line_hor.Draw()
#    list_TLine2.append(line_hor)
#
#gPad.SetRightMargin(0.3)
#c2.cd(2)
#noiseocc_map_hvcmos.Draw("colz")
#noiseocc_map_FE.GetZaxis().SetTitleOffset(1.3)
#
#list_TLine=[]
#
#for i in range(0,24):
#    line_vert =  TLine(i-0.5,11.5,i-0.5,47.5)
#    line_vert.Draw()
#    list_TLine.append(line_vert)
#
#for i in range(11,48):
#    line_hor =  TLine(-0.5,i+0.5,23.5,i+0.5)
#    line_hor.Draw()
#    list_TLine.append(line_hor)
#
#gPad.SetRightMargin(0.3)
#c2.Modified()
#gPad.Update()
#
##----------------------- Writing output files
##outFile.cd()
###c1.Write()
##noiseocc_map_hvcmos.Write()
##c2.Write()
#
ROOT.gApplication.Run()
#
