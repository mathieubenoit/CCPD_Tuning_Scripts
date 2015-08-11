# VERSION 4
# analyses noise occupancy and redo noise occupancy map in FEI4 frame AND HVCMOS frame
# calling this macro: python NoiseOcc_suggTDAC_1VNC.py OUTPUTROOTfile.root NoiseOccFILE.root
# WARNING ONLY FOR VNC1 ON 

import sys
#sys.path.append("/usr/local/lib/root")
import ROOT
from ROOT import gROOT, TCanvas, TF1, TGraph, TLegend, TMath, TMultiGraph, TFile, TH2F, TH1F, TDirectory,TH2D, TH1D, TLatex, gPad, TH2I, TLine, TAxis, kTRUE, TKey, THStack, kRed, kOrange, kBlue, kGreen, gStyle
from array import array
import math

#------------------

#rootFileName = sys.argv[1]
#print sys.argv
#outFile = TFile(rootFileName,"RECREATE")
#
#outFile.cd()
#---------------- Recreate noise occupancy scan FE frame
#-------------------------------------------------------

data_dic={}
noise_dic={}

name = sys.argv[1]

print 'name ',name 

noiseocc_map_FE = TH2F("noiseocc_map_FE","noise occ map FE;#Col;#Row",6,0.5,6.5,16,146.5,162.5)
noiseocc_map_hvcmos = TH2F("noiseocc_map_hvcmos","noise occ map hvcmos;#Row;#Col",12,-0.5,11.5,24,23.5,47.5)


fileTDAC = ROOT.TFile(name)
canvasTDAC = fileTDAC.Get("pixcan")
k = TKey()
k = fileTDAC.GetKey("pixcan")
histoTDAC = k.ReadObj().FindObject("Occup_0_00_MA")


for col in range(147,163):
    for row in range(1,7):
        noise=histoTDAC.GetBinContent(row,col)
        noiseocc_map_FE.Fill(row,col,noise)
        noise_dic["row"+str(row)+"_col"+str(col)+""]=noise
#        data_dic["TDAC"+str(numbTDAC)+"_row"+str(row)+"_col"+str(col)+""]=noise

#print 'ok'
#noiseocc_map_FE.SetNameTitle("noiseocc_"+name+"","noiseocc_"+name+"")
#noiseocc_map_FE.Write()
#print 'ok'


#c1=TCanvas("NoiseOccMap_FE","NoiseOccMap_FE", 400,300)
#c1.cd()
#noiseocc_map_FE.Draw("colztext")
#c1.SetRightMargin(0.36)
#noiseocc_map_FE.GetZaxis().SetTitleOffset(1.3)
#list_TLine3=[]

#for i in range(0,13):
#    line_vert =  TLine(i+0.5,141.5,i+0.5,165.5)
#    line_vert.Draw()
#    list_TLine3.append(line_vert)

#for i in range(141,166):
#    line_hor =  TLine(0.5,i+0.5,12.5,i+0.5)
#    line_hor.Draw()
#    list_TLine3.append(line_hor)

#c1.Modified()


#def ChangeFrame(FErow,FEcol,value):
#    for vert in range(12):
#        for hor in range(6):
#            for i in range(4):



#-----Definition changeframe FE --> HVCMOS map ------
def changeframe(FEcol,FErow):
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
	return  (HVCMOScol,HVCMOSrow)


for col in range(147,163):
    for row in range(1,7):
        colhvframe = changeframe(row,col)[0]
        rowhvframe = changeframe(row,col)[1]
        noise=histoTDAC.GetBinContent(row,col)
        noiseocc_map_hvcmos.Fill(rowhvframe,colhvframe,noise)


c2=TCanvas("NoiseOccMap_HVCMOS","NoiseOccMap_HVCMOS", 1200,600)
#c2.cd()
c2.Divide(2,1)
c2.cd(1)
gStyle.SetOptStat("e")
gStyle.SetStatY(0.95)
gStyle.SetStatX(0.9)
noiseocc_map_FE.Draw("colztext")
noiseocc_map_FE.GetZaxis().SetTitleOffset(1.2)
noiseocc_map_FE.GetZaxis().SetLabelSize(0.025)
noiseocc_map_FE.GetZaxis().SetLabelOffset(0.005)
noiseocc_map_FE.GetYaxis().SetTitleOffset(1.25)

list_TLine2=[]

for i in range(0,7):
    line_vert =  TLine(i+0.5,146.5,i+0.5,162.5)
    line_vert.Draw()
    list_TLine2.append(line_vert)

for i in range(146,163):
    line_hor =  TLine(0.5,i+0.5,6.5,i+0.5)
    line_hor.Draw()
    list_TLine2.append(line_hor)

gPad.SetRightMargin(0.15)
c2.cd(2)
noiseocc_map_hvcmos.Draw("colztext")
noiseocc_map_hvcmos.GetZaxis().SetTitleOffset(1.2)
noiseocc_map_hvcmos.GetZaxis().SetLabelSize(0.025)
noiseocc_map_hvcmos.GetZaxis().SetLabelOffset(0.005)

list_TLine=[]

for i in range(0,12):
    line_vert =  TLine(i-0.5,23.5,i-0.5,47.5)
    line_vert.Draw()
    list_TLine.append(line_vert)

for i in range(23,48):
    line_hor =  TLine(-0.5,i+0.5,11.5,i+0.5)
    line_hor.Draw()
    list_TLine.append(line_hor)

gPad.SetRightMargin(0.15)
c2.Modified()
gPad.Update()

#----------------------- Writing output files
#outFile.cd()
##c1.Write()
#noiseocc_map_hvcmos.Write()
#c2.Write()

ROOT.gApplication.Run()

