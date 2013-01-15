# ###############################################
# main analysis script with MLE
# ###############################################


#!/usr/bin/env python
from ROOT import *
from Functions import *
from Parameters import *
from pyWIMP.DMModels.base_model import BaseVariables
from pyWIMP.DMModels.flat_model import FlatModel
import pyWIMP.DMModels.wimp_model as wimp_model

from DetectorClass import *


#WIMP_Masses = [8, 10, 12, 15, 20, 25, 30]
WIMP_Masses = [6, 8, 10, 15, 20, 30]
rmin = 0; rmax = 1; rstart = 0.5


# choice of detector and energybinning
DetectorName = 'ID3'
EnergyBinning = (100, 0, 10)
ID = Detector(DetectorName, EnergyBinning)
mass_of_detector = ID.GetMass() #kg

m_chi_initial = 10 #just to initialize pyWIMP


# times
TimeBins = ID.GetTimeBinning(); minTime = TimeBins[0]; maxTime = TimeBins[-1]


# energies
energyBins, minEnergy, maxEnergy = EnergyBinning[0], EnergyBinning[1], EnergyBinning[2]



# parameter container
basevars = BaseVariables(minTime, maxTime, minEnergy, maxEnergy, True, time_offset=-2./365.25) #offset to start at 1st Mar 2009
time = basevars.get_time(); energy = basevars.get_energy()


# wimp modulation model & pdf
wm = wimp_model.WIMPModel(basevars, m_chi_initial, mass_of_detector, constant_quenching=True, nucl_recoil=True)


# efficiency hists
gamma_efficiency_hist = ID.GetTotalEfficiency('gamma')
gamma_efficiency_datahist = RooDataHist('gamma_efficiency_datahist', 'Total Efficiency DataHist Gamma',RooArgList(time, energy), gamma_efficiency_hist)
gamma_efficiency_pdf = RooHistPdf('gamma_efficiency_pdf', 'Total Efficiency Gamma', RooArgSet(time, energy), gamma_efficiency_datahist)

neutron_efficiency_hist = ID.GetTotalEfficiency('neutron')
neutron_efficiency_datahist = RooDataHist('neutron_efficiency_datahist', 'Total Efficiency DataHist Neutron',RooArgList(time, energy), neutron_efficiency_hist)
neutron_efficiency_pdf = RooHistPdf('neutron_efficiency_pdf', 'Total Efficiency Neutron', RooArgSet(time, energy), neutron_efficiency_datahist)


# define background pdf
background_hist = gamma_efficiency_hist
background_pdf = gamma_efficiency_pdf



#-----------------------------------------------------------------------------------------------------
# kdata loop

# empty data set
data = RooDataSet(DetectorName+'_dataset', DetectorName+' events', RooArgSet(energy, time))

kdatafile = 'Run12_ID3+6_bckg_with_subrecords.root'

events = 0

infile = KDataReader(kdatafile)
event = infile.GetEvent()
entries = infile.GetEntries()
for entry in range(entries):
  infile.GetEntry(entry)

  stamp = event.GetStampTime()

  for entry2 in range(event.GetNumBolos()):

    #bolo records
    bolo = event.GetBolo(entry2)
    EnergyRecoil = bolo.GetEnergyRecoil()
    EnergyIon = bolo.GetEnergyIonFiducial()
    EnergyHeat = bolo.GetEnergyHeat(1)
    Detector = bolo.GetDetectorName()
    EventFlag = bolo.GetEventFlag()
    VoltageFlag = bolo.GetVoltageFlag()

    #samba records
    samba = bolo.GetSambaRecord()
    UnixTime = samba.GetNtpDateSec()
    RunName = samba.GetRunName()

    if Detector == ID.GetName() and EventFlag == 2:
     Energy = ID.GetNuclearRecoilEnergy(EnergyHeat, VoltageFlag)

     if ID.IsGoodEvent(UnixTime, Energy) == True:
	ID.AddEvent(UnixTime, Energy, EnergyIon)

	energy.setVal(Energy)
	time.setVal(unixtime_to_year(UnixTime))

	data.add(RooArgSet(energy, time))

	events += 1


infile.Close()
print "KData:",events,"events read"

#-----------------------------------------------------------------------------------------------------


WIMP_PDF = []
F_NORM = []
WIMP_HIST = []
WIMP_HIST = []
SIGNAL_DATAHIST = []
SIGNAL_HIST = []
SIGNAL_PDF = []
RATIO = []
FINAL_PDF = []
R, R_MAX, R_90 = [], [], []
XS = []
N_CHI, N_SIG, N_OBS = [], [], []
NLL = []
MINUIT = []

counter = 0
# loop over wimp masses
for m_chi in WIMP_Masses:

  WIMP_PDF.append(wm.get_WIMP_model_with_escape_vel(m_chi))
  F_NORM.append(wm.get_normalization().getVal())
  N_CHI.append(WIMP_PDF[-1].expectedEvents(RooArgSet(time, energy))) # number of expected events per nucleus cross section [pb] for detector mass in time and energy range


  # create and fill histogram
  WIMP_HIST.append(TH2D('pywimp_hist_%GeV' % m_chi, 'pyWIMP Hist;Time (years); Energy (keV); evts/pb', TimeBins.size-1, TimeBins.flatten('C'), energyBins, minEnergy, maxEnergy))
  pywimp_hist = WIMP_HIST[-1]
  for xbin in range(1,pywimp_hist.GetNbinsX()+1):
    t_low = pywimp_hist.GetXaxis().GetBinLowEdge(xbin)
    t_high = pywimp_hist.GetXaxis().GetBinLowEdge(xbin+1)
    t_mean = 0.5*(t_low+t_high)
    time.setVal(t_mean)
    for ybin in range(1, pywimp_hist.GetNbinsY()+1):
      e_low = pywimp_hist.GetYaxis().GetBinLowEdge(ybin)
      e_high = pywimp_hist.GetYaxis().GetBinLowEdge(ybin+1)
      e_mean = 0.5*(e_low+e_high)
      energy.setVal(e_mean)
      pdf_val = WIMP_PDF[-1].getVal() # value in middle of time and energy bin (!) only good for fine binning
      pywimp_hist.SetBinContent(xbin, ybin, pdf_val)

  #WIMP_DATAHIST.append(RooDataHist('pywimp_datahist', 'pywimp_datahist',RooArgList(time, energy), pywimp_hist))
  #WIMP_HIST_PDF.append(RooHistPdf('pywimphist_pdf', 'pywimphist_pdf', RooArgSet(time, energy), WIMP_DATAHIST[-1]))


  SIGNAL_HIST.append(pywimp_hist.Clone('wimp_signal_%GeV' % m_chi))
  SIGNAL_HIST[-1].Multiply(neutron_efficiency_hist)
  SIGNAL_DATAHIST.append(RooDataHist('wimp_signal_datahist_%GeV' % m_chi, 'wimp signal datahist',RooArgList(time, energy), SIGNAL_HIST[-1]))
  SIGNAL_PDF.append(RooHistPdf('wimp_signal_pdf_%GeV' % m_chi, 'wimp signal pdf', RooArgSet(time, energy), SIGNAL_DATAHIST[-1]))


  N_SIG.append(SIGNAL_HIST[-1].Integral('WIDTH')) # number of events per nucleus cross section after application of efficiencies

  # RooFit signal ratio
  RATIO.append(RooRealVar('ratio_%GeV' % m_chi,'signal fraction',rstart,rmin,rmax))

  # efficiency pdfs
  FINAL_PDF.append(RooAddPdf('final_pdf_%GeV' % m_chi,'r*signal+(1-r)*bckgd',RooArgList(SIGNAL_PDF[-1], background_pdf),RooArgList(RATIO[-1])))


  #-----------------------------------------------------------------------------------------------------
  ## normal fit procedure, not needed when using nll method below
  #FINAL_PDF[-1].fitTo(data, RooFit.Save(kTRUE), RooFit.PrintEvalErrors(2), RooFit.Minos(kTRUE))

  # Construct function object representing log(L)
  NLL.append(RooNLLVar('nll','nll',FINAL_PDF[-1],data, RooFit.PrintEvalErrors(2)))
  MINUIT.append(RooMinuit(NLL[-1]))
  MINUIT[-1].fit('hvr')


  # calculation of cross section
  #R.append(RATIO[-1].getVal())
  R.append(RATIO[-1].getAsymErrorHi())
  R_MAX.append(RATIO[-1].getError())
  R_90.append(R[-1] + 1.64 * R_MAX[-1])
  N_OBS.append(R_90[-1] * events)
  XS.append((N_OBS[-1]/N_SIG[-1])*F_NORM[-1])

  counter += 1


# print table with all values
i = 0
print "fit results:"
print '{0:8} | {1:8} | {2:8} | {3:8} | {4:8} | {5:8} | {6:8} | {7:8}'.format('m_chi','r','r_max','r_90','N_max','N_obs','N_chi','xs')
for mass_of_wimp in WIMP_Masses:
  print '{0:8} | {1:8.2g} | {2:8.2g} | {3:8.2g} | {4:8.2g} | {5:8.2g} | {6:8.2g} | {7:8.1e}'.format(mass_of_wimp,R[i],R_MAX[i],R_90[i],N_OBS[i],N_CHI[i],N_SIG[i],XS[i])
  i += 1