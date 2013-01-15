#!/usr/bin/env python
from ROOT import *
from Parameters import *
from Functions import *
import numpy as np


class Detector:
  def __init__(self, name, energy_binning):
    self.fName = name

    self.fMass = self._AssignMass()

    self.fBadDays = []

    self.fEventList = [[], [], []] #time & energy list of added events

    self.fRunType = 'Bckgd'

    self.fVoltageFlagList = GetVoltageFlagList(self.fName, self.fRunType)

    #self.fERAConstants = GetERAConstants(self.fName)

    self.fEnergyBinning = energy_binning

    self.fGoodUnixPeriods = self._ReadInGoodPeriods()

    self.fGoodUnixPeriodsGap = self._InsertGaps(self.fGoodUnixPeriods) #fill gaps between periods

    self.fYearBins = self._CalcYearBins(self.fGoodUnixPeriodsGap) #get numpy array with time bins in years

    # histograms
    self.fHeatBaseline = TH1D(self.fName+'_heat_baseline', self.fName+' heat baseline;Time (years);Energy (keV_{nr})', \
                                  self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fHeatThreshold = TH1D(self.fName+'_heat_threshold', self.fName+' heat threshold;Time (years);Energy (keV_{nr})', \
                                   self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fFiducialTopBaseline = TH1D(self.fName+'_fiducial_top_baseline', self.fName+' fiducial top baseline;Time (years);Energy (keV_{nr})', \
                                         self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fFiducialBottomBaseline = TH1D(self.fName+'_fiducial_bottom_baseline', self.fName+' fiducial bottom baseline;Time (years);Energy (keV_{nr})', \
                                            self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fVetoTopBaseline = TH1D(self.fName+'_veto_top_baseline', self.fName+' veto top baseline;Time (years);Energy (keV_{nr})', \
                                     self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fVetoBottomBaseline = TH1D(self.fName+'_veto_bottom_baseline', self.fName+' veto bottom baseline;Time (years);Energy (keV_{nr})', \
                                        self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fGuardTopBaseline = TH1D(self.fName+'_guard_top_baseline', self.fName+' guard top baseline;Time (years);Energy (keV_{nr})', \
                                      self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fGuardBottomBaseline = TH1D(self.fName+'_guard_bottom_baseline', self.fName+' guard bottom baseline;Time (years);Energy (keV_{nr})', \
                                         self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fFiducialMeanBaseline = TH1D(self.fName+'_fiducial_mean', self.fName+' fiducial mean baseline;Time (years);Energy (keV_{nr})', \
                                         self.fYearBins.size-1, self.fYearBins.flatten('C'))

    self.fVoltage = TH1D(self.fName+'_voltage', self.fName+' voltage;Time (years);Voltage (V)', \
                             self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fLivetimeEfficiency = TH1D(self.fName+'_livetime_efficiency', self.fName+' Livetime Efficiency;Time (years) ;Efficiency', \
                                        self.fYearBins.size-1, self.fYearBins.flatten('C'))
    self.fTriggerEfficiency = TH2D(self.fName+'_trigger_efficiency', self.fName+' Trigger Efficiency;Time (years); Energy (keV_{nr}); Efficiency', \
	                                 self.fYearBins.size-1, self.fYearBins.flatten('C'), self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fFiducialEfficiencyGamma = TH2D(self.fName+'_fiducial_efficiency_gamma', self.fName+' Fiducial Efficiency Gamma;Time (years); Energy (keV_{nr}); Efficiency', \
	                                     self.fYearBins.size-1, self.fYearBins.flatten('C'), self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fFiducialEfficiencyNeutron = TH2D(self.fName+'_fiducial_efficiency_neutron', self.fName+' Fiducial Efficiency Neutron;Time (years); Energy (keV_{nr}); Efficiency', \
	                                         self.fYearBins.size-1, self.fYearBins.flatten('C'), self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fTotalEfficiencyGamma = TH2D(self.fName+'_total_efficiency_gamma', self.fName+' Total Efficiency Gamma;Time (years); Energy (keV_{nr}); Efficiency', \
                                     self.fYearBins.size-1, self.fYearBins.flatten('C'), self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fTotalEfficiencyNeutron = TH2D(self.fName+'_total_efficiency_neutron', self.fName+' Total Efficiency Neutron;Time (years); Energy (keV_{nr}); Efficiency', \
	                                 self.fYearBins.size-1, self.fYearBins.flatten('C'), self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])

    self.fIntegratedTriggerEfficiency = TH1D(self.fName+'_integrated_trigger_efficiency', self.fName+' Integrated Trigger Efficiency;Energy (keV_{nr});Efficiency', \
                                      self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fIntegratedFiducialEfficiencyGamma = TH1D(self.fName+'_integrated_fiducial_efficiency_gamma', self.fName+' Integrated Fiducial Efficiency Gamma;Energy (keV_{nr});Efficiency', \
                                      self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fIntegratedTotalEfficiencyGamma = TH1D(self.fName+'_integrated_total_efficiency_gamma', self.fName+' Integrated Total Efficiency Gamma;Energy (keV_{nr});Efficiency', \
                                      self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fIntegratedFiducialEfficiencyNeutron = TH1D(self.fName+'_integrated_fiducial_efficiency_neutron', self.fName+' Integrated Fiducial Efficiency Neutron;Energy (keV_{nr});Efficiency', \
                                      self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])
    self.fIntegratedTotalEfficiencyNeutron = TH1D(self.fName+'_integrated_total_efficiency_neutron', self.fName+' Integrated Total Efficiency Neutron;Energy (keV_{nr});Efficiency', \
                                      self.fEnergyBinning[0], self.fEnergyBinning[1], self.fEnergyBinning[2])


    self._FillHistograms(self.fGoodUnixPeriods)
    self._CalcTriggerEfficiency()
    self._CalcFiducialEfficiency('gamma')
    self._CalcFiducialEfficiency('neutron')
    self._CalcTotalEfficiency('gamma')
    self._CalcTotalEfficiency('neutron')
    self._CalcIntegratedEnergyEfficiency('trigger', 'gamma')
    self._CalcIntegratedEnergyEfficiency('fiducial', 'gamma')
    self._CalcIntegratedEnergyEfficiency('total', 'gamma')
    self._CalcIntegratedEnergyEfficiency('trigger', 'neutron')
    self._CalcIntegratedEnergyEfficiency('fiducial', 'neutron')
    self._CalcIntegratedEnergyEfficiency('total', 'neutron')



  def _AssignMass(self):
    mass = Masses[self.fName]['Mass']
    return mass['Fiducial']


  def _ReadInGoodPeriods(self):
    inpath = 'Run12PeriodInformation/'

    UnixStart, UnixEnd, Run, Heat_Baseline, Heat_Threshold, Fiducial_Top, Fiducial_Bottom, Veto_Top, Veto_Bottom, Guard_Top, Guard_Bottom = [], [], [], [], [], [], [], [], [], [], []

    infile = open(inpath+self.fName+'_Bckgd.txt','r')

    ELEC = EricsLowEnergyCuts[self.fName]
    for line in infile:
      unixstart, unixend, run, periodflag, heat_baseline, heat_threshold, fiducial_top, fiducial_bottom, veto_top, veto_bottom, guard_top, guard_bottom = line.split()

      startday = unixtime_to_day(int(unixstart))
      endday = unixtime_to_day(int(unixend))
      if GetVoltageFlag(self.fVoltageFlagList, str(run)) >= 0 \
        and startday not in self.fBadDays \
        and endday not in self.fBadDays \
        and ELEC['Coll1']['Min'] < float(fiducial_top) < ELEC['Coll1']['Max'] \
        and ELEC['Coll2']['Min'] < float(fiducial_bottom) < ELEC['Coll2']['Max'] \
        and ELEC['Veto1']['Min'] < float(veto_top) < ELEC['Veto1']['Max'] \
        and ELEC['Veto2']['Min'] < float(veto_bottom) < ELEC['Veto2']['Max'] \
        and ELEC['Guard1']['Min'] < float(guard_top) < ELEC['Guard1']['Max'] \
        and ELEC['Guard2']['Min'] < float(guard_bottom) < ELEC['Guard2']['Max'] \
        and ELEC['Heat']['Min'] < float(heat_baseline) < ELEC['Heat']['Max']:
	  UnixStart.append(int(unixstart))
	  UnixEnd.append(int(unixend))
	  Run.append(str(run))
	  Heat_Baseline.append(float(heat_baseline))
	  Heat_Threshold.append(float(heat_threshold))
	  Fiducial_Top.append(float(fiducial_top))
	  Fiducial_Bottom.append(float(fiducial_bottom))
	  Veto_Top.append(float(veto_top))
	  Veto_Bottom.append(float(veto_bottom))
	  Guard_Top.append(float(guard_top))
	  Guard_Bottom.append(float(guard_bottom))
    infile.close()

    return [UnixStart, UnixEnd, Run, Heat_Baseline, Heat_Threshold, Fiducial_Top, Fiducial_Bottom, Veto_Top, Veto_Bottom, Guard_Top, Guard_Bottom]


  def _InsertGaps(self, InList): #put gaps in between periods if not adjacent
    OutList = [ [ InList[i][0] ]  for i in range(len(InList))] #copy first row of values to new list
    for i in range(1,len(InList[0])):
      unixstart = InList[0][i]
      unixstop = InList[1][i]
      lastunixstop = InList[1][i-1]
      if unixstart != lastunixstop: #if gap between two periods
	OutList[0].append(lastunixstop) #insert last unixtime end as start
	OutList[1].append(unixstart) #insert next unixtime start as end
	for j in range(2,len(InList)): OutList[j].append(0.) #insert '0' for other values
      OutList[0].append(unixstart)
      OutList[1].append(unixstop)
      for j in range(2,len(InList)): OutList[j].append( InList[j][i] )
    return OutList


  def _CalcYearBins(self, InList):
    UnixTimeStart = InList[0]
    OutList = []
    for i in range (len(UnixTimeStart)):
      unixtime = UnixTimeStart[i]
      year = unixtime_to_year(unixtime)
      OutList.append(year)
    return np.array(OutList, dtype=np.float)


  def _FillHistograms(self, InList):
    for i in range(len(InList[0])):
      unixtimestart = InList[0][i]
      unixtimeend = InList[1][i]
      run = InList[2][i]
      heat_baseline = InList[3][i]
      heat_threshold = InList[4][i]
      fiducial_top = InList[5][i]
      fiducial_bottom = InList[6][i]
      veto_top = InList[7][i]
      veto_bottom = InList[8][i]
      guard_top = InList[9][i]
      guard_bottom = InList[10][i]

      voltageflag = GetVoltageFlag(self.fVoltageFlagList, run)
      #voltage = self.fERAConstants['gVolts'][voltageflag]
      voltage = 6.4

      time = unixtime_to_year((unixtimeend+unixtimestart)/2)

      self.fHeatBaseline.Fill(time, heat_baseline)
      self.fHeatThreshold.Fill(time, heat_threshold)
      self.fFiducialTopBaseline.Fill(time, fiducial_top)
      self.fFiducialBottomBaseline.Fill(time, fiducial_bottom)
      self.fVetoTopBaseline.Fill(time, veto_top)
      self.fVetoBottomBaseline.Fill(time, veto_bottom)
      self.fGuardTopBaseline.Fill(time, guard_top)
      self.fGuardBottomBaseline.Fill(time, guard_bottom)

      fiducial_mean = (fiducial_top * fiducial_bottom)/sqrt(pow(fiducial_top, 2) + pow(fiducial_bottom, 2))
      self.fFiducialMeanBaseline.Fill(time, fiducial_mean)

      self.fVoltage.Fill(time, voltage)
      self.fLivetimeEfficiency.Fill(time, 1)


  def _CalcTriggerEfficiency(self):
    for xbin in range(1,self.fTriggerEfficiency.GetNbinsX()+1):
      heat_baseline_ee = self.fHeatBaseline.GetBinContent(xbin)
      heat_threshold_ee = self.fHeatThreshold.GetBinContent(xbin)
      voltage = self.fVoltage.GetBinContent(xbin)

      heat_baseline_nr = GetEnergyRecoilFromEstimator(heat_baseline_ee, voltage)
      heat_threshold_nr = GetEnergyRecoilFromEstimator(heat_threshold_ee, voltage)

      heat_sigma_nr = heat_baseline_nr / (2*sqrt(2*log(2)))

      EfficiencyCurve = ThresholdEfficiency
      EfficiencyCurve.SetParameter(0, heat_threshold_nr)
      EfficiencyCurve.SetParameter(1, heat_sigma_nr)

      for ybin in range(1, self.fTriggerEfficiency.GetNbinsY()+1):
	mean_energy = self.fTriggerEfficiency.GetYaxis().GetBinCenter(ybin)
	value = EfficiencyCurve.Eval(mean_energy)

	if value >= 0 and heat_baseline_ee != 0 and heat_threshold_ee != 0:
	  efficiency = value
	else:
	  efficiency = 0

	self.fTriggerEfficiency.SetBinContent(xbin, ybin, efficiency)
	self.fTriggerEfficiency.SetBinError(xbin, ybin, 0)
    return None


  def _CalcFiducialEfficiency(self, qtype):
    if qtype == 'gamma':
      fid_eff_hist = self.fFiducialEfficiencyGamma
    elif qtype == 'neutron':
      fid_eff_hist = self.fFiducialEfficiencyNeutron


    for xbin in range(1,fid_eff_hist.GetNbinsX()+1): #time bin loop
      fiducial_baseline_ee = self.fFiducialMeanBaseline.GetBinContent(xbin)
      voltage = self.fVoltage.GetBinContent(xbin)

      E = pow(((2*fiducial_baseline_ee)/0.16),(1/1.18))
      q = Q_Lindhard.Eval(E)

      if fiducial_baseline_ee != 0:
	if qtype == 'gamma':
	  factor = GetEnergyRecoilFromEstimator(2*fiducial_baseline_ee, voltage)/(2*fiducial_baseline_ee)
	elif qtype == 'neutron':
	  factor = 1/q
      else:
	factor = 0

      fiducial_sigma = fiducial_baseline_ee / (2*sqrt(2*log(2)))

      # define fiducial efficiency curve for given time bin with according baseline values
      EfficiencyCurve = ThresholdEfficiency
      EfficiencyCurve.SetParameter(0, 2*fiducial_baseline_ee*factor)
      EfficiencyCurve.SetParameter(1, fiducial_sigma*factor)

      for ybin in range(1, fid_eff_hist.GetNbinsY()+1): #energy bin loop
	mean_energy = fid_eff_hist.GetYaxis().GetBinCenter(ybin)

	value = EfficiencyCurve.Eval(mean_energy)
	if value >= 0 and fiducial_baseline_ee != 0:
	  efficiency = value
	else:
	  efficiency = 0.0

	fid_eff_hist.SetBinContent(xbin, ybin, efficiency)
	fid_eff_hist.SetBinError(xbin, ybin, 0.0)
    return None


  def _CalcTotalEfficiency(self, qtype):
    if qtype == 'gamma':
      fid_eff_hist = self.fFiducialEfficiencyGamma
      total_eff_hist = self.fTotalEfficiencyGamma
    elif qtype == 'neutron':
      fid_eff_hist = self.fFiducialEfficiencyNeutron
      total_eff_hist = self.fTotalEfficiencyNeutron

    for xbin in range(1,total_eff_hist.GetNbinsX()+1):
      for ybin in range(1, self.fTriggerEfficiency.GetNbinsY()+1):
	livetime_eff = self.fLivetimeEfficiency.GetBinContent(xbin)
	trigger_eff = self.fTriggerEfficiency.GetBinContent(xbin, ybin)
	fiducial_eff = fid_eff_hist.GetBinContent(xbin, ybin)

	total_eff = livetime_eff * trigger_eff * fiducial_eff

	total_eff_hist.SetBinContent(xbin, ybin, total_eff)
	total_eff_hist.SetBinError(xbin, ybin, 0)
    return None


  def _CalcIntegratedEnergyEfficiency(self, eff, qtype): #calculate average energy efficiency hist, integrated over time
    if eff == 'trigger':
      inhist = self.fTriggerEfficiency
      outhist = self.fIntegratedTriggerEfficiency
    elif eff == 'fiducial':
      if qtype == 'gamma':
	inhist = self.fFiducialEfficiencyGamma
	outhist = self.fIntegratedFiducialEfficiencyGamma
      if qtype == 'neutron':
	inhist = self.fFiducialEfficiencyNeutron
	outhist = self.fIntegratedFiducialEfficiencyNeutron
    elif eff == 'total':
      if qtype == 'gamma':
	inhist = self.fTotalEfficiencyGamma
	outhist = self.fIntegratedTotalEfficiencyGamma
      if qtype == 'neutron':
	inhist = self.fTotalEfficiencyNeutron
	outhist = self.fIntegratedTotalEfficiencyNeutron

    StartTime = inhist.GetXaxis().GetBinLowEdge(1)
    EndTime = inhist.GetXaxis().GetBinLowEdge(inhist.GetNbinsX()+1)
    TotalTime = EndTime-StartTime

    for ybin in range(1, inhist.GetNbinsY()+1):
      Weighted_Efficiency_Sum = 0
      for xbin in range(1,inhist.GetNbinsX()+1):
	timewidth = inhist.GetXaxis().GetBinWidth(xbin)
	efficiency = inhist.GetBinContent(xbin, ybin)
	weighted_efficiency = timewidth * efficiency
	Weighted_Efficiency_Sum += weighted_efficiency
      average_efficiency = Weighted_Efficiency_Sum / TotalTime
      outhist.SetBinContent(ybin, average_efficiency)
    return None


  def WriteEventList(self):
    Emin = self.fEnergyBinning[1]
    Emax = self.fEnergyBinning[2]
    outfile = open('%s_%d-%d_keV_{nr}.txt' % (self.fName, Emin, Emax), 'w')
    for i in range(len(self.fEventList[0])):
      energy = self.fEventList[0][i]
      time = self.fEventList[1][i]
      outfile.write("%s %s \n" % (energy, time))
    outfile.close()
    return True


  def IsGoodEvent(self, UnixTime, Energy):
    Emin = self.fEnergyBinning[1]
    Emax = self.fEnergyBinning[2]
    if Emin <= Energy <= Emax:
      for entry in range(len(self.fGoodUnixPeriods[0])): #loop over the list with good periods
	if UnixTime >= self.fGoodUnixPeriods[0][entry] and UnixTime <= self.fGoodUnixPeriods[1][entry]:
	  return True
      return false
    else:
      return False


  def AddEvent(self, time, energy1, energy2):
    self.fEventList[0].append(time)
    self.fEventList[1].append(energy1)
    self.fEventList[2].append(energy2)


  def GetName(self):
    return self.fName


  def GetMass(self):
    return self.fMass


  def GetLivetime(self): #in days
    return self.fLivetimeEfficiency.Integral('WIDTH')*365


  def GetExposure(self): # in kg*days
    return self.GetLivetime()*self.fMass


  def GetLivetimeEfficiency(self):
    return self.fLivetimeEfficiency


  def GetFiducialEfficiency(self, qtype):
    if qtype == 'gamma': return self.fFiducialEfficiencyGamma
    elif qtype == 'neutron': return self.fFiducialEfficiencyNeutron


  def GetTriggerEfficiency(self, qtype):
    if qtype == 'gamma': return self.fTriggerEfficiency
    elif qtype == 'neutron': return self.fTriggerEfficiency


  def GetTotalEfficiency(self, qtype):
    if qtype == 'gamma': return self.fTotalEfficiencyGamma
    elif qtype == 'neutron': return self.fTotalEfficiencyNeutron


  def GetNuclearRecoilEnergy(self, energy_ee, VoltageFlag):
    #voltage = self.fERAConstants['gVolts'][VoltageFlag]
    voltage = 6.4
    energy = GetEnergyRecoilFromEstimator(energy_ee, voltage)
    return energy


  def GetAllBaselines(self):
    outlist = [
      self.fHeatBaseline,
      self.fHeatThreshold,
      self.fFiducialTopBaseline,
      self.fFiducialBottomBaseline,
      self.fVetoTopBaseline,
      self.fVetoBottomBaseline,
      self.fGuardTopBaseline,
      self.fGuardBottomBaseline
    ]
    return outlist


  def GetTimeBinning(self):
    return self.fYearBins


  def GetIntegratedEnergyEfficiency(self, name, qtype):
    if name == 'trigger':
      return self.fIntegratedTriggerEfficiency
    elif name == 'fiducial':
      if qtype == 'gamma': return self.fIntegratedFiducialEfficiencyGamma
      elif qtype == 'neutron': return self.fIntegratedFiducialEfficiencyNeutron
    elif name == 'total':
      if qtype == 'gamma': return self.fIntegratedTotalEfficiencyGamma
      elif qtype == 'neutron': return self.fIntegratedTotalEfficiencyNeutron


  def GetWeightedAverage(self, name):
    if name == 'heat': hist = self.fHeatBaseline
    elif name == 'threshold': hist = self.fHeatThreshold
    elif name == 'fiducialtop': hist = self.fFiducialTopBaseline
    elif name == 'fiducialbottom': hist = self.fFiducialBottomBaseline
    elif name == 'fiducialmean': hist = self.fFiducialMeanBaseline
    elif name == 'vetotop': hist = self.fVetoTopBaseline
    elif name == 'vetobottom': hist = self.fVetoBottomBaseline
    elif name == 'guardtop': hist = self.fGuardTopBaseline
    elif name == 'guardbottom': hist = self.fGuardBottomBaseline

    TotalTime = self.fLivetimeEfficiency.Integral('WIDTH')
    Weighted_Value_Sum = 0
    for xbin in range(1,hist.GetNbinsX()+1):
      time = hist.GetXaxis().GetBinWidth(xbin)
      value = hist.GetBinContent(xbin)
      weighted_value = time * value
      Weighted_Value_Sum += weighted_value
    average_value = Weighted_Value_Sum / TotalTime
    return average_value
