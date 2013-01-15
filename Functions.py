# ###############################################
# collection of functions
# ###############################################


import datetime as dt
import calendar
from ROOT import TF1, exp, log
from Parameters import *
from couchdbkit import Server, Database


# global
global inpath; inpath = 'Run12PeriodInformation/'


# date and time binning functions
March1st = dt.datetime(2009, 3, 1, 0, 0, 0)
OffsetJules = dt.timedelta(hours=-1)
global dayzero; DayZero = March1st + OffsetJules


def unixtime_to_year(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return (delta.days/365.24)+(delta.seconds/(365.24*86400))


def unixtime_to_day(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return delta.days


def unixtime_to_dayhour(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return delta.days + delta.seconds/86400.


def unixtime_to_hour(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  delta = date - DayZero
  return float(delta.seconds)/3600


def unixtime_to_date(unixtime):
  date = dt.datetime.utcfromtimestamp(unixtime)
  return date


def day_to_date(days):
  delta = dt.timedelta(days)
  date = DayZero + delta
  return date


def date_to_day(date):
  delta = date - DayZero
  return delta.days


def day_to_unixtime(day):
  date = day_to_date(day)
  unixtime = calendar.timegm(date.timetuple())
  return unixtime

def date_to_unixtime(date):
  unixtime = calendar.timegm(date.timetuple())
  return unixtime


# root functions


# Lindhard quenching relation (nuclear recoil)
Q_Lindhard = TF1('lindhard_quenching', '[0]*(x^[1])', 0, 100)
Q_Lindhard.SetParName(0, 'a')
Q_Lindhard.SetParName(1, 'b')
Q_Lindhard.FixParameter(0, 0.16)
Q_Lindhard.FixParameter(1, 0.18)
Q_Lindhard.SetNpx(1000)
Q_Lindhard.SetTitle('Lindhard Quenching Nuclear Recoil;E_{Recoil} [keV];Q(E_{Rec})')


# Recoil energy estimator for nuclear recoils
ERecEstimator = TF1('recoil_energy_estimator', '(x/(1+[0]/[1]))*(1+[0]/[1]*0.16*x^0.18)', 0, 100)
ERecEstimator.SetParName(0, 'Voltage')
ERecEstimator.SetParName(1, 'Creation Potential')
ERecEstimator.FixParameter(1, 3.0)
ERecEstimator.SetNpx(1000)
ERecEstimator.SetTitle('E_{Rec} estimator from E_{Heat};E_{Rec} [keV_{nr}];E_{Heat} [keV_{ee}]')


# threshold efficiency from threshold energy and threshold resolution
ThresholdEfficiency = TF1('threshold_efficiency', '0.5*(1+ROOT::Math::erf(((x-[0])/([1]*sqrt(2)))))', 0, 100)
ThresholdEfficiency.SetNpx(1000)
ThresholdEfficiency.SetParName(0, 'Threshold Energy')
ThresholdEfficiency.SetParName(1, 'Energy Resolution')
ThresholdEfficiency.SetTitle('Threshold Efficiency;E_{Heat} [keV];Efficiency')


# experimental ionization efficiency curve for ID3
IonizationEfficiency = TF1('ionization_efficiency', '0.95*(1-exp([0]*(x-[1])))', 0, 100)
IonizationEfficiency.SetNpx(1000)
IonizationEfficiency.SetParName(0, 'Parameter0')
IonizationEfficiency.SetParName(1, '2xFWHM')
IonizationEfficiency.SetTitle('Ionization Efficiency;E_{Ion} [keV];Efficiency')


def GetEnergyRecoilFromEstimator(energy_ee, voltage):
  function = ERecEstimator
  function.SetParameter(0, voltage)
  function.FixParameter(1, 3.0) #electron-hole creation potential
  energy_recoil = function.GetX(energy_ee)
  return energy_recoil


def GetVoltageFlagList(Detector, Runtype):
  RunName, VoltageFlag = [], []
  infile = open(inpath+Detector+'_voltage_'+Runtype+'.txt', 'r')
  for line in infile:
    runname, voltage = line.split()
    RunName.append(str(runname))
    VoltageFlag.append(int(voltage))
  infile.close()
  return [RunName, VoltageFlag]


def GetVoltageFlag(InList, RunName):
  try:
    flag = InList[1][InList[0].index(RunName)]
    return flag
  except ValueError:
    print "error finding run",RunName
    return 0


def GetSambaNumber(Detector):
  if Detector in ['ID5', 'ID403', 'FID401']:
    return 'S1'
  elif Detector in ['ID3', 'ID401', 'FID402']:
    return 'S2'
  elif Detector in ['ID4', 'ID6', 'ID404']:
    return 'S3'
  elif Detector in ['ID2', 'ID402', 'ID405']:
    return 'S4'


def GetERAConstants(Detector):
  s = Server('http://edwdbik.fzk.de:5984')
  db = s['analysis']
  vr = db.view('constants/run12', include_docs=True)
  constantDoc = {}
  for row in vr:
    constantDoc[row['key']] = row['doc']
  return constantDoc[Detector]