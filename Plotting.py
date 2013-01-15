# ###############################################
# only code snippets to run after analysis script
# ###############################################



# #####################################################################
# plot natural log likelihood value (NLL) for range of ratio parameters
mass_of_wimp = 5

n = WIMP_Masses.index(mass_of_wimp)
ratio = RATIO[n]
nll = NLL[n]

c_nll = TCanvas('c_nll', 'c_nll', 800, 600)
rframe = ratio.frame()
nll.plotOn(rframe)
rframe.Draw()
# #####################################################################



# ##########################################
# plot signal for different energies
Masses = [6, 8, 10, 15, 20, 30]
Hist_List = []
c2 = TCanvas()
for mass in Masses:
  n = WIMP_Masses.index(mass)
  print mass, n
  hist = SIGNAL_HIST[n].ProjectionY('proj_'+str(n))
  hist.Scale(N_SIG[n]/hist.Integral())
  Hist_List.append(hist)
  if n==0: hist.Draw('HIST')
  else: hist.Draw('HISTSAME')

  hist.GetXaxis().SetTitle('Energy [keV_{nr}]')
  hist.GetYaxis().SetTitle('evts/pb/0.1 keV')

# additional background pdf
c2.Update()

background = background_hist.ProjectionY('background')
background.Scale(gPad.GetUymax())
background.SetLineColor(kGreen)

rightmax = 1.0
background.Draw('SAME')

axis1 = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),0,rightmax,510,"+L")
axis1.SetLineColor(kGreen)
axis1.SetLabelColor(kGreen)
axis1.SetTitleColor(kGreen)
axis1.SetTitle('Efficiency')
axis1.Draw()
# ##########################################



# ##########################################
# plot signal, background and data in energy
mass_of_wimp = 5

n = WIMP_Masses.index(mass_of_wimp)
signal = SIGNAL_PDF[n]
background = gamma_efficiency_pdf
final = FINAL_PDF[n]

c_energy = TCanvas('c_energy', 'c_energy', 800, 600)
eframe = energy.frame()
data.plotOn(eframe, RooFit.Binning(15))
background.plotOn(eframe, RooFit.LineColor(kGreen))
signal.plotOn(eframe, RooFit.LineColor(kRed))
data.plotOn(eframe, RooFit.Binning(15))
eframe.Draw()
# ##########################################



# #####################################################################
# plot exclusion curve
M_CHI_Eric = [8, 10, 12, 15, 20, 25, 30]
M_CHI_MLE = WIMP_Masses
XS_Eric = [1.17e-4, 1.35e-5, 4.17e-6, 1.55e-6, 7.34e-7, 5.47e-7, 4.85e-7]
XS_MLE = XS

graph_Eric = TGraph(len(XS_Eric))
graph_MLE = TGraph(len(XS_MLE))

for i in range(len(XS_Eric)):
  graph_Eric.SetPoint(i, M_CHI_Eric[i], XS_Eric[i])
for i in range(len(XS_MLE)):
  graph_MLE.SetPoint(i, M_CHI_MLE[i], XS_MLE[i])

c = TCanvas()
gPad.SetLogy()

graph1 = graph_Eric
graph2 = graph_MLE

graph1.Draw('ALP')

graph1.GetYaxis().SetTitle('\sigma_{\chi - nucleon} [pb]')
graph1.GetYaxis().CenterTitle()
graph1.GetXaxis().SetTitle('WIMP mass [GeV]')
graph1.GetXaxis().CenterTitle()
graph1.SetLineColor(kBlue)
graph1.SetLineWidth(3)

graph1.GetXaxis().SetLabelSize(0.05)
graph1.GetXaxis().SetTitleSize(0.05)
graph1.GetXaxis().SetTitleOffset(0.9)
graph1.GetYaxis().SetLabelSize(0.05)
graph1.GetYaxis().SetTitleSize(0.05)
graph1.GetYaxis().SetTitleOffset(0.9)

graph2.Draw('LP')
graph2.SetLineColor(kRed)
graph2.SetLineWidth(3)
# #####################################################################
