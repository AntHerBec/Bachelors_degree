from __future__ import print_function
from Selector import Selector
from Plotter import Plotter
import ROOT as r
import time as time

r.gROOT.SetBatch(1) # To work only by batch (i.e. through terminal, w/o windows)

# Create a selector for ttbar sample
#mySelector = Selector('ttbar')
# Obtain the histogram of invariant mass of the muons from the selector
#h = mySelector.GetHisto('InvMass')
# Do whatever with it...
#print 'The histogram InvMass for the ttbar sample has %i entries' %h.GetEntries()


# Create an object plotter thar contains all MC samples and uses data:
MCsamples = ['qcd', 'wjets', 'ww', 'wz', 'zz', 'dy', 'single_top', 'ttbar']

t0 = time.time()
plot = Plotter(MCsamples, 'data')
t1 = time.time()

tiempo_plotter = t1 - t0

print('\n', 'Tiempo de ejecución: {tiempo:.2f} s'.format(tiempo = tiempo_plotter), '\n')

# Set colors for each process... We can use colors defined in ROOT
colors = [r.kGray, r.kBlue-1, r.kTeal-1, r.kTeal+1, r.kTeal+4, r.kAzure-8, r.kOrange+1, r.kRed+1]
plot.SetColors(colors)

# Set other plotting options (Legend position, size, titles, etc...)
plot.SetXtitle('p_{T}^{#mu} (GeV)')
plot.SetYtitle('Events')
plot.SetTitle('')

# Draw the stack plot for data and simulation
plot.Stack('MuonPt_ttbar')

plot.SetXtitle('Number of jets')
plot.Stack('NJet_ttbar')


plot.SetXtitle('Number of Muons')
plot.Stack('NMuon_ttbar')

plot.SetXtitle('b-tagging')
plot.Stack('btag_ttbar')


# Print the contributions for each background and observed data
plot.PrintCounts('NJet_ttbar')
plot.XS('NJet_ttbar')
