from __future__ import print_function
import os
import warnings as wr
import ROOT as r

class Selector:
    ''' Class to do an event selection'''
    ### =============================================
    ### Constructor
    def __init__(self, filename = ''):
        ''' Initialize a new Selector by giving the name of a sample.root file '''
        self.name = filename
        self.filename = self.name
        if self.filename[-5:] != '.root': self.filename += '.root'
        if not os.path.exists(self.filename): self.filename = 'files/' + self.filename
        if not os.path.exists(self.filename):
            if (self.name != ''): wr.warn("[Selector::constructor] WARNING: file {f} not found".format(f = self.name))
        else:
            self.CreateHistograms()
            self.Loop()
        return

    ### =============================================
    ### Attributes
    # General
    histograms = []
    name = ""
    filename = ""
    incert_est_cuadrado = 0 # voy sumando los pesos al cuadrado


    # Variables de histogramas

    muon_pt_ttbar = -99         # Sólo exigimos un muón, no es necesario que sea lista
    btag_ttbar = []             # Varios jets, debe ser una lista
    nmuon_ttbar = -99
    njet_ttbar = -99

    # Pesos y variables útiles para imponer los criterios

    buenos_muones = -99         # Contador para muones con pT > 28
    buenos_jets = -99           # Contador para jets con pT > 20
    bjets_buenPt = -99          # Contador para b-jets con pT > 20
    weight = -99


    ### =============================================
    ### Methods
    def CreateHistograms(self):
        ''' CREATE YOUR HISTOGRAMS HERE '''

        self.histograms.append(r.TH1F(self.name + '_MuonPt_ttbar',     ';p_{T}^{#mu} (GeV);Events', 20, 0, 200))
        self.histograms.append(r.TH1F(self.name + '_NJet_ttbar',     ';Number of Jets;Events', 7, 0, 7))
        self.histograms.append(r.TH1F(self.name + '_btag_ttbar',     ';btag;Events', 32, -2, 6))
        self.histograms.append(r.TH1F(self.name + '_NMuon_ttbar',     ';#eta;Events', 5, 0, 5))
        return


    def GetHisto(self, name):
        ''' Use this method to access to any stored histogram '''
        for h in self.histograms:
            n = h.GetName()
            if self.name + '_' + name == n: return h
        wr.warn("[Selector::GetHisto] WARNING: histogram {h} not found.".format(h = name))
        return r.TH1F()

 
    def FillHistograms(self):

        self.GetHisto('MuonPt_ttbar').Fill(self.muon_pt_ttbar,    self.weight)
        self.GetHisto('NJet_ttbar').Fill(self.njet_ttbar,          self.weight)
        self.GetHisto('NMuon_ttbar').Fill(self.nmuon_ttbar,            self.weight)
        for i in self.btag_ttbar:
            self.GetHisto('btag_ttbar').Fill(i,              self.weight)

        return


    def Loop(self):
        ''' Main method, used to loop over all the entries '''
        f = r.TFile.Open(self.filename)
        tree = f.events

        nEvents = tree.GetEntries()

        limitante = int(nEvents)             # Para comprobaciones, trunca el bucle al llegar al número de eventos indicado

        print("Opening file {f} and looping over {n} events...".format(f = self.filename, n = nEvents))

        for event in tree:

            limitante -= 1
            if limitante == 0: break

            self.buenos_muones = 0        # Establezco los contadores a 0
            self.buenos_jets = 0
            self.bjets_buenPt = 0

            if not event.triggerIsoMu24: continue  # Aplica el trigger a los eventos

            if event.NMuon != 0:

                for i in range(event.NMuon):
                    muon = r.TLorentzVector()
                    muon.SetPxPyPzE(event.Muon_Px[i], event.Muon_Py[i], event.Muon_Pz[i], event.Muon_E[i])

                    if muon.Pt() >= 28:
                        self.buenos_muones += 1

                if self.buenos_muones != 1: continue   # Selecciona eventos con 1 muón bueno

                self.nmuon_ttbar = event.NMuon
                self.muon_pt_ttbar = muon.Pt()

                if event.NElectron != 0: continue      # Selecciona eventos donde el leptón sea un muón


                self.btag_ttbar = []
                contador_btags = 0

                for i in range(event.NJet):
                    jet = r.TLorentzVector()
                    jet.SetPxPyPzE(event.Jet_Px[i], event.Jet_Py[i], event.Jet_Pz[i], event.Jet_E[i])
                    if jet.Pt() >= 20:
                        self.buenos_jets += 1
                        if event.Jet_btag[i] >= 2:        # Establece el criterio de btag considerado para jets: de 1 a 5 cada vez mas probable, -1 es que ha fallado
                            self.btag_ttbar.append(event.Jet_btag[i])
                            contador_btags += 1

                if self.buenos_jets not in [3,4]: continue    # Criterio equivalente a 'buenos_muones' pero con jets ademas del b-tag
                if contador_btags not in [1,2]: continue      # Criterio para la cantidad de jets buenos con b-tagging suficiente

                self.njet_ttbar = event.NJet

                self.weight = event.EventWeight if not self.name == 'data' else 1
                self.incert_est_cuadrado += self.weight**2
                self.FillHistograms()

                '''
                -------------------------------------------------------------------------------------------------------------------------------------------------------------
                ## Intento de aplicar criterio de masa invariante de al menos 1 pareja de jets cercana a W,
                ## pero no cuento con información de generación acerca de la energía de los jets y no se me
                ## ocurre cómo asignarle un valor sin recurrir a variables de reconstrucción para recalcular
                ## la aceptancia
                -------------------------------------------------------------------------------------------------------------------------------------------------------------
                lista_jets_Minv = []                  # Para seleccionar eventos con alguna pareja de jets de masa invariante similar a la del W
                contador_Minv = 0

                for i in range(event.NJet):
                    jet = r.TLorentzVector()
                    jet.SetPxPyPzE(event.Jet_Px[i], event.Jet_Py[i], event.Jet_Pz[i], event.Jet_E[i])
                    if jet.Pt() >= 20:
                        self.buenos_jets += 1
                        if event.Jet_btag[i] >= 2:        # Establece el criterio de btag considerado para jets: de 1 a 5 cada vez mas probable, -1 es que ha fallado
                            self.btag.append(event.Jet_btag[i])
                            contador_btags += 1

                if self.buenos_jets not in [3,4]: continue
                if contador_btags not in [1,2]: continue
                
                for i,j in combinations(lista_jets_Minv, 2):    # Hace directamente las combinaciones posibles de jets en la lista, sin necesidad de recurrir a dobles bucles
                    suma = i + j
                    Minv = suma.M()
                    if abs(Minv - 80.4) < 30: contador_Minv += 1

                if contador_Minv == 0: continue
                
                '''

        return
