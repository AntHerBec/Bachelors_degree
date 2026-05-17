from __future__ import print_function
import warnings as wr
import os
from Selector import Selector
import ROOT as r
import numpy as np

class Plotter:
    ''' Class to draw histograms and get info from Selector'''
    ### =============================================
    ### Constructor
    def __init__(self, backgrounds, data = '', path = "./results"):
        ''' Initialize a new plotter... give a list with all names of MC samples and the name of the data sample '''
        self.data = data
        self.backgrounds = backgrounds
        self.savepath = path
        counter = 0
        for p in self.backgrounds:
            self.listOfSelectors.append(Selector(p))
            self.colors.append(counter+1)
            counter += 1

        if self.data != '': 
            self.dataSelector = Selector(self.data)
        return

    ### =============================================
    ### Attributes
    savepath = "."
    listOfSelectors = []
    dataSelector = Selector()
    data = ''
    colors = []

    # Default parameters
    fLegX1, fLegY1, fLegX2, fLegY2 = 0.75, 0.55, 0.89, 0.89
    LegendTextSize  = 0.035
    xtitle = ''
    ytitle = ''
    title  = ''


    ### =============================================
    ### Methods
    def SetLegendPos(self, x1, y1, x2, y2):
        ''' Change the default position of the legend'''
        self.fLegX1 = x1
        self.fLegY1 = y1
        self.fLegX2 = x2
        self.fLegY2 = y2


    def SetLegendSize(self, t = 0.065):
        ''' Change the default size of the text in the legend'''
        self.LegendTextSize = t


    def SetSavePath(self, newpath):
        '''Change where plots and text dumps are going to be saved. By default: ./results'''
        self.savepath = newpath


    def SetColors(self, col):
        ''' Set the colors for each MC sample '''
        self.colors = col


    def SetTitle(self, tit):
        ''' Set title of the plot '''
        self.title = tit


    def SetXtitle(self, tit):
        ''' Set title of X axis '''
        self.xtitle = tit


    def SetYtitle(self, tit):
        ''' Set title of Y axis '''
        self.ytitle = tit


    def GetHisto(self, process, name):
        ''' Returns histogram 'name' for a given process '''
        for s in self.listOfSelectors:
            if name not in s.name: continue
            h = s.GetHisto(name)
            return h

        wr.warn("[Plotter::GetHisto] WARNING: histogram {h} for process {p} not found!".format(h = name, p = process))
        return r.TH1F()


    def GetEvents(self, process, name):
        ''' Returns the integral of a histogram '''
        return self.GetHisto(process, name).Integral()


    def Stack(self, name):
        ''' Produce a stack plot for a variable given '''
        if (isinstance(name, list)):
            for nam in name: self.Stack(nam)
            return

        c = r.TCanvas('c_'+name, 'c', 10, 10, 800, 600)
        l = r.TLegend(self.fLegX1, self.fLegY1, self.fLegX2, self.fLegY2)
        l.SetTextSize(self.LegendTextSize)
        l.SetBorderSize(0)
        l.SetFillColor(10)

        hstack = r.THStack('hstack_' + name, "hstack")
        counter = 0
        for s in self.listOfSelectors:
            h = s.GetHisto(name)
            h.SetFillColor(self.colors[counter])
            h.SetLineColor(0)
            hstack.Add(h)
            l.AddEntry(h, s.name, "f")
            counter += 1
        hstack.Draw("hist")

        if self.title  != '': hstack.SetTitle(self.title)
        if self.xtitle != '': hstack.GetXaxis().SetTitle(self.xtitle)
        if self.ytitle != '': hstack.GetYaxis().SetTitle(self.ytitle)
        hstack.GetYaxis().SetTitleOffset(1.35)
        Max = hstack.GetStack().Last().GetMaximum();

        if self.data != '':
            hdata = self.dataSelector.GetHisto(name)
            hdata.SetMarkerStyle(20)
            hdata.SetMarkerColor(1)
            hdata.Draw("pesame")
            MaxData = hdata.GetMaximum()
            if(Max < MaxData): Max = MaxData
            l.AddEntry(hdata, self.dataSelector.name, 'p')
        l.Draw("same")
        hstack.SetMaximum(Max * 1.1)
        create_folder(self.savepath)
        c.Print(self.savepath + "/" + name + '.png', 'png')
        c.Print(self.savepath + "/" + name + '.pdf', 'pdf')
        c.Close()
        return


    def PrintCounts(self, name):
        ''' Print the number of events for each sample in a given histogram '''
        if (isinstance(name, list)):
            for nam in name: self.PrintEvents(nam)
            return

        print("\nPrinting number of events for histogram {h}:".format(h = name))
        print('----------------------------------------------------')
        total = 0.

        print('-------------------------------------')
        print("Number of events applying selection:")
        print('-------------------------------------')

        for s in self.listOfSelectors:
            h = s.GetHisto(name)
            print("{nam}: {num}".format(nam = s.name, num = h.Integral()))
            total += h.Integral()


        print('-------------------------------------')
        print('Expected (MC): {tot}'.format(tot = total))
        print('-------------------------------------')

        if self.data != '':
            hdata = self.dataSelector.GetHisto(name)

            print('-------------------------------------')
            print('Observed (data): {tot}'.format(tot = hdata.Integral()))
            print('-------------------------------------\n')

        return

    def XS(self, name):

        ### Cálculo de la aceptancia ###

        f = r.TFile("files/ttbar.root")
        t = f.Get("events")

        contador_total_ttbar = 0      # Contador de los sucesos ttbar totales (más general que establecer directamente los 36941)
        contador_ttbar_buenMC = 0     # Contador de los sucesos que cumplen los criterios fiduciales

        for event in t:

            contador_total_ttbar += 1

            # Definición con la información de generación y criterios para el muón (y para estar en el canal leptónico adecuado)
            if abs(event.MCleptonPDGid) != 13: continue

            muon = r.TLorentzVector(event.MClepton_px, event.MClepton_py, event.MClepton_pz, 0)  # Me da igual el valor de E para Pt, eta y phi
            if muon.Pt() < 28 or abs(muon.Eta()) >= 2.4: continue

            # Definición con la información de generación y criterios para los jets
            jet1 = r.TLorentzVector(event.MChadronicWDecayQuark_px, event.MChadronicWDecayQuark_py,
                                     event.MChadronicWDecayQuark_pz, 0)
            jet2 = r.TLorentzVector(event.MChadronicWDecayQuarkBar_px, event.MChadronicWDecayQuarkBar_py,
                                     event.MChadronicWDecayQuarkBar_pz, 0)
            bjet_h = r.TLorentzVector(event.MChadronicBottom_px, event.MChadronicBottom_py, event.MChadronicBottom_pz,
                                       0)
            bjet_l = r.TLorentzVector(event.MCleptonicBottom_px, event.MCleptonicBottom_py, event.MCleptonicBottom_pz,
                                       0)

            contador_buenos_jets = 0

            lista_jets = [jet1, jet2, bjet_h, bjet_l]

            for i in lista_jets:
                if i.Pt() >= 20 and abs(i.Eta()) < 2.4: contador_buenos_jets += 1

            if contador_buenos_jets < 3: continue
            contador_ttbar_buenMC += 1

        aceptancia = contador_ttbar_buenMC / contador_total_ttbar

        # ---------------------------------------------------------------------------------------------------------------

        ### Cálculo de la eficiencia de b-tagging ###

        def eta_y_phi(px, py, pz): # Función para calcular pseudorrapidez y ángulo phi azimutal para un TLorentzVector

            v = r.TLorentzVector(px, py, pz, 0)  # para eta y phi no me importa el valor de E

            eta = v.Eta()
            phi = v.Phi()

            return [eta, phi]

        contador_bjets_buenos = 0                 # Contador de b-jets (adecuadamente reconstruidos)
        contador_btag = 0                         # Contador esos b-jets que además son tageados correctamente

        for event in t:

            px_h, py_h, pz_h = event.MChadronicBottom_px, event.MChadronicBottom_py, event.MChadronicBottom_pz
            px_l, py_l, pz_l = event.MCleptonicBottom_px, event.MCleptonicBottom_py, event.MCleptonicBottom_pz

            eta_bjet_hadronico, phi_bjet_hadronico = eta_y_phi(px_h, py_h, pz_h)
            eta_bjet_leptonico, phi_bjet_leptonico = eta_y_phi(px_l, py_l, pz_l)

            if event.NJet != 0:

                for i in range(event.NJet):              # Bucle para determinar la distancia en el plano (eta,phi) entre los b-jets y los jets reconstruidos

                    eta_bjet, phi_bjet = eta_y_phi(event.Jet_Px[i], event.Jet_Py[i], event.Jet_Pz[i])
                    delta_eta_hadronico, delta_eta_leptonico = eta_bjet - eta_bjet_hadronico, eta_bjet - eta_bjet_leptonico
                    delta_phi_hadronico, delta_phi_leptonico = phi_bjet - phi_bjet_hadronico, phi_bjet - phi_bjet_leptonico

                    for j in [delta_phi_hadronico, delta_phi_leptonico]:               # Para que phi esté en [-pi, pi]
                        if j < - np.pi: j += 2 * np.pi
                        if j > np.pi: j -= 2 * np.pi

                    delta_R_hadronico = np.sqrt(delta_eta_hadronico ** 2 + delta_phi_hadronico ** 2)
                    delta_R_leptonico = np.sqrt(delta_eta_leptonico ** 2 + delta_phi_leptonico ** 2)

                    if delta_R_hadronico < 0.1 or delta_R_leptonico < 0.1:          # Criterio de cercanía en el plano
                        contador_bjets_buenos += 1
                        if event.Jet_btag[i] >= 2:
                            contador_btag += 1

        eficiencia_btagging = contador_btag / contador_bjets_buenos * 0.86          # Se multiplica además por el factor 0.86 especificado en el guión (para b-tagging >=2)


        # -----------------------------------------------------------------------------------------------------------------

        ### Cálculo de la eficiencia del trigger ###

        nbins = 250
        bin_final = 250

        h1 = r.TH1F("Pt_trig", "Momento Transverso (con trigger)", nbins, 0, bin_final)
        h2 = r.TH1F("Pt_notrig", "Momento Transverso (sin trigger)", nbins, 0, bin_final)
        h2.SetLineColor(r.kRed)

        t.Draw("sqrt(Muon_Px*Muon_Px + Muon_Py*Muon_Py)>>Pt_notrig")
        t.Draw("sqrt(Muon_Px*Muon_Px + Muon_Py*Muon_Py)>>Pt_trig", "triggerIsoMu24 == 1")

        h1.GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")
        h1.GetYaxis().SetTitle("Cuentas")
        h2.GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")
        h2.GetYaxis().SetTitle("Cuentas")

        r.gStyle.SetOptStat(0)

        canvas = r.TCanvas("Histogramas", "Histogramas de MET")

        canvas.Divide(1, 2)

        canvas.cd(1)
        h2.SetTitle("Momento Transverso (con trigger y sin trigger)")
        h2.SetTitleSize(0.1)
        h2.Draw("HIST")
        h1.Draw("HIST,SAME")
        leg = r.TLegend(0.7, 0.75, 0.5, 0.6)
        leg.SetTextSize(0.05)
        leg.AddEntry(h1, "Con trigger")
        leg.AddEntry(h2, "Sin trigger")
        leg.Draw()

        canvas.cd(2)
        ratio = h1.Clone("ratio")
        ratio.SetTitle("Cociente entre no trig y trig")
        ratio.SetTitleSize(0.1)
        ratio.SetLineColor(r.kBlack)
        ratio.Divide(h2)
        ratio.Draw("HIST")
        leg2 = r.TLegend(0.45, 0.3, 0.6, 0.2)
        leg2.SetTextSize(0.05)
        leg2.AddEntry(ratio, "Cociente")
        leg2.Draw()

        canvas.cd(1).SetLeftMargin(0.15)
        canvas.cd(2).SetLeftMargin(0.15)

        canvas.SaveAs("results/ComparacionCompleta_trigger.pdf")                 # Histogramas completos con trigger y sin trigger
        canvas.Update()

        nbins = 150
        bin_final = 150

        h3 = r.TH1F("Pt_trig_2", "Momento Transverso (con trigger)", nbins, 0, bin_final)
        h4 = r.TH1F("Pt_notrig_2", "Momento Transverso (sin trigger)", nbins, 0, bin_final)
        h4.SetLineColor(r.kRed)

        t.Draw("sqrt(Muon_Px*Muon_Px + Muon_Py*Muon_Py)>>Pt_notrig_2", "sqrt(Muon_Px*Muon_Px + Muon_Py*Muon_Py) >= 28")
        t.Draw("sqrt(Muon_Px*Muon_Px + Muon_Py*Muon_Py)>>Pt_trig_2", "triggerIsoMu24 == 1 && sqrt(Muon_Px*Muon_Px + Muon_Py*Muon_Py) >= 28")

        h3.GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")
        h3.GetYaxis().SetTitle("Cuentas")
        h4.GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")
        h4.GetYaxis().SetTitle("Cuentas")

        r.gStyle.SetOptStat(0)

        canvas2 = r.TCanvas("Histogramas2", "Histogramas de MET")

        canvas2.Divide(1, 2)

        canvas2.cd(1)
        h4.SetTitle("Momento Transverso (con trigger y sin trigger)")
        h4.SetTitleSize(0.8)
        h4.Draw("HIST")
        h3.Draw("HIST,SAME")
        leg3 = r.TLegend(0.7, 0.75, 0.5, 0.6)
        leg3.SetTextSize(0.05)
        leg3.AddEntry(h3, "Con trigger")
        leg3.AddEntry(h4, "Sin trigger")
        leg3.Draw()

        canvas2.cd(2)
        ratio2 = h3.Clone("ratio")
        ratio2.SetTitle("Cociente entre no trig y trig")
        ratio2.SetTitleSize(0.8)
        ratio2.SetLineColor(r.kBlack)
        ratio2.Divide(h4)
        ratio2.Draw("HIST")
        leg4 = r.TLegend(0.45, 0.3, 0.6, 0.2)
        leg4.SetTextSize(0.05)
        leg4.AddEntry(ratio2, "Cociente")
        leg4.Draw()

        canvas2.cd(1).SetLeftMargin(0.15)
        canvas2.cd(2).SetLeftMargin(0.15)

        canvas2.SaveAs("results/Eficiencia_trigger.pdf")                 # Histogramas en el rango adecuado de valores de pT con trigger y sin trigger
        canvas2.Update()

        eficiencia_trigger = ratio2.Integral() / nbins * bin_final / (bin_final - 28)    # Cálculo de la eficiencia del trigger

        if (isinstance(name, list)):
            for nam in name: self.XS(nam)
            return

        # variables para nº eventos por proceso (fondos y ttbar por separado) y nº eventos total
        backgrounds = ["wjets", "qcd", "ww", "wz", "zz", "dy", "single_top"]
        bkg = {i: 0. for i in backgrounds}
        ttbar_events = 0.
        total_bkg = 0.


        # variable para incertidumbre estadística (al cuadrado)

        incertidumbre_estadistica_cuadrado_total = 0.

        ### Eficiencia total y luminosidad ###

        eficiencia_muones = 0.99
        eficiencia = eficiencia_trigger * eficiencia_muones * eficiencia_btagging
        luminosidad = 50


        ### Cálculo del número de eventos de fondo y señal ###

        for s in self.listOfSelectors:
            h = s.GetHisto(name)
            if s.name != 'ttbar':
                total_bkg += h.Integral()
                bkg[s.name] += h.Integral()

            else: ttbar_events += h.Integral()


        print('\n---------------------------------------------------')
        print('Esperado (MC, no ttbar): {tot:.2f}\n'.format(tot = total_bkg))

        if self.data != '':
            hdata = self.dataSelector.GetHisto(name)
            datos = hdata.Integral()

        print('Número de eventos ttbar (MC): {ttbar:.2f}\n'.format(ttbar = ttbar_events))

        signal_events = datos - total_bkg # Resto a los eventos observados el fondo simulado

        print('Número de eventos ttbar /señal (observados): {ttbar:.2f}'.format(ttbar = signal_events))
        print('-----------------------------------------------------')

        xs = signal_events / (luminosidad * aceptancia * eficiencia)           # Valor nominal de la sección eficaz

        ### Cálculo de la Incertidumbre Estadística ###

        for s in self.listOfSelectors:
            incertidumbre_estadistica_cuadrado_total += s.incert_est_cuadrado

        if self.data != '':
            incertidumbre_estadistica_cuadrado_total += self.dataSelector.incert_est_cuadrado

        incertidumbre_estadistica_total = np.sqrt(incertidumbre_estadistica_cuadrado_total)       # Valor real de la incertidumbre estadística total

        xs_inc_est = incertidumbre_estadistica_total / (luminosidad * aceptancia * eficiencia)

        ### Cálculo de las Incertidumbres Sistemáticas ###

        ## Normalización de eventos

        # La diferencia en valor absoluto entre la sección eficaz nominal y la modificada es igual modificando el valor
        # del fondo sumando o restando el porcentaje correspondiente. Por lo tanto, sólo se calcula una contribución

        inc_bkg = [] # Calcula las incertidumbres para cada contribución de fondo modificado

        # porcentaje de incertidumbre en normalización (aparece en el guión)
        porcentaje = {"wjets": 0.5, "qcd": 1, "ww": 0.5, "wz": 0.5, "zz": 0.5, "dy": 0.15, "single_top": 0.3}

        for i in backgrounds:

            # incertidumbre en xs debida a la normalización
            xs_mod_bkg = (signal_events + bkg[i]*porcentaje[i]) / (luminosidad * aceptancia * eficiencia)
            inc_bkg.append(abs(xs_mod_bkg - xs))


        inc_bkg = np.array(inc_bkg)

        incertidumbre_norm = np.sqrt(sum(inc_bkg ** 2))        # Cálculo de la incertidumbre total por normalizacion como suma en cuadratura de las individuales

        ## Incertidumbres para aceptancia, eficiencias y luminosidad (fórmulas en el guión o deducidas en el informe)

        inc_acept = np.sqrt(aceptancia*(1 - aceptancia)/contador_total_ttbar)
        inc_btagging = np.sqrt(eficiencia_btagging*(1 - eficiencia_btagging)/contador_bjets_buenos)
        inc_muones = 0.01
        inc_trig = (1 - eficiencia_trigger)/2
        inc_lumi = luminosidad*0.1

        denominador = np.array([aceptancia, eficiencia_btagging, eficiencia_muones, eficiencia_trigger, luminosidad])     # Ya que están todas en el denominador
        incertidumbres_denominador = np.array([inc_acept, inc_btagging, inc_muones, inc_trig, inc_lumi])

        xs_mod_suma = xs*denominador/(denominador + incertidumbres_denominador)
        xs_mod_resta = xs*denominador/(denominador - incertidumbres_denominador)

        inc_xs_suma = abs(xs_mod_suma - xs)               # incertidumbre en xs debida a las magnitudes del denominador
        inc_xs_resta = abs(xs_mod_resta - xs)

        # Incertidumbres sistemáticas completas variando arriba y abajo, sin incluir luminosidad

        incertidumbre_sist_total_suma = np.sqrt(incertidumbre_norm**2 + sum(inc_xs_suma[:-1]**2))
        incertidumbre_sist_total_resta = np.sqrt(incertidumbre_norm**2 + sum(inc_xs_resta[:-1]**2))


        ### Prints finales ###

        encabezado_tabla_normalizacion = ["Fondo", "Incertidumbre", "Propagación a xs (pb) (+,-)"]                # Tabla con valores de la incertidumbre por normalización
        datos_tabla_normalizacion = [
            ["W + jets", porcentaje["wjets"]*bkg["wjets"], inc_bkg[0]],
            ["QCD", porcentaje["qcd"]*bkg["qcd"], inc_bkg[1]],
            ["WW", porcentaje["ww"]*bkg["ww"], inc_bkg[2]],
            ["WZ", porcentaje["wz"]*bkg["wz"], inc_bkg[3]],
            ["ZZ", porcentaje["zz"]*bkg["zz"], inc_bkg[4]],
            ["DY", porcentaje["dy"]*bkg["dy"], inc_bkg[5]],
            ["Single top", porcentaje["single_top"]*bkg["single_top"], inc_bkg[6]],
            ["Norm_total", incertidumbre_norm*(luminosidad * aceptancia * eficiencia), incertidumbre_norm]
        ]


        encabezado_tabla_denom = ["", "Valor", "Incertidumbre", "Propagación a xs (pb) (+)", "Propagación a xs (pb) (-)"]   # Tabla con valores de las incertidumbres sistemáticas restantes (en el denom)
        datos_tabla_denom = [
            ["Aceptancia", denominador[0], incertidumbres_denominador[0], inc_xs_suma[0], inc_xs_resta[0]],
            [f"\u03b5 (b-tagging)", denominador[1], incertidumbres_denominador[1], inc_xs_suma[1], inc_xs_resta[1]],
            [f"\u03b5 (muones)", denominador[2], incertidumbres_denominador[2], inc_xs_suma[2], inc_xs_resta[2]],
            [f"\u03b5 (trigger)", denominador[3], incertidumbres_denominador[3], inc_xs_suma[3], inc_xs_resta[3]],
            ["Luminosidad", denominador[4], incertidumbres_denominador[4], inc_xs_suma[4], inc_xs_resta[4]]
        ]



        formato_encabezado_norm = "| {:<15} | {:<15} | {:<28} |"
        formato_fila_norm = "| {:<15} | {:>15.4f} | {:>28.4f} |"

        formato_encabezado_denom = "| {:<15} | {:<15} | {:<15} | {:<25} | {:<25} |"
        formato_fila_denom = "| {:<15} | {:>15.4f} | {:>15.4f} | {:>25.4f} | {:>25.4f} |"
        separador = "-" * 68


        print("\nIncertidumbres Sistemáticas (Normalización):" + "\n" + separador)
        print(formato_encabezado_norm.format(*encabezado_tabla_normalizacion))
        print(separador)

        for row in datos_tabla_normalizacion:
            print(formato_fila_norm.format(*row))

        print(separador + "\n")

        separador = "-" * 111

        print("Resto de Incertidumbres Sistemáticas:" + "\n" + separador)
        print(formato_encabezado_denom.format(*encabezado_tabla_denom))
        print(separador)

        for row in datos_tabla_denom:
            print(formato_fila_denom.format(*row))

        print(separador + "\n")

        # Print de las incertidumbres estadísticas totales, y valor nominal de la sección eficaz con su incertidumbre asociada a cada fuente

        print('Incertidumbre estadística total:   {inc_est:.2f}\n'.format(inc_est = incertidumbre_estadistica_total))
        print('Incertidumbre sistemática total (sin luminosidad):   +{inc_sist_suma:.2f}, - {inc_sist_resta:.2f}\n'.format(inc_sist_suma = incertidumbre_sist_total_suma, inc_sist_resta = incertidumbre_sist_total_resta))
        print(f'Se ha obtenido una sección eficaz de {xs:.2f} \u00b1 {xs_inc_est:.2f} (stat) + {incertidumbre_sist_total_suma:.2f} - {incertidumbre_sist_total_resta:.2f} (sist) +'
              f'{inc_xs_suma[4]:.2f} - {inc_xs_resta[4]:.2f} (lumi) pb\n')

        return


    def SaveCounts(self, name, overridename = ""):
        ''' Save in a text file the number of events for each sample in a given histogram '''
        if (isinstance(name, list)):
            for nam in name: self.SaveCounts(nam, overridename = overridename)
            return

        filename = "yields_{h}".format(h = name) if (overridename == "") else overridename
        create_folder(self.savepath)
        outfile = open(self.savepath + "/" + filename + (".txt" if ".txt" not in overridename else ""), "w")

        thelines = []

        thelines.append("Number of events for histogram {h}:\n".format(h = name))
        thelines.append("----------------------------------------------------\n")
        total = 0.
        for s in self.listOfSelectors:
            h = s.GetHisto(name)
            thelines.append("{nam}: {num}\n".format(nam = s.name, num = h.Integral()))
            total += h.Integral()
        thelines.append('Expected (MC): {tot}\n'.format(tot = total))
        thelines.append('------------------------------\n')
        if self.data != '':
            hdata = self.dataSelector.GetHisto(name)
            thelines.append('Observed (data): {tot}\n'.format(tot = hdata.Integral()))
            thelines.append('------------------------------\n')

        outfile.writelines(thelines)
        outfile.close()
        return


def create_folder(path):
    if not os.path.exists(path): os.system("mkdir -p " + path)
    return
