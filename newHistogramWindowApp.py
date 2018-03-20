
###################################################################################################
# New Histogram window Application
# Within the main graph window, the user has an option to create a Mass Correlation Histogram 
# When the user clicks 'Create Histogram' in Main graph window, a popup window to load 
# the datasets needed appears.  The user selects the appropriate data files, (raw data from the atom
# probing experiments) and then clicks 'create graph' button 
# This calls the constructore for the new histogram window app with the selected data - resulting is
# a window with an image of the graph and the various graphing and plotting funtions
####################################################################################################


# -*- coding: utf-8 -*-
from PyQt4 import QtGui, QtCore  # Import the PyQt4 module we'll need
import PyQt4.uic
import sys  # We need sys so that we can pass argv to QApplication
import os
import numpy as np
import cPickle as pickle
from datetime import datetime # allows time stamping
import itertools
import re
import collections
import CorrHistDesign

from CorrHistObj import *

import matplotlib
matplotlib.rcParams['backend.qt4'] = 'PyQt4'  # Make sure the correct qt backend is set (not PySide)
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

# Class to define the correlation histogram window
# an instance of this window is called in the maingraph window when the launch button
# is pressed

class NewHistogramWindowApp(QtGui.QDialog,CorrHistDesign.Ui_Dialog):
    '''this class defines a new window object which displays the correlation histogram.
        it contains the attributes of the corrHistObj created, and contains the functionality to
        create the correlation histogram and display it in a window'''
        
    # Constructor
    # should this inherit the methods from mainAPTApp or MainGraphWindow??
    def __init__(self,mainGraphWindowObj):
        
        super(self.__class__, self).__init__()
        self.setupUi(self)  # This is defined in design.py file automatically
        
        # set Histogram window attribute array to data object 
        self.mainGraphWindowObj = mainGraphWindowObj
        self.pyTableFileName = self.mainGraphWindowObj.mainAPTAppObj.pyTableFileName
        
        # Create a Corr Hist Object instance from .h5 file
        #self.aCorrHistObj = mainGraphWindowObj.aCorrHistObj
        #edit to take in a filename determined in create_histogram_window
        self.aCorrHistObj = mainGraphWindowObj.mainAPTAppObj.mainGraphWindow.aCorrHistObj
        
        # Set window attributes to object attributes
        self.nBins = self.aCorrHistObj.nBins #number of bins on histogram
        self.maxDa = self.aCorrHistObj.maxDa #max Da displayed on histogram
        self.bwCorrHist = self.aCorrHistObj.bwCorrHist
        #self.binEdgeArray = self.aCorrHistObj.binEdgeArray
        self.f_mtc = self.aCorrHistObj.f_mtc #first of multiples array, repeats included
        self.s_mtc = self.aCorrHistObj.s_mtc #all following hits in a given event
        self.H = self.aCorrHistObj.H # raw count
        self.h_log = self.aCorrHistObj.h_log # log count
        
    
        
        # Add the canvas for the graph and the toolbar
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.widgetGraph)
        self.axes = self.fig.add_subplot(111)
        
        

        vbox = self.verticalLayoutGraph
        vbox.addWidget(self.canvas)  # the matplotlib canvas
        vbox.addWidget(self.mpl_toolbar)
        
            
        # Check default boxes to display initial graph
        self.LogRadioButton.setChecked(True)

        self.LogRadioButton.clicked.connect(self.replotgraph) #connection to log color scale
        self.BinRadioButton_2.clicked.connect(self.replotgraph) #connection to binary color scale
        self.DissTracksCheckBox.clicked.connect(self.replotgraph) #enables or disables the feature to add dissociation tracks
        #self.DissTracksCheckBox.clicked.connect(self.populate_peak_list)
        self.doubleSpinBox.valueChanged.connect(self.change_bw) #connection to the function which controls slider binwidth and replots the graph
        self.listWidget_2.itemClicked.connect(self.choose_species) # connection from first peak list to method that populates species table
        self.listWidget_3.itemClicked.connect(self.choose_first_daughter) #connections from listWidget_3 to lineEdit_6
        self.listWidget.itemClicked.connect(self.choose_species_2) #connection from listWidget_2 to listWidget_4
        self.listWidget_4.itemClicked.connect(self.choose_second_daughter) #connection from listWidget_4 (possible species for y) to lineEdit_7
        self.checkBox_OnOff.clicked.connect(self.replotgraph)
        
        #connections to sorting radio buttons
        self.radioButton.clicked.connect(self.sort_atomlist_1)
        self.radioButton_2.clicked.connect(self.sort_atomlist_1)
        self.radioButton_6.clicked.connect(self.sort_atomlist_1)
        self.radioButton_9.clicked.connect(self.sort_atomlist_2)
        self.radioButton_10.clicked.connect(self.sort_atomlist_2)
        self.radioButton_12.clicked.connect(self.sort_atomlist_2)
        self.radioButton_5.clicked.connect(self.sort_atomlist_1)
        self.radioButton_11.clicked.connect(self.sort_atomlist_2)
        self.radioButton_3.clicked.connect(self.sort_atomlist_1)
        self.radioButton_4.clicked.connect(self.sort_atomlist_1)
        self.radioButton_7.clicked.connect(self.sort_atomlist_2)
        self.radioButton_8.clicked.connect(self.sort_atomlist_2)
        self.im = self.axes.imshow(self.h_log, origin='low', cmap = 'jet')
        #self.h2, = self.axes.plot(np.random.randn(50))
        self.updatestatus()
        #self.replotgraph()
    
    def updatestatus(self):
        # If dissociation tracks are enabled, enable the group box with features
        # -- populate peak lists (listWidget_1, listWidget_2) with data from peakList window
        # else, reset all lines in the groupbox,and disable the groupbox
        if self.DissTracksCheckBox.isChecked():
            self.groupBox_4.setEnabled(True)
            if self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet is not None:
                temp_list = self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet.peakListSorted
                for i in temp_list:
                    da = i.peakLabel
                    self.listWidget.addItem(da)
                    self.listWidget_2.addItem(da)
            
        else:
            self.listWidget.clear()
            self.listWidget_2.clear()
            self.listWidget_3.clear()
            self.listWidget_4.clear()
            self.lineEdit_6.clear()
            self.lineEdit_7.clear()
            self.groupBox_4.setEnabled(False)         
            
        
    def choose_species(self):
        self.replotgraph()
        self.groupBox.setEnabled(True)
        self.listWidget_3.clear()
        self.lineEdit_6.clear()
        self.lineEdit.clear()
        self.radioButton.setChecked(True)
        aPeakLabel = str(self.listWidget_2.currentItem().text())
        # self.projected_peaks = self.mainGraphWindowObj.mainAPTAppObj.peakInfoWindow.currentPossibleSpecies.projectedPeaksList
        # self.currentIdx = self.listWidget_2.currentItem().row()
        peakInfoWindow = self.mainGraphWindowObj.mainAPTAppObj.peakInfoWindow
        aRangedPeakSet = self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet
        aPeak = (aRangedPeakSet.peakMTCDict[aPeakLabel])
        #peakInfoWindow.updateData(aPeak)   
        self.projected_peaks = aPeak.projectedPeaksList
        
        list_of_lists = []
        for i in range(len(self.projected_peaks)):
            # molecule = str(self.projected_peaks[i].peakIDObj.peakIDStr)
            #self.m_1 = (self.projected_peaks[i].mtc)
            # self.listWidget_3.addItem(molecule) 
            
            if float(self.projected_peaks[i].abundance) >= 0.04:
                molecule = str(self.projected_peaks[i].peakIDObj.peakIDStr)
                #self.listWidget_3.addItem(molecule) 
                molecule_list = self.projected_peaks[i].molIsotopeList
                list_of_lists.append((molecule_list[0].isotopeAtomsList , molecule))  
                    # for mol in molecule_list:
                    #     list_of_lists.append((mol.isotopeAtomsList,molecule))
                    #     print '--',mol.molecule.formula, mol.isotopeAtomsList
            
                print
                print "NEW MOLECULE: ", self.projected_peaks[i].peakIDObj.molecule.formula, "mtc: ", self.projected_peaks[i].mtc, "abundance: ", self.projected_peaks[i].abundance
                print "a.ions dictionary: ",self.projected_peaks[i].peakIDObj.molecule.ionsDict
                print
                print "b.atoms list: "
                for mol in self.projected_peaks[i].molIsotopeList:
                    print '--',mol.molecule.formula, mol.isotopeAtomsList
                print
                elements_list =  self.projected_peaks[i].peakIDObj.molecule.ionsDict.keys()
                print "c.element isotope values:"
                for element in elements_list:
                    if element in self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet.allElementsDict: 
                        print '--',element, self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet.allElementsDict[element].isotopes
                # list_item = str("%.3f" % float(self.projected_peaks[i].mtc))+ ' (Da)  '+str("%.3f" % float(self.projected_peaks[i].abundance*100.0))+' % (Abund.)'
                # self.listWidget_3.addItem(list_item)
        list_of_lists.sort(key =lambda x: len(x[0]))
        print list_of_lists
        self.listWidget_3.clear()
        for item in list_of_lists:
            self.listWidget_3.addItem(str(item[1]))
                
                
    
    def choose_species_2(self):
        self.replotgraph()
        self.groupBox_2.setEnabled(True)
        self.listWidget_4.clear()
        self.lineEdit_7.clear()
        self.lineEdit.clear()
        self.radioButton_9.setChecked(True)
        aPeakLabel = str(self.listWidget.currentItem().text())
        # self.projected_peaks = self.mainGraphWindowObj.mainAPTAppObj.peakInfoWindow.currentPossibleSpecies.projectedPeaksList
        # self.currentIdx = self.listWidget_2.currentItem().row()
        peakInfoWindow = self.mainGraphWindowObj.mainAPTAppObj.peakInfoWindow
        aRangedPeakSet = self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet
        aPeak = (aRangedPeakSet.peakMTCDict[aPeakLabel])
        #peakInfoWindow.updateData(aPeak)   
        self.projected_peaks_2 = aPeak.projectedPeaksList
        
        list_of_lists = []
        for i in range(len(self.projected_peaks_2)):
              
            if float(self.projected_peaks_2[i].abundance) >= 0.04:
                molecule = str(self.projected_peaks_2[i].peakIDObj.peakIDStr)
                molecule_list = self.projected_peaks_2[i].molIsotopeList
                list_of_lists.append((molecule_list[0].isotopeAtomsList , molecule)) 
            #self.m_1 = (self.projected_peaks[i].mtc)
                #self.listWidget_4.addItem(molecule) 
                print
                print "NEW MOLECULE: ", self.projected_peaks_2[i].peakIDObj.molecule.formula, "mtc: ", self.projected_peaks_2[i].mtc, "abundance: ", self.projected_peaks_2[i].abundance
                print "a.ions dictionary: ",self.projected_peaks_2[i].peakIDObj.molecule.ionsDict
                print
                print "b.atoms list: "
                for mol in self.projected_peaks_2[i].molIsotopeList:
                    print '--',mol.molecule.formula, mol.isotopeAtomsList
                print
                elements_list =  self.projected_peaks_2[i].peakIDObj.molecule.ionsDict.keys()
                print "c.element isotope values:"
                for element in elements_list:
                    if element in self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet.allElementsDict: 
                        print '--',element, self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet.allElementsDict[element].isotopes
                # list_item = str("%.3f" % float(self.projected_peaks[i].mtc))+ ' (Da)  '+str("%.3f" % float(self.projected_peaks[i].abundance*100.0))+' % (Abund.)'
                # self.listWidget_3.addItem(list_item)
        list_of_lists.sort(key =lambda x: len(x[0]))
        print list_of_lists
        self.listWidget_4.clear()
        for item in list_of_lists:
            self.listWidget_4.addItem(str(item[1]))
                            
                
                
    def choose_first_daughter(self):
        self.lineEdit_6.clear()
        molecule = str(self.listWidget_3.currentItem().text())
        for i in range(len(self.projected_peaks)):
            if self.projected_peaks[i].peakIDObj.peakIDStr == molecule:
                print "MOLECULE: ",molecule, ', MTC: ',self.projected_peaks[i].mtc
                self.molecule_1 = self.projected_peaks[i].peakIDObj.peakIDStr
                self.m_1 = self.projected_peaks[i].mtc 
                self.lineEdit_6.setText(str("%.3f" % float(self.m_1)))
                break
        self.replotgraph()
                
    def choose_second_daughter(self):
        self.lineEdit_7.clear()
        molecule_2 = str(self.listWidget_4.currentItem().text())
        for i in range(len(self.projected_peaks_2)):
            if self.projected_peaks_2[i].peakIDObj.peakIDStr == molecule_2:
                print "MOLECULE: ",molecule_2, ', MTC: ',self.projected_peaks_2[i].mtc
                self.molecule_2 = self.projected_peaks_2[i].peakIDObj.peakIDStr
                self.m_2 = self.projected_peaks_2[i].mtc 
                self.lineEdit_7.setText(str("%.3f" % float(self.m_2)))
                break
        self.replotgraph()
                    
    def fill_atoms_list(self): 
        self.lineEdit_3.clear()
        self.lineEdit_4.clear()
        self.listWidget.clear()
        self.listWidget_2.clear()
        i = self.listWidget_1.currentRow()
        self.m_p = float(self.projected_peaks[i].mtc)
        self.replotgraph()
        print
        print "atoms list for: ",self.projected_peaks[i].peakIDObj.molecule.formula
        for isotope in self.projected_peaks[i].molIsotopeList:
            print '--',isotope.molecule.formula, isotope.isotopeAtomsList
            self.listWidget_2.addItem(str(isotope.isotopeAtomsList))
    
    def fill_dissociations_list(self):
        self.listWidget.clear()
        self.lineEdit_4.clear()
        self.lineEdit_3.clear()
        i = self.listWidget_1.currentRow()
        j = self.listWidget_2.currentRow()
        self.current_atoms_list = self.projected_peaks[i].molIsotopeList[j].isotopeAtomsList
        self.find_all_reactions(self.current_atoms_list)
        
    def assign_daughter_ions(self):
        k = self.listWidget.currentRow()
        self.reaction = self.reactions_list[k]
        molecule = 0
        molecule_2 = 0
        self.m_1 = 0
        self.m_2 = 0
        print
        print "dissociation reaction: "
        print self.reaction[0],'+',self.reaction[1]
        for atom in self.reaction[0]:
            element, iso_num = atom.split("_")
            atom_mtc = self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet.allElementsDict[element].isotopes[int(iso_num)]
            #print "element+ isotope: ", element, iso_num, "mass-to-charge: ", mtc[0] 
            molecule += atom_mtc[0]
            
        mol_formula_1 = self.chem_formula(self.reaction[0])
            
        for atom in self.reaction[1]:
            element, iso_num = atom.split("_")
            atom_mtc = self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet.allElementsDict[element].isotopes[int(iso_num)]
            #print "element+ isotope: ", element, iso_num, "mass-to-charge: ", mtc[0] 
            molecule_2 += atom_mtc[0]
        mol_formula_2 = self.chem_formula(self.reaction[1])
            
        if molecule >= molecule_2:
            self.m_1 = molecule
            self.m_2 = molecule_2
            self.form_1 = mol_formula_1
            self.form_2 = mol_formula_2
        else:
            self.m_1 = molecule_2
            self.m_2 = molecule
            self.form_1 = mol_formula_2
            self.form_2 = mol_formula_1

               
        self.lineEdit_3.setText(self.form_1+", "+str("%.3f" % float(self.m_1)))
        self.lineEdit_4.setText(self.form_2+", "+str("%.3f" % float(self.m_2)))
            
        self.replotgraph()
        
    def chem_formula(self,atom_list):
        formula = []
        for i in atom_list:
            atom, iso_num = i.split("_")
            formula.append(atom)
        #print formula
        mol = [(i,str(formula.count(i))) for i in set(formula)]
        for y in range(len(mol)):
            if mol[y][1] == '1':
                mol[y] = (mol[y][0],'')
        self.mol_formula = ''.join(map(''.join,mol))
        return self.mol_formula
                 
                
    def change_bw(self):
        da =  float(self.doubleSpinBox.value())
        self.bwCorrHist = da
        num_bins = int(np.round(self.maxDa/da))
        
        # remake histogram everytime the value is changed
        # here use the import new da to adjust the number of bins and plot accordingly
        self.H, self.xedges, self.yedges = np.histogram2d(self.f_mtc, self.s_mtc, bins=(num_bins,num_bins), range = [[0,125],[0,125]])
        self.h_log = np.log(self.H +.1)
        self.replotgraph()
        
        
    def replotgraph(self): 
              
        # Create the histogram
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        
        #potentially add a conditional for the diffstrackscheckbox here....
        self.updatestatus()
        temp_table = self.mainGraphWindowObj.mainAPTAppObj.peakListWindow.tablePeakSummary
        current_index = self.mainGraphWindowObj.mainAPTAppObj.peakListWindow.currentIdx
        

        # if bool(self.listWidget_1.item(0)) == True:
        #     print self.m_p
        #     self.axes.plot(self.m_p/self.bwCorrHist,self.m_p/self.bwCorrHist, 'ro', markersize = 2)
        #     
        # if bool(self.listWidget.item(0)) == True:
        #     print self.m_p,self.m_1,self.m_2
        #     t = range(0,5000)
        #     x = [(self.m_1*(1-((i/5000.)*(1.-(self.m_1/self.m_p))))**-1.)/self.bwCorrHist for i in t]
        #     y = [(self.m_2*(1-((i/5000.)*(1.-(self.m_2/self.m_p))))**-1.)/self.bwCorrHist for i in t]
        #     self.axes.plot(x,y,'r')

        if self.LogRadioButton.isChecked():
            self.im = self.axes.imshow(self.h_log, origin='low', cmap = 'jet')
            #self.im.set_data(self.h_log)
            #self.fig.colorbar(self.im)
            self.xScaleFactor =self.bwCorrHist
            self.yScaleFactor = self.xScaleFactor
            self.ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*self.xScaleFactor))
            self.ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*self.yScaleFactor))
            self.axes.xaxis.set_major_formatter(self.ticks_x)
            self.axes.yaxis.set_major_formatter(self.ticks_y)
            
        else:  
            norm = matplotlib.colors.Normalize(vmin=0, vmax=1)   
            self.im = self.axes.imshow(self.H, origin='low', cmap = 'hot',norm = norm)
            #self.im.set_data(self.H)
            #self.fig.colorbar(self.im)
            self.xScaleFactor =self.bwCorrHist
            self.yScaleFactor = self.xScaleFactor
            self.ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*self.xScaleFactor))
            self.ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*self.yScaleFactor))
            self.axes.xaxis.set_major_formatter(self.ticks_x)
            self.axes.yaxis.set_major_formatter(self.ticks_y)
            
        if (bool(self.lineEdit_6.text()) == True) and (bool(self.lineEdit_7.text()) == True) and self.checkBox_OnOff.isChecked():
            print "graph point ", self.m_1, self.m_2
            
            self.dissociation_equation()
            self.m_p = (self.m_1 + self.m_2)/self.charge_state
            print "mtc of parent ion: ", self.m_p
            t = range(0,5000)
            x = [(self.m_1*(1-((i/5000.)*(1.-(self.m_1/self.m_p))))**-1.)/self.bwCorrHist for i in t]
            y = [(self.m_2*(1-((i/5000.)*(1.-(self.m_2/self.m_p))))**-1.)/self.bwCorrHist for i in t]
            
            self.axes.plot(self.m_1/self.bwCorrHist,self.m_2/self.bwCorrHist, 'ro', markersize = 2)
            self.axes.plot(x,y,'r')
            

        self.canvas.draw()
        
    def dissociation_equation(self):
        mol_1 = str(self.listWidget_3.currentItem().text())
        mol_2 = str(self.listWidget_4.currentItem().text())
        print self.lineEdit_6.text()
        print self.lineEdit_7.text()
        
        self.charge_state = 0
        self.charge_state += int(mol_1.split('+')[1])
        self.charge_state += int(mol_2.split('+')[1])
        split_molecule =(re.findall(r'([A-Z][a-z]*)(\d*)', mol_1))
        split_mol_2 = (re.findall(r'([A-Z][a-z]*)(\d*)', mol_2))
        
        molecule_list = []
        molecule_string = ''
        for tuple in range(len(split_molecule)):
            if split_molecule[tuple][1] == '':
                split_molecule[tuple] = (split_molecule[tuple][0],'1')
            piece = ''
            for count in range(int(split_molecule[tuple][1])):
                piece += split_molecule[tuple][0]
                molecule_string += split_molecule[tuple][0]
                molecule_list.append(split_molecule[tuple][0])
        for tuple in range(len(split_mol_2)):
            if split_mol_2[tuple][1] == '':
                split_mol_2[tuple] = (split_mol_2[tuple][0],'1')
            piece = ''
            for count in range(int(split_mol_2[tuple][1])):
                piece += split_mol_2[tuple][0]
                molecule_string += split_mol_2[tuple][0]
                molecule_list.append(split_mol_2[tuple][0])
        
        f_mol = [(i,str(molecule_list.count(i))) for i in set(molecule_list)]
        for y in range(len(f_mol)):
                if f_mol[y][1] == '1':
                    f_mol[y] = (f_mol[y][0],'')
        formula = ''.join(map(''.join,f_mol)) 
        print "new molecule: ", formula
        print "Charge State:", self.charge_state
        
        dissociation_equation = formula+' --> '+mol_1+' + '+mol_2
        self.lineEdit.setText(dissociation_equation) 
        
        
    def find_all_reactions(self,molecule_list):
        reactions = []  
        flags = [False] * len(molecule_list)
        while True:
            a = [molecule_list[i] for i, flag in enumerate(flags) if flag]
            b = [molecule_list[i] for i, flag in enumerate(flags) if not flag]
            if a != [] and b != []:
                reactions.append([a,b])
            for i in xrange(len(molecule_list)):
                flags[i] = not flags[i]
                if flags[i]:
                    break
            else:
                break
        reactions.sort()
    #     remove_duplicates = list(sublists for sublists,_ in itertools.groupby(sublists))
        self.reactions_list = []
        for item in reactions:
            if sorted(item) not in self.reactions_list:
                self.reactions_list.append(sorted(item))
        for item in self.reactions_list:
            f_mol = self.chem_formula(item[0])
            s_mol = self.chem_formula(item[1])
            self.listWidget.addItem(f_mol+'+'+s_mol)
            #self.listWidget.addItem(str(item))
            
    def sort_atomlist_1(self):
        aPeakLabel = str(self.listWidget_2.currentItem().text())
        peakInfoWindow = self.mainGraphWindowObj.mainAPTAppObj.peakInfoWindow
        aRangedPeakSet = self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet
        aPeak = (aRangedPeakSet.peakMTCDict[aPeakLabel])
        #self.lineEdit_2.clear()
        
        #Sort by number of atoms
        if self.radioButton.isChecked():
            projPeakList = sorted(self.projected_peaks,key=lambda x: x.peakIDObj.molecule.numAtoms)
            self.listWidget_3.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_3.addItem(str(item.peakIDObj.peakIDStr))
                
        # sort by number of molecule species
        elif self.radioButton_2.isChecked():
            projPeakList = sorted(self.projected_peaks,key=lambda x: len(x.peakIDObj.molecule.ionsDict.keys()))
            self.listWidget_3.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_3.addItem(str(item.peakIDObj.peakIDStr))
                
        #sort by mtc distance  
        elif self.radioButton_6.isChecked():
            projPeakList = sorted(self.projected_peaks,key=lambda x: abs(x.mtc-aPeak.centroidMTC))
            self.listWidget_3.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_3.addItem(str(item.peakIDObj.peakIDStr))
                
        #sort by number of peak matches        
        elif self.radioButton_5.isChecked():
            projPeakList = sorted(self.projected_peaks,key=lambda x: x.peakIDObj.numPosMatchPeaks,reverse=True)
            self.listWidget_3.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_3.addItem(str(item.peakIDObj.peakIDStr))
        
        # sort by number of speficied element
        elif self.radioButton_3.isChecked():
            oneElementStr = str(self.lineEdit_2.text())
            if len(oneElementStr) > 0:
                if len(oneElementStr.split(',')) > 1:
                    infoStr = "Please enter only one element symbol."
                    QtGui.QMessageBox.information(self, "Ok", infoStr)
                    self.radioButton.setChecked(True)
                    self.sort_atomlist_1()
                elif oneElementStr not in aRangedPeakSet.allElementsDict.keys():
                    infoStr = "Input name '"+oneElementStr+"' not found in element symbol dictionary."
                    QtGui.QMessageBox.information(self, "Ok", infoStr)
                    self.radioButton.setChecked(True)
                    self.sort_atomlist_1()
                else:
                    projPeakList = sorted(self.projected_peaks,key=lambda x: x.peakIDObj.getNumSpecificElement(oneElementStr),reverse=True)
                    self.listWidget_3.clear()
                    for item in projPeakList:
                        if item.abundance >= 0.04:
                            self.listWidget_3.addItem(str(item.peakIDObj.peakIDStr))
            else:
                infoStr = "Please enter an element symbol."
                QtGui.QMessageBox.information(self, "Ok", infoStr)
                self.radioButton.setChecked(True)
                self.sort_atomlist_1()
        
        # sort by element match?
        elif self.radioButton_4.isChecked():
            multiElementStr = str(self.lineEdit_5.text())
            if len(multiElementStr) > 0:
                errorFlag = False
                multiElementStrList = multiElementStr.split(',')
                print multiElementStrList
                for oneElementStr in multiElementStrList:
                    if oneElementStr not in aRangedPeakSet.allElementsDict.keys():
                        infoStr = "Input name '"+oneElementStr+"' not found in element symbol dictionary."
                        QtGui.QMessageBox.information(self, "Ok", infoStr)
                        self.radioButton.setChecked(True)
                        errorFlag = True
                        break
                if errorFlag == False:
                    projPeakList = sorted(self.projected_peaks,key=lambda x: x.peakIDObj.calculateMatchWithElementList(multiElementStrList),reverse=True)
                    self.listWidget_3.clear()
                    for item in projPeakList:
                        if item.abundance >= 0.04:
                            self.listWidget_3.addItem(str(item.peakIDObj.peakIDStr))
            else:
                infoStr = "Please enter at least one element symbol."
                QtGui.QMessageBox.information(self, "Ok", infoStr)
                self.radioButton.setChecked(True)
                self.sort_atomlist_1()     
            
            
    def sort_atomlist_2(self):
        
        aPeakLabel = str(self.listWidget.currentItem().text())
        peakInfoWindow = self.mainGraphWindowObj.mainAPTAppObj.peakInfoWindow
        aRangedPeakSet = self.mainGraphWindowObj.mainAPTAppObj.rangedPeakSet
        aPeak = (aRangedPeakSet.peakMTCDict[aPeakLabel])
        
        #sort by number of atoms
        if self.radioButton_9.isChecked():
            projPeakList = sorted(self.projected_peaks_2,key=lambda x: x.peakIDObj.molecule.numAtoms)
            self.listWidget_4.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_4.addItem(str(item.peakIDObj.peakIDStr))
                
        #sort by number of molecule species        
        elif self.radioButton_10.isChecked():
            projPeakList = sorted(self.projected_peaks_2,key=lambda x: len(x.peakIDObj.molecule.ionsDict.keys()))
            self.listWidget_4.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_4.addItem(str(item.peakIDObj.peakIDStr))
         
        #sort by mtc distance               
        elif self.radioButton_12.isChecked():
            projPeakList = sorted(self.projected_peaks_2,key=lambda x: abs(x.mtc-aPeak.centroidMTC))
            self.listWidget_4.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_4.addItem(str(item.peakIDObj.peakIDStr))
             
        #sort by peak matches           
        elif self.radioButton_11.isChecked():
            projPeakList = sorted(self.projected_peaks_2,key=lambda x: x.peakIDObj.numPosMatchPeaks,reverse=True)
            self.listWidget_4.clear()
            for item in projPeakList:
                if item.abundance >= 0.04:
                    self.listWidget_4.addItem(str(item.peakIDObj.peakIDStr))
        
        # sort by number of speficied element
        elif self.radioButton_8.isChecked():
            oneElementStr = str(self.lineEdit_3.text())
            if len(oneElementStr) > 0:
                if len(oneElementStr.split(',')) > 1:
                    infoStr = "Please enter only one element symbol."
                    QtGui.QMessageBox.information(self, "Ok", infoStr)
                    self.radioButton_9.setChecked(True)
                    self.sort_atomlist_2()
                elif oneElementStr not in aRangedPeakSet.allElementsDict.keys():
                    infoStr = "Input name '"+oneElementStr+"' not found in element symbol dictionary."
                    QtGui.QMessageBox.information(self, "Ok", infoStr)
                    self.radioButton_9.setChecked(True)
                    self.sort_atomlist_2()
                else:
                    projPeakList = sorted(self.projected_peaks_2,key=lambda x: x.peakIDObj.getNumSpecificElement(oneElementStr),reverse=True)
                    self.listWidget_4.clear()
                    for item in projPeakList:
                        if item.abundance >= 0.04:
                            self.listWidget_4.addItem(str(item.peakIDObj.peakIDStr))
            else:
                infoStr = "Please enter an element symbol."
                QtGui.QMessageBox.information(self, "Ok", infoStr)
                self.radioButton_9.setChecked(True)
                self.sort_atomlist_2()
        
        # sort by element match?
        elif self.radioButton_7.isChecked():
            multiElementStr = str(self.lineEdit_4.text())
            if len(multiElementStr) > 0:
                errorFlag = False
                multiElementStrList = multiElementStr.split(',')
                print multiElementStrList
                for oneElementStr in multiElementStrList:
                    if oneElementStr not in aRangedPeakSet.allElementsDict.keys():
                        infoStr = "Input name '"+oneElementStr+"' not found in element symbol dictionary."
                        QtGui.QMessageBox.information(self, "Ok", infoStr)
                        self.radioButton_9.setChecked(True)
                        errorFlag = True
                        break
                if errorFlag == False:
                    projPeakList = sorted(self.projected_peaks_2,key=lambda x: x.peakIDObj.calculateMatchWithElementList(multiElementStrList),reverse=True)
                    self.listWidget_4.clear()
                    for item in projPeakList:
                        if item.abundance >= 0.04:
                            self.listWidget_4.addItem(str(item.peakIDObj.peakIDStr))
            else:
                infoStr = "Please enter at least one element symbol."
                QtGui.QMessageBox.information(self, "Ok", infoStr)
                self.radioButton_9.setChecked(True)
                self.sort_atomlist_2()
            


################################################################################
