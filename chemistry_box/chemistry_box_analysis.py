#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.colors as plc
from matplotlib.animation import FuncAnimation
import tools
import numpy as np
from cycler import cycler
plt.style.use("default")
plt.rcParams.update({'font.size': 12})
plt.rc("axes", axisbelow = True) #Makes sure gridlines are drawn behind the markers


def PlotConcentrations(time,data,xlabel,ylabel,labels):
    for i in range(0,25,1):
        fig, axis = plt.subplots()
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        custom_cycler = cycler(color=colors)
        axis.set_prop_cycle(custom_cycler)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.set_xlim([0,5])
        plt.plot(time[:,0],data[:,i],label=labels[i])
        plt.title(labels[i])
    #plt.legend(loc="lower right")
    #plt.savefig("test.pdf")

species = ["O3",
           "O1D",
           "OH",
           "REST",
           "NO2",
           "NO",
           "CH2O",
           "HO2",
           "CO",
           "CO2",
           "CH4",
           "CH3O2",
           "isoprene",
           "RO2",
           "MVK",
           "H2O2",
           "HNO3",
           "NO3",
           "N2O5",
           "SO2",
           "H2SO4",
           "H2SO4_P",
           "alpha-pinene",
           "HNO3_P",
           "ELVOC"
           ]
    
conc = tools.ReadGeneralData("output/concentrations.dat")
time = tools.ReadGeneralData("output/time.dat")
PlotConcentrations(time, conc, "time [days]", "concentration [molecules*cm$^-3$*s$^-1$]", species)