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


def Plot(time,time_old,data1,data2,data3,xlabel,ylabel,labels,title):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    #axis.set_xlim([3,5])
    plt.plot(time_old[:,0], data1[:,0],label=labels[0])
    plt.plot(time_old[:,0], data2[:,0],label=labels[1])
    plt.plot(time[:,0], data3[:,0],label=labels[2])
    plt.legend()
    plt.title(title)
    plt.grid(visible=True)

def PlotConc(diameter,data1,data2,data3,xlabel,ylabel,labels,title):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.plot(diameter[0,:]*10**(9), data1[24,:]*10**(-6),label=labels[0])
    plt.plot(diameter[0,:]*10**(9), data2[24,:]*10**(-6),label=labels[1])
    plt.plot(diameter[0,:]*10**(9), data3[24,:]*10**(-6),label=labels[2])
    plt.grid(visible=True)
    plt.title(title)

time = tools.ReadGeneralData("output/time.dat")
time_old = tools.ReadGeneralData("old_data/Only_nucleation/time.dat")
diameter = tools.ReadGeneralData("output/diameter.dat")

PN_nuc = tools.ReadGeneralData("old_data/Only_nucleation/PN.dat")
PV_nuc = tools.ReadGeneralData("old_data/Only_nucleation/PV.dat")
conc_nuc = tools.ReadGeneralData("old_data/Only_nucleation/particle_conc.dat")

PN_nuc_cond = tools.ReadGeneralData("old_data/nucleation_condensation/PN.dat")
PV_nuc_cond = tools.ReadGeneralData("old_data/nucleation_condensation/PV.dat")
conc_nuc_cond = tools.ReadGeneralData("old_data/nucleation_condensation/particle_conc.dat")

PN_nuc_cond_coag = tools.ReadGeneralData("output/PN.dat")
PV_nuc_cond_coag = tools.ReadGeneralData("output/PV.dat")
conc_nuc_cond_coag = tools.ReadGeneralData("output/particle_conc.dat")

labels = ["Only Nucleation","Nucleation + Condensation", "Nuc+Cond+Coag"]
Plot(time,time_old,PN_nuc,PN_nuc_cond,PN_nuc_cond_coag,"Time [days]","Total PN [cm$^{-3}$]",labels,"PN")
Plot(time,time_old,PV_nuc,PV_nuc_cond,PV_nuc_cond_coag,"Time [days]","Total PV [um$^3$$\cdot$cm$^{-3}$]",labels,"PV")
PlotConc(diameter,conc_nuc,conc_nuc_cond,conc_nuc_cond_coag,"Diameter [nm]","$\Delta$N [cm$^{-3}$]",labels,"Particle Number Distribution After 24 hrs")