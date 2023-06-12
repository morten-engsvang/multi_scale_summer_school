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


def Plot(time,data1,data2,data3,xlabel,ylabel,labels,title):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    #axis.set_xlim([3,5])
    plt.plot(time[:,0], data1[:,0],label=labels[0])
    plt.plot(time[:,0], data2[:,0],label=labels[1])
    plt.plot(time[:,0], data3[:,0],label=labels[2])
    plt.legend()
    plt.title(title)
    plt.grid(visible=True)


time = tools.ReadGeneralData("output/time.dat")

PN_nuc = tools.ReadGeneralData("old_data/Only_nucleation/PN.dat")
PV_nuc = tools.ReadGeneralData("old_data/Only_nucleation/PV.dat")

PN_nuc_cond = tools.ReadGeneralData("old_data/nucleation_condensation/PN.dat")
PV_nuc_cond = tools.ReadGeneralData("old_data/nucleation_condensation/PV.dat")

PN_nuc_cond_coag = tools.ReadGeneralData("output/PN.dat")
PV_nuc_cond_coag = tools.ReadGeneralData("output/PV.dat")

labels = ["Only Nucleation","Nucleation + Condensation", "Nuc+Cond+Coag"]
Plot(time,PN_nuc,PN_nuc_cond,PN_nuc_cond_coag,"Time [days]","Total PN [cm$^{-3}$]",labels,"PN")
Plot(time,PV_nuc,PV_nuc_cond,PV_nuc_cond_coag,"Time [days]","Total PV [um$^3$$\cdot$cm$^{-3}$]",labels,"PV")