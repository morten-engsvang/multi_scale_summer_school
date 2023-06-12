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


def Plot(time,data,xlabel,ylabel,label):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    #axis.set_xlim([3,5])
    plt.plot(time[:,0], data[:,0],label=label)
    plt.legend()
    plt.grid(visible=True)


time = tools.ReadGeneralData("output/time.dat")
PN = tools.ReadGeneralData("output/PN.dat")
PV = tools.ReadGeneralData("output/PV.dat")

Plot(time,PN,"Time [days]","Total PN [cm$^{-3}$]","Only Nucleation")
Plot(time,PV,"Time [days]","Total PV [um$^3$$\cdot$cm$^{-3}$]","Only Nucleation")