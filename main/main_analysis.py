#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import tools
import numpy as np
from cycler import cycler
plt.style.use("default")
plt.rcParams.update({'font.size': 12})
plt.rc("axes", axisbelow = True) #Makes sure gridlines are drawn behind the markers


def Plot(hh,data):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel("Wind speed [m/s]")
    axis.set_ylabel("Height [m]")
    #axis.set_ylabel("Binding free energy " + unit)
    #plt.grid(visible=True)
    #axis.xaxis.set_minor_locator(MultipleLocator(1))
    #axis.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30])
    plt.plot(data[0,:],np.transpose(hh), label = "Time = 0 h")
    plt.plot(data[0,:],np.transpose(hh), label = "Time = 4 h")
    plt.plot(data[10,:],np.transpose(hh), label = "Time = 10 h")
    plt.plot(data[30,:],np.transpose(hh), label = "Time = 30 h")
    plt.plot(data[50,:],np.transpose(hh), label = "Time = 50 h")
    plt.plot(data[70,:],np.transpose(hh), label = "Time = 70 h")
    plt.plot(data[90,:],np.transpose(hh), label = "Time = 90 h")
    
    plt.legend()
    #plt.savefig("three_body.pdf", format="pdf", bbox_inches="tight")
    
    
hh = tools.ReadGeneralData("output/hh.dat")
uwind = tools.ReadGeneralData("output/uwind.dat")
vwind = tools.ReadGeneralData("output/vwind.dat")
Plot(hh,uwind)
Plot(hh,vwind)