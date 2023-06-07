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


def Plot(hh,data,xlabel,ylabel):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    #axis.set_ylabel("Binding free energy " + unit)
    #plt.grid(visible=True)
    #axis.xaxis.set_minor_locator(MultipleLocator(1))
    #axis.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30])
    plt.plot(data[0,:],np.transpose(hh), label = "Time = 0 h")
    plt.plot(data[4,:],np.transpose(hh), label = "Time = 4 h")
    plt.plot(data[8,:],np.transpose(hh), label = "Time = 8 h")
    plt.plot(data[10,:],np.transpose(hh), label = "Time = 10 h")
    plt.plot(data[40,:],np.transpose(hh), label = "Time = 40 h")
    plt.plot(data[80,:],np.transpose(hh), label = "Time = 80 h")
    plt.plot(data[120,:],np.transpose(hh), label = "Time = 120 h")
    
    plt.legend()
    #plt.savefig("three_body.pdf", format="pdf", bbox_inches="tight")
    
    
hh = tools.ReadGeneralData("output/hh.dat")
uwind = tools.ReadGeneralData("output/uwind.dat")
vwind = tools.ReadGeneralData("output/vwind.dat")
theta = tools.ReadGeneralData("output/theta.dat")
Plot(hh,uwind,"Wind speed [m/s]","Height [m]")
Plot(hh,vwind,"Wind speed [m/s]","Height [m]")
Plot(hh,theta,"Potential Temperature [K]","Height [m]")