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


def Plot(hh,data,xlabel,ylabel,day):
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
    if day == 1:
        plt.plot(data[0,:],np.transpose(hh), label = "Time: day 1, hr 0")
        plt.plot(data[4,:],np.transpose(hh), label = "Time: day 1, hr 4")
        plt.plot(data[8,:],np.transpose(hh), label = "Time: day 1, hr 8")
        plt.plot(data[12,:],np.transpose(hh), label = "Time: day 1, hr 12")
        plt.plot(data[16,:],np.transpose(hh), label = "Time: day 1, hr 16")
        plt.plot(data[20,:],np.transpose(hh), label = "Time: day 1, hr 20")
    elif day == 5:
        plt.plot(data[96,:],np.transpose(hh), label = "Time: day 5, hr 0")
        plt.plot(data[100,:],np.transpose(hh), label = "Time: day 5, hr 4")
        plt.plot(data[104,:],np.transpose(hh), label = "Time: day 5, hr 8")
        plt.plot(data[108,:],np.transpose(hh), label = "Time: day 5, hr 12")
        plt.plot(data[112,:],np.transpose(hh), label = "Time: day 5, hr 16")
        plt.plot(data[116,:],np.transpose(hh), label = "Time: day 5, hr 20")
    
    plt.legend()
    #plt.savefig("three_body.pdf", format="pdf", bbox_inches="tight")
    
    
hh = tools.ReadGeneralData("output/hh.dat")
uwind = tools.ReadGeneralData("output/uwind.dat")
vwind = tools.ReadGeneralData("output/vwind.dat")
theta = tools.ReadGeneralData("output/theta.dat")
Plot(hh,uwind,"Wind speed [m/s]","Height [m]",1)
Plot(hh,uwind,"Wind speed [m/s]","Height [m]",5)
Plot(hh,vwind,"Wind speed [m/s]","Height [m]",1)
Plot(hh,vwind,"Wind speed [m/s]","Height [m]",5)
Plot(hh,theta,"Potential Temperature [K]","Height [m]",1)
Plot(hh,theta,"Potential Temperature [K]","Height [m]",5)