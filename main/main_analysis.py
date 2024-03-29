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


def Plot(hh,data,xlabel,ylabel,day,title):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
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
    plt.title(title)
    plt.legend()
    #plt.savefig("three_body.pdf", format="pdf", bbox_inches="tight")

def PlotEmissions(time,data1,data2,xlabel,ylabel,labels):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_xlim([3,5])
    
    plt.plot(time[:,0],data1[:,0], label=labels[0])
    plt.plot(time[:,0],data2[:,0], label=labels[1])
    plt.legend(loc="lower right")
    #plt.savefig("test.pdf")

def PlotHeatMap(hh,time,data,xlabel,ylabel,title):
    fig, axis = plt.subplots()
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    height = np.zeros(len(hh[0])-1)
    i = 0
    for element in height:
        half_height = (hh[0,i+1] + hh[0,i]) / 2
        #print(half_height)
        height[i] = half_height
        i += 1
    plt.pcolormesh(time[:,0],height,np.transpose(data),vmax=175,vmin=-1,cmap="jet")
    plt.colorbar()
    plt.title(title)
    
def PlotHeatMap2(diameter,time,data,xlabel,ylabel,title):
    fig, axis = plt.subplots()
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_yscale('log')
    axis.set_xlim([3,5])
    plt.pcolormesh(time[:,0],diameter[0,:]*10**(9),np.log10(np.transpose(data)*10**(-6)),vmax=3.2,vmin=0,cmap="jet")
    plt.colorbar()
    plt.title(title)
    
def PlotHeatMap3(hh,time,data,xlabel,ylabel,title,vmax,vmin):
    fig, axis = plt.subplots()
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    #axis.set_yscale('log')
    axis.set_xlim([3,5])
    plt.pcolormesh(time[:,0],hh[0,:],np.transpose(data),cmap="jet",vmax=vmax,vmin=vmin)
    plt.colorbar()
    plt.title(title)
    
def PlotHeatMapNorm(hh,time,data,xlabel,ylabel,title):
    fig, axis = plt.subplots()
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    height = np.zeros(len(hh[0])-1)
    i = 0
    for element in height:
        half_height = (hh[0,i+1] + hh[0,i]) / 2
        #print(half_height)
        height[i] = half_height
        i += 1
    
    plt.pcolormesh(time[:,0],height,np.transpose(data),norm=plc.PowerNorm(vmin=-1,vmax=300,gamma=0.5),cmap="jet")
    plt.colorbar(ticks=[-1,0,5,50,100,200,300])
    plt.title(title)
    
def PlotHeatMapNorm2(hh,time,data,xlabel,ylabel,mini,maxi,title):
    fig, axis = plt.subplots()
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    plt.pcolormesh(time[:,0],hh[0,:],np.transpose(data),norm=plc.PowerNorm(vmin=mini,vmax=maxi,gamma=0.3),cmap="jet")
    plt.colorbar()
    plt.title(title)
    axis.set_xlim([3,5])
    
def PlotConcentrations(time,data,xlabel,ylabel,species):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_xlim([3,5])
    plt.plot(time[:,0], data[:,1], label = "10 m")
    plt.plot(time[:,0], data[:,5], label = "50 m")
    plt.title(species)
    plt.legend()

def PlotConcDistr(diameter,data1,data2,data3,data4,xlabel,ylabel,title):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_xscale('log')
    axis.set_yscale('log')
    plt.plot(diameter[0,:]*10**(9), data1[120,:]*10**(-6), label="Only nuc")
    plt.plot(diameter[0,:]*10**(9), data2[120,:]*10**(-6), label="Nuc+Cond")
    plt.plot(diameter[0,:]*10**(9), data3[120,:]*10**(-6), label="Nuc+Cond+Coag")
    plt.plot(diameter[0,:]*10**(9), data4[120,:]*10**(-6), label="Nuc+Cond+Coag+Dep")
    axis.set_xlim([1,10**3])
    axis.set_ylim([1,10**6])
    plt.grid(visible=True,which="both")
    plt.legend()
    plt.title(title)

def PlotTimeSeries(time,data1,data2,data3,xlabel,ylabel,title):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    plt.plot(time[:,0],data1[:,1], label = "Nuc+Cond")
    plt.plot(time[:,0],data2[:,1], label = "Nuc+Cond+Coag")
    plt.plot(time[:,0],data3[:,1], label = "Nuc+Cond+Coag+Dep")
    axis.set_xlim([3,5])
    plt.title(title)
    plt.legend()
    
def PlotHeightProfile(hh,data1,data2,data3,xlabel,ylabel,title):
    fig, axis = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    custom_cycler = cycler(color=colors)
    axis.set_prop_cycle(custom_cycler)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    plt.plot(data1[-1,:], hh[0,:] , label = "Nuc+Cond")
    plt.plot(data2[-1,:], hh[0,:], label = "Nuc+Cond+Coag")
    plt.plot(data3[-1,:], hh[0,:],  label = "Nuc+Cond+Coag+Dep")
    plt.title(title)
    plt.legend()

hh = tools.ReadGeneralData("output/hh.dat")
uwind = tools.ReadGeneralData("output/uwind.dat")
vwind = tools.ReadGeneralData("output/vwind.dat")
theta = tools.ReadGeneralData("output/theta.dat")
K_m = tools.ReadGeneralData("output/K_m.dat")
K_h = tools.ReadGeneralData("output/K_h.dat")
richard = tools.ReadGeneralData("output/richard.dat")
time = tools.ReadGeneralData("output/time.dat")
diameter = tools.ReadGeneralData("output/diameter.dat")
aero_conc_nuc = tools.ReadGeneralData("old_data/nuc/aerosol_conc_1.dat")
aero_conc_nuc_cond = tools.ReadGeneralData("old_data/nuc_cond/aerosol_conc_1.dat")
aero_conc_nuc_cond_coag = tools.ReadGeneralData("old_data/nuc_cond_coag/aerosol_conc_1.dat")
aero_conc_dep = tools.ReadGeneralData("output/aerosol_conc_1.dat")

time_old = tools.ReadGeneralData("old_data/nuc_cond/time.dat")
diameter_old = tools.ReadGeneralData("old_data/nuc_cond/diameter.dat")

PN_nuc_cond = tools.ReadGeneralData("old_data/nuc_cond/PN.dat")
PN_nuc_cond_coag = tools.ReadGeneralData("old_data/nuc_cond_coag/PN.dat")
PN_dep = tools.ReadGeneralData("output/PN.dat")

PM_nuc_cond = tools.ReadGeneralData("old_data/nuc_cond/PM.dat")
PM_nuc_cond_coag = tools.ReadGeneralData("old_data/nuc_cond_coag/PM.dat")
PM_dep = tools.ReadGeneralData("output/PM.dat")


emission_isoprene = tools.ReadGeneralData("output/emission_isoprene.dat")
emission_monoterpene = tools.ReadGeneralData("output/emission_monoterpene.dat")

alpha_pinene = tools.ReadGeneralData("output/alpha_pinene.dat")
alpha_pinene_old = tools.ReadGeneralData("old_data/nuc_cond_coag/alpha_pinene.dat")
isoprene = tools.ReadGeneralData("output/isoprene.dat")
isoprene_old = tools.ReadGeneralData("old_data/nuc_cond_coag/isoprene.dat")
oh_radical = tools.ReadGeneralData("output/oh_radical.dat")
oh_radical_old = tools.ReadGeneralData("old_data/nuc_cond_coag/oh_radical.dat")
ho2_radical = tools.ReadGeneralData("output/ho2_radical.dat")
ho2_radical_old = tools.ReadGeneralData("old_data/nuc_cond_coag/ho2_radical.dat")
h2so4 = tools.ReadGeneralData("output/h2so4.dat")
h2so4_old = tools.ReadGeneralData("old_data/nuc_cond_coag/h2so4.dat")
elvoc = tools.ReadGeneralData("output/elvoc.dat")
elvoc_old = tools.ReadGeneralData("old_data/nuc_cond_coag/elvoc.dat")

Plot(hh,uwind,"Wind speed [m/s]","Height [m]",1,"U-wind")
Plot(hh,uwind,"Wind speed [m/s]","Height [m]",5,"U-wind")
Plot(hh,vwind,"Wind speed [m/s]","Height [m]",1,"V-wind")
Plot(hh,vwind,"Wind speed [m/s]","Height [m]",5,"V-wind")
Plot(hh,theta,"Potential Temperature [K]","Height [m]",1,"Potential Temperature")
Plot(hh,theta,"Potential Temperature [K]","Height [m]",5,"Potential Temperature")
PlotHeatMap(hh,time,K_m,"Time [days]","Height [m]","K_m")
PlotHeatMap(hh,time,K_h,"Time [days]","Height [m]","K_h")
PlotHeatMap(hh,time,richard,"Time [days]","Height [m]","Ri")
PlotHeatMapNorm(hh,time,richard,"Time [days]","Height [m]","Ri")
labels = ["isoprene","monoterpene"]
PlotEmissions(time, emission_isoprene, emission_monoterpene, "Time [days]", "Emission Rate [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]",labels)
PlotConcentrations(time, alpha_pinene_old, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "alpha-pinene nuc+cond+coag")
PlotConcentrations(time, alpha_pinene, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "alpha-pinene")
PlotConcentrations(time, isoprene_old, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "isoprene nuc+cond+coag")
PlotConcentrations(time, isoprene, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "isoprene")
PlotConcentrations(time, oh_radical_old, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "OH nuc+cond+coag")
PlotConcentrations(time, oh_radical, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "OH")
PlotConcentrations(time, ho2_radical_old, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "HO2 nuc+cond+coag")
PlotConcentrations(time, ho2_radical, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "HO2")
PlotConcentrations(time, h2so4_old, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "H2SO4 nuc+cond+coag")
PlotConcentrations(time, h2so4, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "H2SO4")
PlotConcentrations(time, elvoc_old, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "ELVOC nuc+cond+coag")
PlotConcentrations(time, elvoc, "time [days]", "concentration [molecules$\cdot$cm$^{-3}\cdot$s$^{-1}$]", "ELVOC")
title = "OH (" + str(np.min(oh_radical)) + ", " + str(np.max(oh_radical)) + ")"
PlotHeatMapNorm2(hh, time, oh_radical, "time [days]", "height [m]", 0, 2*10**(6), title)
title = "HO2 (" + str(np.min(ho2_radical)) + ", " + str(np.max(ho2_radical)) + ")"
PlotHeatMapNorm2(hh, time, ho2_radical, "time [days]", "height [m]", 0, 1.2*10**(8), title)
title = "H2SO4 (" + str(np.min(h2so4)) + ", " + str(np.max(h2so4)) + ")"
PlotHeatMapNorm2(hh, time, h2so4, "time [days]", "height [m]", 0, 3.5*10**(7), title)
title = "ELVOC (" + str(np.min(elvoc)) + ", " + str(np.max(elvoc)) + ")"
PlotHeatMapNorm2(hh, time, elvoc, "time [days]", "height [m]", 0, 1.05*10**(7), title)
title = "ALPHA (" + str(np.min(alpha_pinene)) + ", " + str(np.max(alpha_pinene)) + ")"
PlotHeatMapNorm2(hh, time, alpha_pinene, "time [days]", "height [m]", 0, 2.4*10**(9), title)
title = "ISOPRENE (" + str(np.min(isoprene)) + ", " + str(np.max(isoprene)) + ")"
PlotHeatMapNorm2(hh, time, isoprene, "time [days]", "height [m]", 0, 1.6*10**(9), title)


PlotConcDistr(diameter_old,aero_conc_nuc,aero_conc_nuc_cond,aero_conc_nuc_cond_coag,aero_conc_dep,"Diameter [nm]","N [cm$^{-3}$]","Particle Size Distribution for the first layer after simulation")
PlotTimeSeries(time_old, PN_nuc_cond, PN_nuc_cond_coag, PN_dep, "Time [days]", "Total PN [cm$^{-3}$]", "PN in the first model layer")
PlotTimeSeries(time_old, PM_nuc_cond, PM_nuc_cond_coag, PM_dep, "Time [days]", "Total PM [$\mu$g cm$^{-3}$]", "PM in the first model layer")
PlotHeatMap2(diameter_old,time_old,aero_conc_nuc_cond_coag,"Time [days]","Particle diameter [nm]","PN [log$_{10}$(cm$^{-3}$) vs. time, nuc+cond+coag")
PlotHeatMap2(diameter_old,time_old,aero_conc_dep,"Time [days]","Particle diameter [nm]","PN [log$_{10}$(cm$^{-3}$) vs. diameter and time")
PlotHeatMap3(hh,time_old,PN_nuc_cond_coag,"Time [days]","Height [m]","PN vs. height and time",4*10**4,0)
PlotHeatMap3(hh,time_old,PN_dep,"Time [days]","Height [m]","PN vs. height and time",4*10**4,0)
PlotHeatMap3(hh,time_old,PM_nuc_cond_coag,"Time [days]","Height [m]","PM vs. height and time",1.5,1.2)
PlotHeatMap3(hh,time_old,PM_dep,"Time [days]","Height [m]","PM vs. height and time",1.35,1.2)

PlotHeightProfile(hh,PN_nuc_cond,PN_nuc_cond_coag,PN_dep,"PN [(cm$^{-3}$)","Height [m]","")
PlotHeightProfile(hh,PM_nuc_cond,PM_nuc_cond_coag,PM_dep,"PM [$\mu$g cm$^{-3}$]","Heigh [m]","")