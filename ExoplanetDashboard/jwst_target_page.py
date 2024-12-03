#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 19:40:00 2024

@author: jlothringer
"""

import pandas as pd 
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import csv
from datetime import datetime

cs = ['#1f77b4', '#ff7f0e', '#d62728','#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def ra_to_deg(hours,minutes,seconds):
    #return 360*hours/24+minutes/60+seconds/60/60
    return (hours+minutes/60+seconds/60/60)*15

def dec_to_deg(degrees,minutes,seconds):
    return degrees+(sign(degrees)*minutes)/60+(sign(degrees)*seconds)/60/60

data = pd.read_csv('trexolists_jwst.csv')

archived = np.where(data['Status'] == 'Archived')[0]
planned = np.where(data['Status'] != 'Archived')[0]

transits = np.where(data['Event'] == 'Transit')[0]
eclipses = np.where(data['Event'] == 'Eclipse')[0]
phasecurves = np.where(data['Event'] == 'PhaseC')[0]

data_ra = [ra_to_deg(float(data['R.A. 2000'][i][0:2]),float(data['R.A. 2000'][i][3:5]),float(data['R.A. 2000'][i][6:])) for i in range(len(data['R.A. 2000']))]

data_dec = [dec_to_deg(float(data['Dec. 2000'][i][0:3]),float(data['Dec. 2000'][i][4:6]),float(data['Dec. 2000'][i][7:])) for i in range(len(data['R.A. 2000']))]

#planets = pd.read_csv('allinfo-csv.csv')
planets = pd.read_csv('PSCompPars_2024.10.26_05.54.48.csv')
eaot = pd.read_csv('search.csv')
ra_planets_rad = planets['ra']*np.pi/180.
dec_planets_rad = planets['dec']*np.pi/180.

linked_data = []
tsm = []
esm = []
esm5 = []
for i in range(len(data)):
    ra_targ_rad = data_ra[i]*np.pi/180
    dec_targ_rad = data_dec[i]*np.pi/180
    tmp = np.sin(dec_planets_rad)*np.sin(dec_targ_rad) + np.cos(dec_planets_rad)*np.cos(dec_targ_rad)*cos(ra_planets_rad-ra_targ_rad)
    r = np.arccos(tmp)
    #planets.iloc[argmin(r)]
    inds = np.where(r == np.min(r))[0]
    if len(inds) > 1:
        #multi-planet ID
        
        if data.iloc[i]['Event'] != 'PhaseC':
            dur = data.iloc[i]['Hours']
            #overhead ~40%
            Tdur_est = dur*0.65 / 2 + 0.5
            Tdur = planets.iloc[inds]['pl_trandur']
            linked_data.append(inds[argmin(Tdur_est - Tdur)])
        else:
            dur = data.iloc[i]['Hours']
            Per_est = dur*0.65 - 4.5
            Per = planets.iloc[inds]['pl_orbper']*24
            linked_data.append(inds[argmin(Per_est - Per)])
    else:
        linked_data.append(argmin(r))
    
    try:
        eaot_ind = np.where(eaot['Planet_Name'] == planets.iloc[linked_data[i]]['pl_name'])[0][0]
        tsm.append(eaot['SNR_Transmission_K_mag'][eaot_ind])
        esm.append(eaot['SNR_Emission_15_micron'][eaot_ind])
        esm5.append(eaot['SNR_Emission_5_micron'][eaot_ind])
    except:
        print('Having trouble with '+str(i))
        print(data.iloc[i]['Target'])
        tsm.append(0.0)
        esm.append(0.0)
        esm5.append(0.0)
        linked_data[i] = nan
        
#clean TSMs
tsm = np.asarray(tsm)
esm = np.asarray(esm)
esm5 = np.asarray(esm5)
tsm[np.where(np.isnan(tsm))[0]] = -1.0
esm[np.where(np.isnan(esm))[0]] = -1.0
esm5[np.where(np.isnan(esm5))[0]] = -1.0

    
linked_data_dirty = np.array(linked_data)
linked_data = linked_data_dirty[~np.isnan(linked_data_dirty)]
    
SEs = np.where(planets['pl_rade'][linked_data] < 1.7)[0]
Neps = np.where((planets['pl_rade'][linked_data] > 1.7) & (planets['pl_rade'][linked_data] < 5.0))[0]
HJs = np.where(planets['pl_rade'][linked_data] > 5.0)[0]

uniq_planet = len(unique(linked_data))
print('Planets observed:',len(unique(linked_data_dirty[archived])),'<br>')
print('Total hours observed:',round(np.sum(data['Hours'][archived])),'<br>')

uniq_planned = np.array([linked_data_dirty[planned][i] in linked_data_dirty[archived] for i in range(len(planned))])

#print('Total planets planned:',len(np.where(uniq_planned == False)[0]))
print('Additional planets planned: ',uniq_planet - len(unique(linked_data_dirty[archived])),'<br>')
#print('Total hours observed:',np.sum(data['Hours'][planned]))


print('Total rocky planets observed (incl. planned):',len(unique(linked_data[SEs])),'<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(SEs) & set(archived))])),' hours of rocky exoplanet observation <br>')

print('Total Neptune-like planets observed (incl. planned):',len(unique(linked_data[Neps])), '<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(Neps) & set(archived))])),' hours of Neptune-like exoplanet observation <br>')

print('Total gas giants planets observed (incl. planned):',len(unique(linked_data[HJs])),'<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(HJs) & set(archived))])),' hours of gas giant exoplanet observation <br>')

out = []
out.append([len(unique(linked_data_dirty[archived]))])
out.append([round(np.sum(data['Hours'][archived]))])
out.append([uniq_planet - len(unique(linked_data_dirty[archived]))])
out.append([len(unique(linked_data[SEs]))])
out.append([len(unique(linked_data[Neps]))])
out.append([len(unique(linked_data[HJs]))])
out.append([round(np.sum(data['Hours'][list(set(SEs) & set(archived))]))])

out.append([round(np.sum(data['Hours'][list(set(Neps) & set(archived))]))])

out.append([round(np.sum(data['Hours'][list(set(HJs) & set(archived))]))])

with open('Dashboard_Block05.txt', 'w', newline='') as mycsvfile:
    writer = csv.writer(mycsvfile, quoting=csv.QUOTE_NONE)
    writer.writerows(out)


#Temperature histogram
bins = np.linspace(0,4500,20)

tmp = [planets.iloc[unique(linked_data[SEs])]['pl_eqt'],
       planets.iloc[unique(linked_data[Neps])]['pl_eqt'],
       planets.iloc[unique(linked_data[HJs])]['pl_eqt']]

plt.figure()
plt.hist(tmp,bins,stacked=True)

plt.minorticks_on()
plt.legend(labels=['Super-Earths','Neptunes','Gas Giants'],fontsize=12)

plt.xlabel('Equilibrium Temperature (K)',fontsize=15)
plt.ylabel('# of planets observed',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('Teq_hist.png')


#Mass histogram
bins = np.linspace(0,1000,30)
bins = np.logspace(-0.7,3)

tmp = [planets.iloc[unique(linked_data[SEs])]['pl_bmasse'],
       planets.iloc[unique(linked_data[Neps])]['pl_bmasse'],
       planets.iloc[unique(linked_data[HJs])]['pl_bmasse']]

plt.figure()
plt.hist(tmp,bins,stacked=True)
plt.xscale('log')

plt.minorticks_on()
plt.legend(labels=['Super-Earths','Neptunes','Gas Giants'],fontsize=12)

plt.xlabel('Earth Masses',fontsize=15)
plt.ylabel('# of planets observed',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('mass_hist.png')


#Radius histogram
bins = np.linspace(0,23,24)

tmp = [planets.iloc[unique(linked_data[SEs])]['pl_rade'],
       planets.iloc[unique(linked_data[Neps])]['pl_rade'],
       planets.iloc[unique(linked_data[HJs])]['pl_rade']]

plt.figure()
plt.hist(tmp,bins,stacked=True)

plt.minorticks_on()
plt.legend(labels=['Super-Earths','Neptunes','Gas Giants'],fontsize=12)

plt.xlabel('Earth Radii',fontsize=15)
plt.ylabel('# of planets observed',fontsize=15)
plt.tick_params(labelsize=14)
plt.xlim(0,23)
plt.tight_layout()

savefig('radius_hist.png')


#Planet population
figure(figsize=(9,6))
plot(planets['pl_orbper'],planets['pl_bmassj'],'.',color='k',label='All Exoplanets \n(with measured mass)')

plot(planets.iloc[unique(linked_data_dirty[archived])]['pl_orbper'],planets.iloc[unique(linked_data_dirty[archived])]['pl_bmassj'],'o',color='dodgerblue',label='Observed with JWST')

plot(planets.iloc[unique(linked_data_dirty[planned])]['pl_orbper'],planets.iloc[unique(linked_data_dirty[planned])]['pl_bmassj'],'o',color='orange',label='Planned with JWST')


plt.yscale('log')
plt.xscale('log')

plt.xlim(1e-1,1e5)

plt.legend(fontsize=14)
plt.xlabel('Orbital Period (days)',fontsize=15)
plt.ylabel('Planet Mass (Jupiter Masses)',fontsize=15)
plt.tick_params(labelsize=14)
plt.title('Last updated: '+datetime.today().strftime('%Y-%m-%d'), fontsize=12)
plt.tight_layout()

savefig('population.png',dpi=200)


#TSM planets vs. mass

plt.figure()

plt.plot(eaot['Mp'],eaot['SNR_Transmission_K_mag'],'k.',label='All Exoplanets \n(with measured mass)')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmassj'],np.asarray(tsm)[transits],'o',color='dodgerblue',label='Observed with JWST')

plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-2,2e0)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Jupiter Masses)',fontsize=15)
plt.ylabel('Transit Spectroscopy Metric',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('TSM_all.png')

#TSM planets vs mass
#small
plt.figure()

plt.plot(eaot['Mp']*317.8,eaot['SNR_Transmission_K_mag'],'k.',label='All Exoplanets \n(with measured mass)')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmasse'],np.asarray(tsm)[transits],'o',color='dodgerblue',label='Observed with JWST')

#plt.xscale('log')
plt.yscale('log')
plt.ylim(0.02,1)
plt.xlim(0,20)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Earth Masses)',fontsize=15)
plt.ylabel('Transit Spectroscopy Metric',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('TSM_smallplanet.png')

#ESM planets vs. mass

plt.figure()

plt.plot(eaot['Mp'],eaot['SNR_Emission_15_micron'],'k.',label='All Exoplanets \n(with measured mass)')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmassj'],np.asarray(esm)[transits],'o',color='dodgerblue',label='Observed with JWST')

plt.xscale('log')
plt.yscale('log')
#plt.ylim(1e-2,2e0)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Jupiter Masses)',fontsize=15)
plt.ylabel(r'Emission Spectroscopy Metric (1.5$\mu$m)',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('ESM_all.png')

#ESM small planets

plt.figure()

plt.plot(eaot['Mp']*317.8,eaot['SNR_Emission_5_micron'],'k.',label='All Exoplanets \n(with measured mass)')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmasse'],np.asarray(esm)[transits],'o',color='dodgerblue',label='Observed with JWST')

#plt.xscale('log')
plt.yscale('log')
#plt.ylim(1e-2,2e0)
plt.xlim(0,20)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Jupiter Masses)',fontsize=15)
plt.ylabel(r'Emission Spectroscopy Metric (1.5$\mu$m)',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('ESM_smallplanet.png')

#output html
out = []
for i in range(len(data)):
    out.append(['<tr>'])
    out.append(['<td> '+str(data['Target'][i])+' </td>'])
        
    try:
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['pl_name'])+' </td>'])
        
        out.append(['<td> '+str(tsm[i])+' </td>'])
        out.append(['<td> '+str(esm5[i])+' </td>'])

        
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['pl_bmassj'])+' </td>'])
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['pl_bmasse'])+' </td>'])
    
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['pl_radj'])+' </td>'])
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['pl_rade'])+' </td>'])
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['pl_eqt'])+' </td>'])
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['st_spectype'])+' </td>'])
        out.append(['<td> '+str(planets.iloc[int(linked_data_dirty[i])]['sy_vmag'])+' </td>'])
    except:
        if data['Target'][i] == 'GJ-341':
            out.append(['<td> Not in confirmed planets </td>'])
        else:
            out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])
        out.append(['<td> - </td>'])



    out.append(['<td> '+str(data['Event'][i])+' </td>'])
    out.append(['<td> '+str(data['Mode'][i])+' </td>'])   
    out.append(['<td> '+str(data['Status'][i])+' </td>'])

    out.append(['<td> '+str(data['Category'][i])+' </td>'])
    out.append(['<td> '+str(data['Program'][i])+' </td>'])
    out.append(['<td> '+str(data['PI name'][i])+' </td>'])


    out.append(['</tr>'])
    
#b = open('Dashboard_Block2.html','w',newline='')
#a = csv.writer(b,delimiter=' ',quoting=csv.QUOTE_NONE)

with open('Dashboard_Block2.html', 'w', newline='') as mycsvfile:
    writer = csv.writer(mycsvfile, quoting=csv.QUOTE_NONE)
    writer.writerows(out)
#b.close()

#output csv
out = []
out.append(['Target','Planet Name','TSM','ESM (5um)','Mass (M_J)','Mass (M_E)','Radius (R_J)','Radius (R_E)','T_eq','Star SpT','Star Vmag','Event','Mode','Status','Category','Program','PI Name'])
for i in range(len(data)):
    try:
        tmp = [data['Target'][i],planets.iloc[int(linked_data_dirty[i])]['pl_name'],tsm[i],esm5[i],
               planets.iloc[int(linked_data_dirty[i])]['pl_bmassj'],planets.iloc[int(linked_data_dirty[i])]['pl_bmasse'],
               planets.iloc[int(linked_data_dirty[i])]['pl_radj'],planets.iloc[int(linked_data_dirty[i])]['pl_rade'],
               planets.iloc[int(linked_data_dirty[i])]['pl_eqt'],planets.iloc[int(linked_data_dirty[i])]['st_spectype'],
               planets.iloc[int(linked_data_dirty[i])]['sy_vmag'],data['Event'][i],data['Mode'][i],data['Status'][i],
               data['Category'][i],data['Program'][i],data['PI name'][i]]
    except:
        tmp = [data['Target'][i],'-','-','-',
               '-','-',
               '-','-',
               '-','-',
               '-',data['Event'][i],data['Mode'][i],data['Status'][i],
               data['Category'][i],data['Program'][i],data['PI name'][i]]
    out.append(tmp)
        
with open('JWST_ExoDashboard_Table.csv', 'w') as mycsvfile:
    writer = csv.writer(mycsvfile)
    writer.writerows(out)