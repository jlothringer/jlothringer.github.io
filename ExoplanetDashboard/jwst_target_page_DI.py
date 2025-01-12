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

#cs = ['#1f77b4', '#ff7f0e', '#d62728','#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

plt.style.use('seaborn-v0_8-colorblind')
rcParams['text.usetex'] = True

def ra_to_deg(hours,minutes,seconds):
    #return 360*hours/24+minutes/60+seconds/60/60
    return (hours+minutes/60+seconds/60/60)*15

def dec_to_deg(degrees,minutes,seconds):
    return degrees+(sign(degrees)*minutes)/60+(sign(degrees)*seconds)/60/60

def to_julian_date(timestamp_str):
    # Parse the timestamp string
    dt = datetime.strptime(timestamp_str, "%b_%d_%Y_%H:%M:%S")
    
    # Julian Date calculation
    # Based on algorithm from Jean Meeus' "Astronomical Algorithms"
    year = dt.year
    month = dt.month
    day = dt.day + (dt.hour - 12) / 24.0 + dt.minute / 1440.0 + dt.second / 86400.0
    
    if month <= 2:
        year = year - 1
        month = month + 12
    
    a = int(year / 100)
    b = 2 - a + int(a / 4)
    
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b - 1524.0 #not .5 cuz?
    
    return jd

data_tr = pd.read_csv('trexolists_jwst.csv')
data_tr['Groups'] = data_tr['Groups'].astype(str)


data_di = pd.read_csv('diexolists.csv')
data_di['Groups'] = data_di['Groups'].astype(str)

data = pd.concat([data_tr,data_di],axis=0,ignore_index=True)

#data = data_tr

archived = np.where(data['Status'] == 'Archived')[0]
#planned = np.where(data['Status'] != 'Archived')[0]
planned = np.where((data['Status'] == 'FlightReady') | (data['Status'] == 'Implementation'))[0]

transits = np.where(data['Event'] == 'Transit')[0]
eclipses = np.where(data['Event'] == 'Eclipse')[0]
phasecurves = np.where(data['Event'] == 'PhaseC')[0]
di = np.where(data['Event'] == 'DI')[0]

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
contemp_transit = []
contemp_eclipse = []
for i in range(len(data)):
    contemp_transit.append(None)
    contemp_eclipse.append(None)
    ra_targ_rad = data_ra[i]*np.pi/180
    dec_targ_rad = data_dec[i]*np.pi/180
    tmp = np.sin(dec_planets_rad)*np.sin(dec_targ_rad) + np.cos(dec_planets_rad)*np.cos(dec_targ_rad)*cos(ra_planets_rad-ra_targ_rad)
    #calculate distance to every exoplanet (in either degrees or radians)
    r = np.arccos(tmp)
    #If distance is more than 1/1000th of a degree, then we probably didn't find anything,
    #which would be the case if this observation was of a disk or was a *search* for planets
    #And we're not counting that here... would be cool to add info on that tho.
    if min(r) > 1e-5:
        print('Having trouble with '+str(i))
        print(data.iloc[i]['Target'])
        tsm.append(0.0)
        esm.append(0.0)
        esm5.append(0.0)
        linked_data.append(nan)
    else:
        #planets.iloc[argmin(r)]
        #find minimum distance, but don't use argmin because there may be multiple, indicating a multi-planet system
        inds = np.where(r == np.min(r))[0]
        
        if len(inds) > 1:
            #multi-planet ID
            if data.iloc[i]['Event'] == 'DI':
                #just picking the first planet in the system for DI planets...
                linked_data.append(inds[0])
            elif data.iloc[i]['Event'] != 'PhaseC':
                try:
                    #Not everything has start or end phases
                    start_phase = ((to_julian_date(data.iloc[i]['Start date']) - planets.iloc[inds]['pl_tranmid']) / planets.iloc[inds]['pl_orbper']) % 1              
                    end_phase = ((to_julian_date(data.iloc[i]['End date']) - planets.iloc[inds]['pl_tranmid']) / planets.iloc[inds]['pl_orbper']) % 1
                    
                    eclip = np.where((start_phase < 0.5) & (end_phase > 0.5))[0]
                    trans = np.where((start_phase > 0.5) & (end_phase < 0.5))[0]
                    
                    #if trans
                    contemp_transit[i] = [planets['pl_name'].iloc[inds[trans[j]]] for j in range(len(trans))]
                    contemp_eclipse[i] = [planets['pl_name'].iloc[inds[eclip[j]]] for j in range(len(eclip))]
                    
                    #Chose the eclipse or transit of the outermost if there are contemporaneous events
                    #That's the best we can do without looking into the proposal
                    #We choose outermost b/c that's what Natalie is doing and those events are more valuable/rare
                    if data.iloc[i]['Event'] == 'Transit':
                        linked_data.append(inds[trans[0]])
                    elif data.iloc[i]['Event'] == 'Eclipse':
                        linked_data.append(inds[eclip[0]])
                    else:
                        #There's some stuff that isn't labeled by event b/c they have no phase constraint
                        dur = data.iloc[i]['Hours']
                        #overhead ~40%
                        Tdur_est = dur*0.65 / 2 + 0.5
                        Tdur = planets.iloc[inds]['pl_trandur']
                        linked_data.append(inds[argmin(Tdur_est - Tdur)])
                        
                except:
                    #If we can't find the transit or eclipse, estimate via the duration
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
        print(linked_data[i])
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
 
#This doesn't work- data not correspond with right planet   
#SEs = np.where((planets['pl_rade'][linked_data] < 1.7) & (data['Event'] != 'DI'))[0]
#Neps = np.where((planets['pl_rade'][linked_data] > 1.7) & (planets['pl_rade'][linked_data] < 5.0))[0]
#HJs = np.where((planets['pl_rade'][linked_data] > 5.0) & (data['Event'] != 'DI'))[0]
#DIs = np.where(data['Event'] == 'DI')[0]

SEs = np.where((planets['pl_rade'][linked_data_dirty[np.where((data['Event'] != 'DI') & ~np.isnan(linked_data_dirty))[0]]] < 1.7))[0]
Neps = np.where((planets['pl_rade'][linked_data_dirty[np.where((data['Event'] != 'DI') & ~np.isnan(linked_data_dirty))[0]]] > 1.7) & 
                (planets['pl_rade'][linked_data_dirty[np.where((data['Event'] != 'DI') & ~np.isnan(linked_data_dirty))[0]]] < 5.0))[0]
HJs = np.where((planets['pl_rade'][linked_data_dirty[np.where((data['Event'] != 'DI') & ~np.isnan(linked_data_dirty))[0]]] > 5.0))[0]

transiting_planets = planets['pl_name'][linked_data_dirty[np.where((data['Event'] != 'DI') & ~np.isnan(linked_data_dirty))[0]]]
transit_planets_inds = np.where((data['Event'] != 'DI') & ~np.isnan(linked_data_dirty))[0]
DI_planets = planets['pl_name'][linked_data_dirty[np.where((data['Event'] == 'DI') & ~np.isnan(linked_data_dirty))[0]]]
DI_planets_inds = np.where((data['Event'] == 'DI') & ~np.isnan(linked_data_dirty))[0]


uniq_planet = len(unique(linked_data))
print('Planets observed:',len(unique(linked_data_dirty[archived])),'<br>')
print('Total hours observed:',round(np.sum(data['Hours'][archived])),'<br>')

uniq_planned = np.array([linked_data_dirty[planned][i] in linked_data_dirty[archived] for i in range(len(planned))])

#print('Total planets planned:',len(np.where(uniq_planned == False)[0]))
print('Additional planets planned: ',uniq_planet - len(unique(linked_data_dirty[archived])),'<br>')
#print('Total hours observed:',np.sum(data['Hours'][planned]))


print('Total transiting rocky planets observed (incl. planned):',len(unique(linked_data[SEs])),'<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(SEs) & set(archived))])),' hours of rocky exoplanet observation <br>')

print('Total transiting Neptune-like planets observed (incl. planned):',len(unique(linked_data[Neps])), '<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(Neps) & set(archived))])),' hours of Neptune-like exoplanet observation <br>')

print('Total transiting giants planets observed (incl. planned):',len(unique(linked_data[HJs])),'<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(HJs) & set(archived))])),' hours of gas giant exoplanet observation <br>')

print('Total transiting planets observed (incl. planned):',len(unique(transiting_planets)),'<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(transit_planets_inds) & set(archived))])),' hours of transiting exoplanet observation <br>')

print('Total directly imaged planets observed (incl. planned):',len(unique(DI_planets)),'<br>')
print('for a total of ', round(np.sum(data['Hours'][list(set(DI_planets_inds) & set(archived))])),' hours of characterizing directly imaged exoplanets <br>')
print('for a total of ', round(np.sum(data['Hours'][np.where(data['Event'] == 'DI')[0]])),' hours looking at directly imaged exoplanet systems <br>')


#print('Total transiting planets observed (incl. planned):',len(unique(linked_data[HJs])),'<br>')
#print('for a total of ', round(np.sum(data['Hours'][list(set(HJs) & set(archived))])),' hours of gas giant exoplanet observation <br>')

out = []
out.append([round(np.sum(data['Hours'][archived]))])     # hours spent observing exoplanets
out.append([len(unique(linked_data_dirty[archived]))])   # number of planets observed
out.append([uniq_planet - len(unique(linked_data_dirty[archived]))]) #number of planets still to be observed

out.append([round(np.sum(data['Hours'][np.where((data['Event'] != 'DI') & (data['Status'] == 'Archived'))[0]]))]) #hours spent on transiting systems
out.append([round(np.sum(data['Hours'][np.where((data['Event'] == 'DI') & (data['Status'] == 'Archived'))[0]]))]) #hours spent on DI systems
out.append([round(np.sum(data['Hours'].iloc[np.where(data['Target'] == 'TRAPPIST-1')]))]) #hours spent on Trappist-1


out.append([len(unique(linked_data[SEs]))])  # transiting obs of terrestrial planets
out.append([len(unique(linked_data[Neps]))]) # transiting obs of Neptune-like planets
out.append([len(unique(linked_data[HJs]))])  # Transiting obs of giant planets
out.append([round(np.sum(data['Hours'][list(set(SEs) & set(archived))]))]) #hours on transiting terrestrials
out.append([round(np.sum(data['Hours'][list(set(Neps) & set(archived))]))]) #hours on Neps 
out.append([round(np.sum(data['Hours'][list(set(HJs) & set(archived))]))]) #hours on giant planets 

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
plt.ylabel('\# of planets observed',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('Teq_hist.png',dpi=200)


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
plt.ylabel('\# of planets observed',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('mass_hist.png',dpi=200)


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
plt.ylabel('\# of planets observed',fontsize=15)
plt.tick_params(labelsize=14)
plt.xlim(0,23)
plt.tight_layout()

savefig('radius_hist.png',dpi=200)

#Planet population
figure(figsize=(9,6))

plot(planets['pl_orbper'],planets['pl_bmassj'],'.',color='grey',label='All Exoplanets \n(with measured mass)',alpha=0.35)

plot(planets.iloc[unique(linked_data_dirty[archived])]['pl_orbper'],planets.iloc[unique(linked_data_dirty[archived])]['pl_bmassj'],'H',label='Transiting Planets Observed with JWST',markersize=10,markeredgecolor= "black",zorder=999)

plot(planets.iloc[unique(linked_data_dirty[planned])]['pl_orbper'],planets.iloc[unique(linked_data_dirty[planned])]['pl_bmassj'],'H',label='Planned with JWST',markersize=10,markeredgecolor= "black",zorder=998)

plot(planets.iloc[unique(linked_data_dirty[list(set.intersection(set(archived),set(di)))])]['pl_orbper'],planets.iloc[unique(linked_data_dirty[list(set.intersection(set(archived),set(di)))])]['pl_bmassj'],'H',
     label='Imaged Planets Observed with JWST',markersize=10,markeredgecolor= "black",zorder=999,color='#D55E00')


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

plt.plot(eaot['Mp'],eaot['SNR_Transmission_K_mag'],'.',color='grey',label='All Exoplanets \n(with measured mass)')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmassj'],np.asarray(tsm)[transits],'H',label='Observed with JWST')

plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-2,2e0)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Jupiter Masses)',fontsize=15)
plt.ylabel('Transit Spectroscopy Metric',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('TSM_all.png',dpi=200)

#TSM planets vs mass
#small
plt.figure()

plt.plot(eaot['Mp']*317.8,eaot['SNR_Transmission_K_mag'],'.',label='All Exoplanets \n(with measured mass)',color='grey')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmasse'],np.asarray(tsm)[transits],'H',label='Observed with JWST')

#plt.xscale('log')
plt.yscale('log')
plt.ylim(0.02,1)
plt.xlim(0,20)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Earth Masses)',fontsize=15)
plt.ylabel('Transit Spectroscopy Metric',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('TSM_smallplanet.png',dpi=200)

#ESM planets vs. mass

plt.figure()

plt.plot(eaot['Mp'],eaot['SNR_Emission_15_micron'],'.',label='All Exoplanets \n(with measured mass)',color='grey')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmassj'],np.asarray(esm)[transits],'H',label='Observed with JWST',markersize=10,markeredgecolor= "black")

plt.xscale('log')
plt.yscale('log')
#plt.ylim(1e-2,2e0)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Jupiter Masses)',fontsize=15)
plt.ylabel(r'Emission Spectroscopy Metric (1.5$\mu$m)',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('ESM_all.png',dpi=200)

#ESM small planets

plt.figure()

plt.plot(eaot['Mp']*317.8,eaot['SNR_Emission_5_micron'],'.',label='All Exoplanets \n(with measured mass)',color='grey')
plt.plot(planets.iloc[linked_data_dirty[transits]]['pl_bmasse'],np.asarray(esm)[transits],'H',label='Observed with JWST',markersize=10,markeredgecolor= "black")

#plt.xscale('log')
plt.yscale('log')
#plt.ylim(1e-2,2e0)
plt.xlim(0,20)
plt.legend(fontsize=14)
plt.xlabel('Planet Mass (Jupiter Masses)',fontsize=15)
plt.ylabel(r'Emission Spectroscopy Metric (1.5$\mu$m)',fontsize=15)
plt.tick_params(labelsize=14)
plt.tight_layout()

savefig('ESM_smallplanet.png',dpi=200)

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
    

#output full mega table
#all planet properties for each observation
out = []
out.append(['Target','Planet Name','TSM','ESM (5um)','Mass (M_J)','Mass (M_E)','Radius (R_J)','Radius (R_E)','T_eq','Star SpT','Star Vmag','Event','Mode','Status','Category','Program','PI Name'])
out = list(data.keys())
out.extend(['Contemporaneous Transits','Contemporaneous Eclipses'])
out.extend(planets.keys())
out = [out]

for i in range(len(data)):
    try:
        tmp = [data[key][i] for key in data.keys()]
        
        tmp.extend(contemp_transit)
        tmp.extend(contemp_eclipse)
        
        tmp2 = [planets.iloc[int(linked_data_dirty[i])][key] for key in planets.keys()]
        
        tmp.extend(tmp2)
        
    except:
        tmp = [data[key][i] for key in data.keys()]
        tmp.extend('-')
        tmp.extend('-')
        
        tmp2 = ['-' for key in planets.keys()]
        
        tmp.extend(tmp2)

    out.append(tmp)
        
with open('JWST_ExoDashboard_Table_Full.csv', 'w') as mycsvfile:
    writer = csv.writer(mycsvfile)
    writer.writerows(out)
    
#Pop gif prep
from matplotlib.animation import FuncAnimation
import copy
from collections import deque

# Assuming data, planets, and linked_data_dirty are defined elsewhere

# Data preparation (kept from original)
tmp = copy.copy(data['Start date'])
tmp[np.where(tmp == 'X')[0]] = 'Dec_31_2025_23:00:00'
start_dates = pd.to_datetime(tmp, format='%b_%d_%Y_%H:%M:%S')

new = pd.DataFrame(start_dates)
new = new.join(pd.Series(linked_data_dirty, name='ldd'))
new = new.join(pd.Series(data['Target']))
new = new.join(pd.Series(data['Status']))
new = new.join(pd.Series(data['Event']))
new = new.sort_values(by='Start date')

#manually add some DI planet periods
# 1312    HIP 65426 b
# 1345      HR 2562 b
# 1240     HD 95086 b
# 514     HD 106906 b
# 5757      kap And b
# 4826      MWC 758 c - no mass so not plotted
# 4787      LkCa 15 b - no mass so not plotted
# 49         AB Aur b
# 1256     HD 97048 b
invalid_periods = [1312,1345,1240,514,5757,4826,4787,49,1256]
planets['pl_orbper'].iloc[invalid_periods] = 1e5

archived_new = np.where((new['Status'] == 'Archived') & (~np.isnan(new['ldd'])))[0]
planned_new = np.where((new['Status'] == 'FlightReady') | (new['Status'] == 'Implementation'))[0]

# Create figure and initial plot
fig, ax = plt.subplots(figsize=(9, 6))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(1e-1, 1e6)
ax.set_xlabel('Orbital Period (days)', fontsize=16)
ax.set_ylabel('Planet Mass (Jupiter Masses)', fontsize=16)
ax.tick_params(labelsize=14)
plt.suptitle('JWST Exoplanet Observations', fontsize=16)
title = ax.set_title('2021-12-25 12:20:00', fontsize=16)
ax.grid(True, which="both", linestyle='--', linewidth=0.5)
fig.text(0.82, 0.02, 'Credit: Josh Lothringer')
fig.text(0.891, 0.34, '1 Earth Mass')
fig.text(0.891, 0.45, '1 Neptune Mass')
fig.text(0.891, 0.58, '1 Jupiter Mass')

plt.subplots_adjust(right=0.886)

# Background scatter plot
background = ax.plot(
    planets['pl_orbper'],
    planets['pl_bmassj'],
    '.',
    color='grey',
    label='All Exoplanets',
    alpha=0.35,
    zorder=-2
)[0]

# Empty plots for legend
ax.plot([], [], 'H', markersize=10, markeredgecolor='black', label='Transiting Planets Observed with JWST', color='#0072B2')
ax.plot([], [], 'H', markersize=10, markeredgecolor='black', label='Imaged Planets Observed with JWST', color='#D55E00')
ax.plot([], [], 'H', markersize=10, markeredgecolor='black', label='Planned with JWST', color='#009E73')
plt.legend(fontsize=14, loc='lower right')

# Initialize containers for animated elements
current_point = ax.plot([], [], 'H', markersize=10, markeredgecolor='black', color='gold', zorder=3)[0]
observed_points = []
planet_indices_plotted = set()

# Create a deque to store active labels with their lifetimes
class LabelWithLifetime:
    def __init__(self, text_obj, lifetime):
        self.text_obj = text_obj
        self.lifetime = lifetime

active_labels = deque(maxlen=20)
LABEL_LIFETIME = 5

def get_observation_color(event_type):
    return '#D55E00' if event_type == 'DI' else '#0072B2'

def is_valid_point(planet_idx):
    """Check if the point has valid coordinates"""
    if not np.isfinite(planets.iloc[planet_idx]['pl_orbper']):
        print(f"Invalid orbital period for planet index {planet_idx}")
        return False
    if not np.isfinite(planets.iloc[planet_idx]['pl_bmassj']):
        print(f"Invalid mass for planet index {planet_idx}")
        return False
    return True

def init():
    current_point.set_data([], [])
    return [current_point] + [label.text_obj for label in active_labels]

def update(frame):
    artists = [current_point]
    
    # Update existing labels and remove expired ones
    expired_labels = []
    for label in active_labels:
        label.lifetime -= 1
        if label.lifetime <= 0:
            expired_labels.append(label)
            label.text_obj.set_text('')
        artists.append(label.text_obj)
    
    for label in expired_labels:
        active_labels.remove(label)
    
    if frame < len(archived_new):
        idx = archived_new[frame]
        if ~np.isnan(new.ldd.iloc[idx]):
            try:
                planet_idx = int(new.ldd.iloc[idx])
                
                if is_valid_point(planet_idx):
                    # Update current point
                    current_point.set_data(
                        [planets.iloc[planet_idx]['pl_orbper']],
                        [planets.iloc[planet_idx]['pl_bmassj']]
                    )
                    
                    # Create new label with lifetime
                    new_label = ax.text(
                        planets.iloc[planet_idx]['pl_orbper'] * 1.11,
                        planets.iloc[planet_idx]['pl_bmassj'],
                        planets.iloc[planet_idx]['pl_name'],
                        fontsize=12,
                        zorder=3
                    )
                    active_labels.append(LabelWithLifetime(new_label, LABEL_LIFETIME))
                    artists.append(new_label)
                    
                    # Update title
                    title.set_text(str(new['Start date'].iloc[idx]))
                    
                    # Add previous point to observed set (if it exists)
                    if frame > 0:
                        prev_idx = archived_new[frame-1]
                        if ~np.isnan(new.ldd.iloc[prev_idx]):
                            prev_planet_idx = int(new.ldd.iloc[prev_idx])
                            if is_valid_point(prev_planet_idx):
                                # Get color based on event type
                                point_color = get_observation_color(new['Event'].iloc[prev_idx])
                                point = ax.plot(
                                    planets.iloc[prev_planet_idx]['pl_orbper'],
                                    planets.iloc[prev_planet_idx]['pl_bmassj'],
                                    'H',
                                    markersize=10,
                                    markeredgecolor='black',
                                    color=point_color,
                                    zorder=2
                                )[0]
                                observed_points.append(point)
                                planet_indices_plotted.add(prev_planet_idx)
                else:
                    current_point.set_data([], [])
                
                artists.extend(observed_points)
                
            except (ValueError, IndexError) as e:
                print(f"Error processing frame {frame}: {e}")
                current_point.set_data([], [])
    
    elif frame == len(archived_new) + 5:
        # Add the last observed point to blue/orange set
        last_idx = archived_new[-1]
        if ~np.isnan(new.ldd.iloc[last_idx]):
            last_planet_idx = int(new.ldd.iloc[last_idx])
            if is_valid_point(last_planet_idx) and last_planet_idx not in planet_indices_plotted:
                point_color = get_observation_color(new['Event'].iloc[last_idx])
                point = ax.plot(
                    planets.iloc[last_planet_idx]['pl_orbper'],
                    planets.iloc[last_planet_idx]['pl_bmassj'],
                    'H',
                    markersize=10,
                    markeredgecolor='black',
                    color=point_color,
                    zorder=2
                )[0]
                observed_points.append(point)
        
        # Add planned points
        for idx in planned_new:
            if ~np.isnan(new.ldd.iloc[idx]):
                try:
                    planet_idx = int(new.ldd.iloc[idx])
                    if is_valid_point(planet_idx):
                        point = ax.plot(
                            planets.iloc[planet_idx]['pl_orbper'],
                            planets.iloc[planet_idx]['pl_bmassj'],
                            'H',
                            markersize=10,
                            markeredgecolor='black',
                            color='#009E73',
                            zorder=1
                        )[0]
                        observed_points.append(point)
                except (ValueError, IndexError) as e:
                    print(f"Error processing planned point at index {idx}: {e}")
        
        artists.extend(observed_points)
        current_point.set_data([], [])
    
    return artists

print('MAKING GIF')
with plt.ion():
    ani = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(archived_new) + 30,
        interval=125,
        blit=True
    )
    ani.save('jwst_planets.gif', dpi=150, writer='pillow')
    plt.close('all')
plt.close('all')