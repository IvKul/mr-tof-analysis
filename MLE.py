# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:30:22 2020

@author: ikulikov, jkarthein
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import peakutils
import scipy.optimize as opt
import os
import fnmatch
#here the main parameters to play with
z_cut = 5 #maximum counts of ions in one tof
tof_window = [] #to chose the peak to fit
peak_threshold = 0.9 #the peak threshold
peak_min_distance = 70 #the distance between peaks
nr=1 # number of pints for rolling avavereging
s=2 #how much sigma of the peak to take (5)
fit_ranges = [ 30,50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 300] #ranges to fit
A=600 #amplitude of a gaussian peak
sig=100 #sigma of gaussian peak

#sum up all files in the folder. The header depends on number of "useless" raws in the file (all rows before the [DATA])
path = os.getcwd()
file_list = os.listdir(path)
df_list = []
for file in file_list:
    if fnmatch.fnmatch(file, '*.mpa'): 
        filedata = open(file, 'r')
        data_whole = pd.read_csv(filedata, header=170, float_precision='high')
        data_whole[['tof',"sweep", "counts"]] = pd.DataFrame([ x.split() for x in data_whole["[DATA]"].tolist() ])
        data_whole = data_whole.drop('[DATA]', axis=1).astype(float)
        df_list.append(data_whole)
data_whole = pd.concat(df_list)
'''
#or just one file
data_whole = pd.read_csv('Ar-run190.asc', header=11, float_precision='high')
data_whole[['tof',"sweep", "counts"]] = pd.DataFrame([ x.split() for x in data_whole["[DATA]"].tolist() ])
data_whole = data_whole.drop('[DATA]', axis=1).astype(float)
'''
#computing the tof value. caloff and calfact should be taken from the file
caloff=22190995.200000
calfact=0.8
data_whole['tof']=(caloff+data_whole['tof']*calfact)

#z-cut and windows cut
data_whole = data_whole[data_whole.sweep.isin(data_whole.sweep.value_counts()[data_whole.sweep.value_counts() >= z_cut].index)]
if tof_window != []:
    data_whole = data_whole[(data_whole.tof >= tof_window[0]) & (data_whole.tof <= tof_window[1])]
 #ploting the data 
fig, ((ax_x, blank),(ax_0, ax_y)) = plt.subplots(2,2,sharex='col',sharey='row', figsize=(9,9),
                                                 gridspec_kw={'height_ratios':[1,4],
                                                             'width_ratios':[4,1],
                                                             'hspace': 0.05,
                                                             'wspace':0.05})
blank.remove()
data_whole.plot(x='tof', y='sweep', style='o', alpha=0.15, ms=2,ax=ax_0, label='data')
ax_0.set_xlabel('Time of flight, ns', fontsize=24)
ax_0.set_ylabel('Sweep', fontsize=24)

data_tof=data_whole.sort_values(by=['tof']).groupby(['tof']).aggregate(np.sum)
data_tof=data_tof.reset_index()
ax_x.bar(data_tof.tof, data_tof.counts)
ax_x.set_yscale('log')

data_sweeps=data_whole.sort_values(by=['sweep']).groupby(['sweep']).aggregate(np.sum)
data_sweeps=data_sweeps.reset_index()

ax_y.plot(data_sweeps.counts, data_sweeps.sweep)
ax_0.set_ylabel('sweep', fontsize=24)
ax_x.set_ylabel('counts', fontsize=24)

ax_y.xaxis.set_ticks_position('top')
ax_y.xaxis.set_label_position('top')
ax_0.set_xlabel('time of flight (ns)', fontsize=24)
ax_y.set_ylabel('Y projection', fontsize=24)

#finding the peak
x_proj_peak_ind = peakutils.indexes(data_tof['counts'], thres=peak_threshold, min_dist=peak_min_distance)
peak = pd.DataFrame(data_tof, columns=data_tof.columns, index=x_proj_peak_ind)

#ploting the peak

a_list = []
b_list = []
c_list = []
d_list = []
for peak_pos in peak['tof']:
    ax_x.axvline(peak_pos, c='black')
    ax_x.axvline(peak_pos-300, c='black', alpha=0.5)
    ax_x.axvline(peak_pos+300, c='black', alpha=0.5)
    ax_0.axvline(peak_pos, c='black')
    ax_0.axvline(peak_pos-300, c='black', alpha=0.5)
    ax_0.axvline(peak_pos+300, c='black', alpha=0.5)
 #taking +-300 ns of the data (300ns is the maximum windows we usually take)

for peak_pos in peak['tof']:
    for (x,y,z) in zip(data_whole['tof'],data_whole['sweep'],data_whole['counts']):
        if peak_pos-300<=x<=peak_pos+300:
            a=x
            b=y
            c=z
            d=peak_pos
            a_list.append(a)
            b_list.append(b)
            c_list.append(c)
            d_list.append(d)
data_of_the_peak = pd.DataFrame({'tof': a_list, 'sweeps': b_list, 'counts': c_list, 'peak': d_list})
del a_list, b_list, c_list, d_list

#statr the rolling average of the chosen peak and find the standart deviation
std_value=data_of_the_peak['tof'].rolling(nr).mean().std()

#take the data acording to found sigma (+-3sigma coresponds to 99% of the data)
a_list = []
b_list = []
c_list = []
for peak_pos in peak['tof']:
    for (x,y,z) in zip(data_whole['tof'],data_whole['sweep'],data_whole['counts']):
        if peak_pos-s*std_value<=x<=peak_pos+s*std_value:
            a=x
            b=y
            c=z
            a_list.append(a)
            b_list.append(b)
            c_list.append(c)
filtered_data = pd.DataFrame({'tof': a_list, 'sweeps': b_list, 'counts': c_list})
del a_list, b_list, c_list

#ploting the taken peak on X projection. Ploting the box plot to see how much assymetric is your peak
x_proj_filtered = filtered_data.sort_values(by=['tof']).groupby(['tof']).aggregate(np.sum)
x_proj_filtered=x_proj_filtered.reset_index()


ax_x.bar(x_proj_filtered.tof, x_proj_filtered.counts, color='red')
filtered_data.plot(x='tof', y='sweeps', style='o', alpha=0.15, ms=2,ax=ax_0, label='$\pm $'+str(s)+'$\sigma$'+' of data' ,c='red')
plt.tight_layout()
plt.savefig("data_plot.pdf")

fig2, ax2 = plt.subplots()
ax2.boxplot(filtered_data['tof'], vert=False)
ax2.set_xlabel('time of flight (ns)', fontsize=20)
ax2.set_ylabel('', fontsize=20)
ax2.set_title('Box Plot')
plt.tight_layout()
plt.savefig("box_plot.pdf")

results_df_dict = {}
#MLE estimation of a gaussian destribution

def trunc_norm_pdf(x, amplitude, mean, stddev):
    pdf_vals=amplitude*1/(stddev*np.sqrt(2*np.pi)) * np.exp(-((x - mean) / 4 / stddev)**2)
    return pdf_vals

def log_lik_norm(params, *args):
    amplitude,mean, stddev=params
    x,y=args
    pdf_vals = amplitude*1/(stddev*np.sqrt(2*np.pi)) * np.exp(-((x - mean) / 4 / stddev)**2)
    ln_pdf_vals = y*np.log(pdf_vals)-pdf_vals #main magic happens here
    log_lik_val = ln_pdf_vals.sum()
    neg_log_lik_val = -log_lik_val
    return neg_log_lik_val

a_list = []
b_list = []
c_list = []
d_list = []
e_list = []
mu=peak_pos

for i in range(len(x_proj_peak_ind)):
    results_df_dict['peak_{}'.format(i+1)] = []
    for fit_range in fit_ranges:
        # applying fit range by cutting the data frame
        fit_data = filtered_data.tof[(filtered_data.tof < mu + fit_range) &(filtered_data.tof > mu - fit_range)]
        fit_data_counts = filtered_data.counts[(filtered_data.tof < mu + fit_range) &(filtered_data.tof > mu - fit_range)]
        fit_df=pd.concat([fit_data ,fit_data_counts],join='inner', axis=1)
        fit_df=fit_df.sort_values(by=['tof']).groupby(['tof']).aggregate(np.sum)
        tof=fit_df.index.astype(float)#*1e-6
        counts=fit_df['counts'].astype(float)     
        params_init = np.array([A, mu, sig])
        mle_args = (tof,counts)
        results_uncstr_norm = opt.minimize(log_lik_norm, params_init, args=(mle_args))
        A_MLE,mu_MLE, sig_MLE = results_uncstr_norm.x
        vcv_mle = results_uncstr_norm.hess_inv
        stderr_mu_MLE = np.sqrt(vcv_mle[1,1])
        stderr_sig_MLE = np.sqrt(vcv_mle[2,2])
        #print(A_MLE,'mu_MLE=', mu_MLE, ' sig_MLE=', sig_MLE)
        fig3, ax3 = plt.subplots()      
        hist_gaus_fit = trunc_norm_pdf(tof, A_MLE,mu_MLE,sig_MLE)
        ax3.plot(tof, hist_gaus_fit, color='red', label='fit')
        ax3.bar(tof, counts, label='data')
        ax3.set_ylabel('counts', fontsize=24)
        ax3.set_xlabel('time of flight (ns)', fontsize=24)
#        r = (np.abs(counts - hist_gaus_fit)-0.5)**2 #for small count rates (less then 5 ions per tof). It is so called Yates's correction
        r = (counts - hist_gaus_fit)**2 #for normal count rates
        chisq = np.sum(r/counts)/(len(counts)-1)
        title='$\chi^2$= '+str(round(chisq,3)) +'\n$ \mu$='+str(round(mu_MLE,3))+'$\pm$'+str(round(stderr_mu_MLE,3))+'\n$ \sigma$='+str(round(sig_MLE,3))+'$\pm$'+str(round(stderr_sig_MLE,3))
        ax3.set_title(title)
        ax3.legend()
        a=mu_MLE
        b=stderr_mu_MLE
        c=sig_MLE
        d=stderr_sig_MLE
        e=chisq
        a_list.append(a)
        b_list.append(b)
        c_list.append(c)
        d_list.append(d)
        e_list.append(e)
        plt.tight_layout()
        plt.savefig("gaus_MLE" + str(fit_range) + ".pdf")
        plt.cla()

result_gaus = pd.DataFrame({'chi': e_list,'err_sigma': d_list, 'sigma': c_list,'mu': a_list, 'err_mu': b_list,  })
result_gaus.index=fit_ranges
del a_list, b_list, c_list, d_list
fig4, ax4 = plt.subplots(3)
ax4[0].errorbar(x=fit_ranges, y=result_gaus['sigma'],yerr=result_gaus['err_sigma'], fmt='o', elinewidth=2,color='black', label='gaus')
#ax4[0].set_ylim(min(result_gaus['sigma'])-5,max(result_gaus['sigma'])+5)
ax4[1].errorbar(x=fit_ranges, y=result_gaus['mu'],yerr=result_gaus['err_mu'], fmt='o', elinewidth=2,color='black', label='gaus')
#ax4[1].set_ylim(min(result_gaus['mu'])-5,max(result_gaus['mu'])+5)
ax4[2].scatter(x=fit_ranges, y=result_gaus['chi'],color='black', label='gaus')
#ax4[2].set_ylim(min(result_gaus['chi']),max(result_gaus['chi']))
#ax4[0].set_xlabel('range, ns', fontsize=14, fontweight='bold')
ax4[0].set_ylabel('$\sigma$, ns', fontsize=14)
#ax4[1].set_xlabel('range, ns', fontsize=14, fontweight='bold')
ax4[1].set_ylabel('$\mu$, ns', fontsize=14)
ax4[2].set_xlabel('range, ns', fontsize=14)
ax4[2].set_ylabel('$\chi^2$', fontsize=14)
ax4[0].tick_params( axis='x', bottom=True,  top=True, labelbottom=False)
ax4[1].tick_params( axis='x', bottom=True,  top=True, labelbottom=False)
ax4[2].tick_params( axis='x', bottom=True,  top=True, labelbottom=True)
plt.tight_layout()

plt.tight_layout()
plt.savefig("gaus_parameters_MLE.pdf")
result_gaus['chi_test']=np.abs(result_gaus['chi']-1)
best_range=result_gaus[result_gaus['chi_test'] == result_gaus['chi_test'].min()]

plt.show()
