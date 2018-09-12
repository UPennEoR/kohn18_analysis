#!/bin/sh

#  SaulPaperBandpassPlot.py
#  
#  


import numpy as np, matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import glob
from matplotlib import gridspec

#%%

npzfile = glob.glob("*.B.cal.npz")

for npz in npzfile:

    d = np.load(npz)

    NOT_REAL_ANTS="0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 73, 74, 75, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, 90, 91, 92, 93, 94, 95, 98, 99, 100, 101, 102, 103, 106, 107, 108, 109, 110, 111"#HERA-19

    ex_ants = [int(a) for a in NOT_REAL_ANTS.split(', ')] # Turns str to list
    keep_ants = [a for a in range(d['bp'].shape[2]) if a not in ex_ants] # Equivalen to having an if no statement in for loop. You are looping over a list.
    good_ants = [a for a in keep_ants if a not in [81,22,43,80]]

    # manually flagged RFI channels based on Saul's memos
    low = list(range(0,95))
    orbcomm = list(range(375,390))
    others = list(range(695,705))+list(range(759,761))+list(range(768,771))+list(range(830,832))+list(range(850,852))+[912,932,933,934,935]+[960,961]+list(range(990,1023))
    others+=[510,511,512,850,851,852,853,854]
    msk=low+orbcomm+others

    msk_spec = np.zeros((1024))

    good_chans = [c for c in range(1024) if c not in msk]
    freqs = np.linspace(100,200,num=1024)

    for c in range(1024):
        if c in good_chans:
            msk_spec[c] = 1

#--------------------------------------------------
# Create Amplitude East and North Plots vs Freq
#--------------------------------------------------
# East -Pol
    f,ax = plt.subplots(figsize=(15,8))
    ax.fill_between(freqs[201:300],0,3, facecolor='k',alpha=0.2)
    ax.fill_between(freqs[581:680],0,3,facecolor='k',alpha=0.2)

    for i,a in enumerate(good_ants):
        if a==80:
            continue
        if i <= 8:
            ls='-'
        else:
            ls='--'
        msk_gain = np.ma.masked_where(np.abs(d['bp'][0,:,a])*msk_spec==0.,
                                      np.abs(d['bp'][0,:,a])*msk_spec) #forcing python to not plot RFI
        index = np.where(msk_gain >0.5)[0]
        ax.plot(freqs,msk_gain,ls,label=str(a),lw=2)

    ax.set_ylim(0,3)
    ax.set_xlim(115,188)
    plt.legend(loc='upper right')
    plt.xlabel('Frequency [MHz]',size=12)
    plt.ylabel('Amplitude [arb.]',size=12)
    plt.grid()
    f.suptitle("Amplitude:East-Pol Gain Solutions", fontsize=14)#Title centered above all subplots
    f.savefig('amp_east_{}.png'.format(npz))
#    plt.show()
    plt.close()

# North -Pol
    f,ax = plt.subplots(figsize=(15,8))
    ax.fill_between(freqs[201:300],0,3, facecolor='k',alpha=0.2)
    ax.fill_between(freqs[581:680],0,3,facecolor='k',alpha=0.2)

    for i,a in enumerate(good_ants):
        if a==80:
            continue
            if i <= 8:
                ls='-'
        else:
            ls='--'
            msk_gain = np.ma.masked_where(np.abs(d['bp'][1,:,a])*msk_spec==0.,
                                          np.abs(d['bp'][1,:,a])*msk_spec) #forcing python to not plot RFI
                
            index = np.where(msk_gain >0.5)[0]
            ax.plot(freqs,msk_gain,ls,label=str(a),lw=2)

    ax.set_ylim(0,3)
    ax.set_xlim(115,188)
    plt.legend(loc='upper right')
    plt.xlabel('Frequency [MHz]',size=12)
    plt.ylabel('Amplitude [arb.]',size=12)
    plt.grid()
    f.suptitle("Amplitude:North-Pol Gain Solutions", fontsize=14)#Title centered above all subplots
    f.savefig('amp_north_{}.png'.format(npz))
#    plt.show()
    plt.close()

#--------------------------------------------------
# Create Phase East and North Plots vs Freq
#--------------------------------------------------
# East -Pol
    f,ax = plt.subplots(figsize=(15,8))
    
    for i in range(len(good_ants)):
        if i < 10:
            ax.plot(freqs,np.angle(d['bp'][0,:,good_ants[i]],deg=True),'--',label=str(good_ants[i]))
        else:
            ax.plot(freqs,np.angle(d['bp'][0,:,good_ants[i]],deg=True),'.',label=str(good_ants[i]))

    ax.set_xlim(115,188)
    plt.xlabel('Frequency [MHz]',size=12)
    plt.ylabel('Phase [deg]')
    plt.legend(loc='upper right')
    f.suptitle("Phase:East-Pol Gain Solutions")
    f.savefig("phase_east_{}.png".format(npz))
#    plt.show()
    plt.close()

# North -Pol
    f,ax = plt.subplots(figsize=(15,8))

    for i in range(len(good_ants)):
        if i < 10:
            ax.plot(freqs,np.angle(d['bp'][1,:,good_ants[i]],deg=True),'--',label=str(good_ants[i]))
        else:
            ax.plot(freqs,np.angle(d['bp'][1,:,good_ants[i]],deg=True),'.',label=str(good_ants[i]))

    ax.set_xlim(115,188)
    plt.xlabel('Frequency [MHz]',size=12)
    plt.ylabel('Phase [deg]')
    plt.legend(loc='upper right')
    f.suptitle("Phase:North-Pol Gain Solutions")
    f.savefig("phase_north_{}.png".format(npz))
#    plt.show()
    plt.close()
