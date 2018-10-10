import numpy as np, matplotlib.pyplot as plt

d = np.load('2457548.45923.npz')

NOT_REAL_ANTS="0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 73, 74, 75, 76, 77, 78, 79, 82, 83, 84, 85, 86, 87, 90, 91, 92, 93, 94, 95, 98, 99, 100, 101, 102, 103, 106, 107, 108, 109, 110, 111"

ex_ants = [int(a) for a in NOT_REAL_ANTS.split(', ')]
keep_ants = [a for a in range(d['gains'].shape[2]) if a not in ex_ants]
good_ants = [a for a in keep_ants if a not in [81,22,43]]

low = list(range(0,102))
orbcomm = list(range(375,390))
others = list(range(695,705))+list(range(759,761))+list(range(768,771))+list(range(830,832))+list(range(850,852))+list(range(921,1023))
others+=[510,511,512,850,851,852,853,854]
msk=low+orbcomm+others

msk_spec = np.zeros((1024))

good_chans = [c for c in range(1024) if c not in msk]
freqs = np.linspace(100,200,num=1024)

for c in range(1024):
    if c in good_chans:
        msk_spec[c] = 1

f,ax = plt.subplots()
ax.fill_between(freqs[202:301],0,3, facecolor='k',alpha=0.2)
ax.fill_between(freqs[581:681],0,3,facecolor='k',alpha=0.2)

for i,a in enumerate(good_ants):
    if a==80:
        continue
    if i <= 8:
        ls='-'
    else:
        ls='--'
    msk_gain = np.ma.masked_where(np.abs(d['gains'][0,:,a])*msk_spec==0.,
                                  d['gains'][0,:,a])
    
    
    #ax.plot(np.fft.fftshift(np.fft.fftfreq(1024,np.diff(freqs)[0]*1e6))*1e9,
    #        np.abs(np.fft.fftshift(np.fft.ifft(msk_gain))),
    #        ls,label=str(a),lw=2)
    ax.plot(freqs,np.abs(msk_gain),ls,label=str(a),lw=2)

ax.set_ylim(0,2.75)
#ax.set_ylim(-np.pi,np.pi)
#ax.set_xlim(100,200)
#plt.legend()
#plt.xlabel('Frequency [MHz]',size=12)
plt.ylabel('Amplitude [arb.]',size=12)
plt.grid()

plt.show()
