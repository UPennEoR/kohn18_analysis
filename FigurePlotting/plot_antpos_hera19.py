#! /usr/bin/env python
import aipy as a, numpy as n, matplotlib.pyplot as p #pylab as p
import sys, optparse
import IPython
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', family='serif')

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--aspect_neq', action='store_true', help='Do not force equal aspect in x/y plot.')
o.add_option('--nonumbers', action='store_true', help='Do not plot antenna numbers next to symbols.')
o.add_option('--ex_ants', default=[], help='Plot a red "X" over these antennas.')
o.add_option('--dishborders', action='store_true', help="Plot dish borders.")
o.add_option('--uvcover',action='store_true', help='Plot uv coverage also')
o.add_option('--uvcover_hilo',action='store_true', help='Plot uv coverage for high and low bands as well.')
opts, args = o.parse_args(sys.argv[1:])

th = n.arange(0, 2*n.pi, .01)
r = 5.
if len(opts.ex_ants)>0:
    ex_ants = map(int,opts.ex_ants.split(','))
else:
    ex_ants = []

aa = a.cal.get_aa(opts.cal, .1, .1, 1)
antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
antpos = n.array(antpos) * a.const.len_ns / 100.
x,y,z = antpos[:,0], antpos[:,1], antpos[:,2]
x -= n.average(x)
y -= n.average(y)

for ant,(xa,ya,za) in enumerate(zip(x,y,z)):
    hx,hy = r*za*n.cos(th)+xa, r*za*n.sin(th)+ya
    if za > 0: fmt = '#eeeeee'
    else: fmt = '#a0a0a0'
    #p.fill(hx,hy, fmt)
    if za != 0:
        #p.plot(xa,ya,'k.',ms=12)
        circle = p.Circle((xa,ya), radius=7, fill=False)
        p.gca().add_patch(circle)
        if not opts.nonumbers: p.text(xa,ya, str(ant))
        if opts.dishborders:
            circ = p.Circle((xa, ya), 7.3, color='k', fill=False, lw=3)            
        if ant in ex_ants:
            p.plot(xa,ya,'rx',ms=7)
#p.grid()
p.xlim(-40,40)
p.ylim(-40,40)
p.xlabel("East-West Antenna Position (m)")
p.ylabel("North-South Antenna Position (m)")
a = p.gca()
if not opts.aspect_neq: a.set_aspect('equal')
p.savefig('antpos_hera19.pdf',clobber=True)
#p.show()
p.close()

# XXX make these options
lam = 2 #m
lamlo,lamhi,num = 1,4,1024

lamlo_lo,lamhi_lo,num_lo=2.3,2.5,100
lamlo_hi,lamhi_hi,num_hi=1.8,1.9,100

if opts.uvcover:
    for ant1,(xa1,ya1,za1) in enumerate(zip(x,y,z)):
        if ant1 in ex_ants:
            continue
        if za1< 0:
            continue
            
        for ant2,(xa2,ya2,za2) in enumerate(zip(x,y,z)):
            if ant1==ant2:
                continue
            if ant2 in ex_ants:
                continue
            if za2 < 0:
                continue
            dx = (xa1-xa2)/lam
            dy = (ya1-ya2)/lam
            if dx == 0 and dy == 0:
                continue
            dxarr = dx*n.linspace(lamlo,lamhi,num=num)
            dyarr = dy*n.linspace(lamlo,lamhi,num=num)
            
            dxarr_lo = dx*n.linspace(lamlo_lo,lamhi_lo,num=num_lo)
            dyarr_lo = dy*n.linspace(lamlo_lo,lamhi_lo,num=num_lo)
            
            dxarr_hi = dx*n.linspace(lamlo_hi,lamhi_hi,num=num_hi)
            dyarr_hi = dy*n.linspace(lamlo_hi,lamhi_hi,num=num_hi)
            
            p.plot(dxarr,dyarr,'b.',alpha=0.1)
            
            if opts.uvcover_hilo:
                p.plot(dxarr_lo,dyarr_lo,'r.',alpha=0.1)
                p.plot(dxarr_hi,dyarr_hi,'y.',alpha=0.1)
            
    p.grid()
    p.xlabel(r'u [$\lambda$]',size=15)
    p.ylabel(r'v [$\lambda$]',size=15)
    if not opts.uvcover_hilo:
        p.savefig('uvcoverage_hera19.png',clobber=True)
    else:
        p.savefig('uvcoverage_hera19-hilo.png',clobber=True)
    p.close()


