#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
@file Poincare_1.py
@author Yasuhiro Suzuki, National Institute for Fusion Science
@brief A script that makes Poincare plot
@date 20 Apr 2020
@version Initial Version
"""

import math, numpy, sys, json
import matplotlib.pyplot as plt

json_open = open('setting.json', 'r')

json_load = json.load(json_open)

mode = json_load['settings']['mode']
lmin = json_load['settings']['lc_min']
lmax = json_load['settings']['lc_max']
rmin = json_load['dimensions']['rmin']
rmax = json_load['dimensions']['rmax']
zmin = json_load['dimensions']['zmin']
zmax = json_load['dimensions']['zmax']

jobid = sys.argv[1]

fig = plt.figure()

rfile =  jobid + ".plot"

R, Z, Lc = numpy.loadtxt(rfile, delimiter=',', usecols=(0,1,2), unpack=True)

if mode == "2D":

    plt.plot(R, Z, color="red", marker=".", markersize=0.1, linewidth=0.0)

elif mode == "3D":

    plt.scatter(R, Z, marker="o", c=numpy.log10(Lc), s=0.5, vmin=math.log10(lmin), vmax=math.log10(lmax), linewidth=0.0, cmap=plt.cm.jet)
    cb = plt.colorbar()
    cb.set_label(u"$\\mathrm{log}(L_C)$ ", size=15)

plt.xlim(rmin,rmax)
plt.ylim(zmin,zmax)
plt.axes().set_aspect(1)
plt.legend(loc='lower left')

plt.xlabel(u'$R$ [m]', size=15)
plt.ylabel(u'$Z$ [m]', size=15)

plt.title(rfile)

plt.grid(True)

#plt.savefig(jobid + ".png")

plt.tight_layout()
plt.show()
