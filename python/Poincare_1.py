#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
@file Poincare_1.py
@author Yasuhiro Suzuki, National Institute for Fusion Science
@brief
@date 20 Apr 2020
@version Initial Version
"""

import math, re, numpy, sys
import matplotlib.pyplot as plt

"""
@brief This script make Poincare plot
"""

jobid = sys.argv[1]

fig = plt.figure()

rfile =  jobid + ".plot"

R, Z, Lc = numpy.loadtxt(rfile, delimiter=',', usecols=(0,1,2), unpack=True)

maps=[m for m in plt.cm.datad if not m.endswith("_r")]

plt.scatter(R, Z, marker="o", c=numpy.log10(Lc), s=0.5, vmin=1.0, vmax=3.0, linewidth=0.0, cmap=plt.cm.jet)
cb = plt.colorbar()
cb.set_label(u"$\\mathrm{log}(L_C)$ ", size=15)

plt.xlim(2.2,5.6)
plt.ylim(-1.7,1.7)
plt.axes().set_aspect(1)
plt.legend(loc='lower left')

plt.xlabel(u"$R$ [m]", size=15)
plt.ylabel(u"$Z$ [m]", size=15)

plt.title(rfile)

plt.grid(True)

#plt.savefig(jobid + ".png")

plt.tight_layout()
plt.show()
