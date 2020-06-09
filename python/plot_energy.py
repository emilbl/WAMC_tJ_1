import os
import sys
from glob import glob
import numpy as np
from operator import itemgetter

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as mticker

sys.path.insert(0, os.getcwd() + "/../")
from combine_mc_data import combine_mc_data
from bootstrap import bootstrap


####
#### data from exact diagonalization
####
insetData = {
  "tbetas": [1.1, 1.6500000000000001, 2.2, 2.75, 3.3000000000000003, 4.4, 5.5, 6.6000000000000005, 7.7, 8.8, 10.0, 12.0, 14.0, 16.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0],
  "tKs": [-2.5793065429657607, -2.830590331017092, -2.919771863200312, -2.9502514424098094, -2.95602095748353, -2.93959056533772, -2.9140912335903537, -2.8920510280845266, -2.87631409186938, -2.8660837624666624, -2.8593381554866495, -2.8540051009245597, -2.8522605757595545, -2.8519833528010974, -2.852743388338861, -2.8546578603080235, -2.855353339555331, -2.8555590559060007, -2.8556172723243254, -2.8556335268715918, -2.8556380446388188]
}


####
####
####
jobNames = [
  # N=1
  "saved/19-12-03",
  "saved/19-12-04",
  "saved/19-12-11",
  "saved/19-12-12",
  "saved/19-12-17",
  "saved/19-12-22",

  # N=0
  "saved/20-04-07"
]

filter = {
  "L": 20,
  "t": -0.1,
  "J": 0.03
}


####
#### application parameters
####
maxNumBootstrapIteration = 1000

style = {
  "marker":    "o",
  "ms":        10,
  "linestyle": ":"
}

marker = {
  "marker":    "o",
  "color":     colorRed,
  "ms":        5,
  "linestyle": ":",
  "lw":        1,
}

errorMarker = {
  "marker":     "o",
  "ls":         ":",
  "capsize":    4,
  "lw":         1,
  "elinewidth": 1,
  "markersize": 5
}

indicator = {
  "marker": "o",
  "ms":     20,
  "color":  [0.7, 0.7, 0.7]
}


####
#### collect and sort data
####
dataSets = combine_mc_data(jobNames, filter)
dataSets = list(dataSets.values())
dataSets.sort(key=lambda data: data["beta"])


####
#### create figure
####
fig = plt.figure(figsize=(9, 6))
plt.rc('font', size=18)


####
#### extract colors from default color cycle
####
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
colorRed = colors[3]


####
#### loop over data
####
Et      = []
EJ_nn_0 = []
EJ_0    = []
for i, dataSet in enumerate(sorted(dataSets, key=itemgetter("t", "beta"))):
  ####
  #### get parameters
  ####
  beta = dataSet["beta"]
  L    = dataSet["L"]
  t    = dataSet["t"]
  J    = dataSet["J"]
  N    = dataSet["numHoles"]


  ####
  #### kinetic energy
  ####
  if N == 1:
    val, err = bootstrap(lambda x, y: (x + y), [np.array(dataSet["avgExchaEnergy"])[:, 0], np.array(dataSet["avgKinetEnergy"])[:, 1]], dataSet["avgSign"], dataSet["numData"], maxNumBootstrapIteration)
    Et += [{
      "beta": beta,
      "N":    N,
      "val":  val,
      "err":  err
    }]

  ####
  #### NN interaction energy
  ####
  val, err = bootstrap(lambda x: (x), [np.array(dataSet["avgNnInterEnergy"])[:, 0]], dataSet["avgSign"], dataSet["numData"], maxNumBootstrapIteration)
  EJ_nn_0 += [{
    "beta": beta,
    "N":    N,
    "val":  val,
    "err":  err
  }]

  ####
  #### exchange energy
  ####
  val, err = bootstrap(lambda x, y: (x + y), [np.array(dataSet["avgNnInterEnergy"])[:, 0], np.array(dataSet["avgKinetEnergy"])[:, 0]], dataSet["avgSign"], dataSet["numData"], maxNumBootstrapIteration)
  EJ_0 += [{
    "beta": beta,
    "N":    N,
    "val":  val,
    "err":  err
  }]


####
#### plot average sign
####
beta    = [f['beta'] for f in sorted(Et, key=itemgetter("beta")) if f["N"] != 0]
E_K_val = np.array([f['val'] for f in sorted(Et, key=itemgetter("beta"))])
E_K_err = np.array([f['err'] for f in sorted(Et, key=itemgetter("beta"))])
plt.errorbar(np.array(beta) * abs(t), E_K_val / abs(t), yerr=E_K_err / abs(t), **errorMarker, zorder=2)

plt.xscale("log")
plt.gca().xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:.0f}"))
plt.gca().xaxis.set_minor_formatter(mticker.StrMethodFormatter("{x:.0f}"))
plt.xlabel(r"$ \beta t $")
plt.ylabel(r"$ \langle E_\mathrm{K} \rangle / t $")

####
#### indicate beta = 22, 88
####
for _beta in [22, 88]:
  index = beta.index(_beta)
  plt.plot(_beta * abs(t), E_K_val[index] / abs(t), **indicator, zorder=1)

# exact data inset
ins = inset_axes(plt.gca(), width="60%", height="60%", loc=1)
ins.plot(insetData["tbetas"], insetData["tKs"], **marker, zorder=2)

# indicate saturation value
ins.hlines(insetData["tKs"][-1], 0.9, 99, linestyle="--", color=(0.8, 0.8, 0.8), zorder=1)
ins.set_xlim(left=0.9, right=99)

ins.set_xscale("log")
ins.xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:.0f}"))
ins.xaxis.set_minor_formatter(mticker.NullFormatter())
ins.set_xlabel(r"$ \beta t $")
ins.set_ylabel(r"$ \langle E_\mathrm{K} \rangle / t $")


plt.tight_layout()
plt.show()