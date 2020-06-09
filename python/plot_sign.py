import os
import sys
import numpy as np
from operator import itemgetter

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.getcwd() + "/../")
from combine_mc_data import combine_mc_data
from bootstrap import bootstrap


####
#### what data to fetch
####
jobNames = [
  # N=1
  "saved/19-12-03",
  "saved/19-12-04",
  "saved/19-12-11",
  "saved/19-12-12",
  "saved/19-12-17",
  "saved/19-12-22",
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

errorMarker = {
  "marker": "o",
  "ls": ":",
  "capsize": 4,
  "lw": 1,
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
#### loop over data
####
sign = []
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
  #### sign
  ####
  if N == 1:
    val, err = bootstrap(lambda x: (x), [np.array(dataSet["avgSign"])], np.array(dataSet["avgSign"])*0 + 1, dataSet["numData"], maxNumBootstrapIteration)
    sign += [{
      "beta": beta,
      "N":    N,
      "val":  val,
      "err":  err
    }]


####
#### plot average sign
####
plt.title(r"$ \langle s \rangle_\mathrm{B} $")
plt.xlabel(r"$ \beta t $")

# extract betas
beta = [f['beta'] for f in sorted(sign, key=itemgetter("beta")) if f["N"] != 0]

S_val = np.array([f['val'] for f in sorted(sign, key=itemgetter("beta"))])
S_err = np.array([f['err'] for f in sorted(sign, key=itemgetter("beta"))])
plt.errorbar(np.array(beta) * abs(t), S_val, yerr=S_err, **errorMarker, zorder=2)

plt.yscale("log")
plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter("{x}"))


####
#### indicate beta = 22, 88
####
for _beta in [22, 88]:
  index = beta.index(_beta)
  plt.plot(_beta * abs(t), S_val[index], **indicator, zorder=1)

plt.tight_layout()
plt.show()