# Worm-algorithm Monte Carlo

Worm-algorithm Monte Carlo applied to the t-J model for a single carrier

This code is accompanied with the publication "Unbiased description of magnetic polarons in a Mott insulator", Emil Blomquist and Johan Carlström, arXiv:1912.08825 (to be published in Nature Communications Physics)


## Installation

First clone the repository
```
git clone https://github.com/emilbl/WAMC_tJ_1.git
```

Then fetching all required submodules
```
git submodule update --init --recursive
```

Finally compile the code (requires at least `gcc v. 7.5.0`)
```
make
```

## Run simulation

```
./worm.sh -i data/input/tJ.json
```
The output will be stored in `data/jobs`

## Analyze data

The data analysis scripts require python 3. Also, make sure you have `numpy` and `matplotlib` installed (`pip3 install scipy`)

```
python3 python/plot_sign.py
```

```
python3 python/plot_energy.py
```

```
python3 python/plot_spin_correlations.py
```