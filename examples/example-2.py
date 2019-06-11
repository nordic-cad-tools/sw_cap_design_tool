from pyseeman.techlib import ITRS16cap, ITRS16sw, Switch, Capacitor
from pyseeman.core import Topology, plot_regulation
import numpy as np

import matplotlib.pyplot as plt

# create a list of topologies
topologies = [
    Topology("series-parallel", 1, 3),
    Topology("series-parallel", 1, 2),
    Topology("series-parallel", 2, 3),
    Topology("series-parallel", 4, 5),
]

# vout sweep
vout = np.linspace(0, 1.2, 1000).reshape((1, 1000))
plot_regulation(topologies, 1.2, vout, 500e-3, 1e-7, [ITRS16sw], [ITRS16cap])

# vin sweep
hv_switch = Switch(
    tech_name="HV tech",
    dev_name="HV switch",
    area=100e-6 * 100e-6,
    conductance=1e-3,
    gate_rating=5.0,
    drain_rating=5.0,
    gate_cap=10e-15,
    drain_cap=1e-15,
    body_cap=1e-15,
)

hv_cap = Capacitor(
    tech_name="HV tech",
    dev_name="HV cap",
    capacitance=20e-12,
    area=100e-6 * 100e-6,
    bottom_cap=100e-15,
    esr=0,
    rating=5,
)


vin = np.linspace(1, 5, 1000).reshape((1, 1000))
plot_regulation(topologies, vin, 1.2, 100e-3, 1e-6, [hv_switch], [hv_cap])

plt.show()
