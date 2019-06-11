from pyseeman.techlib import ITRS16cap, ITRS16sw
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
ITRS16cap.rating = 3
ITRS16sw.drain_rating = 3
ITRS16sw.gate_rating = 3

vin = np.linspace(1, 5, 1000).reshape((1, 1000))
plot_regulation(topologies, vin, 1.2, 100e-3, 1e-6, [ITRS16sw], [ITRS16cap])

plt.show()
