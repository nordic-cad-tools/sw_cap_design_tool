from pyseeman.techlib import ITRS16cap, ITRS16sw
from pyseeman.core import Topology, plot_opt_contour
import matplotlib.pyplot as plt

# create the ladder topology with a 3:2 ratio
ladder = Topology("Ladder", 2, 3)
# implement the topology
ladder_imp = ladder.implement(
    vin=1.5, switch_techs=[ITRS16sw], cap_techs=[ITRS16cap], comp_metric=1
)
# plot the efficiency contour
# ax1 = plot_opt_contour(ladder_imp, vin=1.5, iout=10e-3, area_cap=1e-6)
ax1 = ladder_imp.plot_opt_contour(iout=10e-3, area_cap=1e-6)
plt.show()
