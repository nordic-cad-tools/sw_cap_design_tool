import numpy as np

class Implementation:

    def __init__(self, *args, **kwargs):
        pass

    def __repr__(self):
        my_string = ""
        for k, v in self.__dict__.items():
            my_string += f"{k} = {v}\n"
        return my_string

def implement_topology(topology, Vin, switchTechs, capTechs, compMetric=1):
    """
    implement_topology: matches components with switches and components in
    the topology.
    implementation = implement_topology(topology, Vin, switchTechs, capTechs, compMetric)
       topology: structure created by generate_topology
       Vin: input voltage of converter
       switchTechs: an array of technology structures available for switch use
       capTechs: an array of technology structures available for cap use
       compMetric: a metric (1=area, 2=loss) used for determining the best
                   component (1=default)
    """

    #Break out components of topology structure
    ratio = topology.ratio
    ac = topology.ac
    ar = topology.ar
    vc = topology.vc * Vin
    vr = topology.vr * Vin
    vcb = topology.vcb * Vin
    vrb = topology.vrb * Vin

    switch_assign = []
    cap_assign = []
    switch_rel_size = []
    cap_rel_size = []

    # Assign Capacitors
    for i in range(ac.shape[0]):
        print(i)
        Mc = 0
        Cc = 0 # cap cost
        for j, cap_tech in enumerate(capTechs):
            if vc[i] <= cap_tech.rating:
                # Cap could work ... let's see if it's good
                # Use area-limited metric, which is usually applicable
                C = cap_tech.area
                M = cap_tech.capacitance*vc[i]**2 / C
                if M > Mc:
                    if Mc == 0:
                        cap_assign.append(cap_tech)
                    else:
                        cap_assign.append(cap_tech)

                    Mc = M
                    Cc = C
        # check to make sure a suitable device exists
        if Mc == 0:
            raise ValueError("No capacitors meet the voltage requirement of: {}".format(vc[i]))

        # determine relative device size
        if ac[i] == 0:
            cap_rel_size = cap_rel_size.append(0)  # avoid divide by 0
        else:
            cap_rel_size.append((ac[i]*vc[i]) / (np.sqrt(Mc)*cap_assign[i].area))


    # Assign Switches
    for i in range(ar.shape[0]):

        Msw = 0
        Csw = 0  # switch cost;
        for j, switch_tech in enumerate(switchTechs):
            if (vr[i] <= switch_tech.drain_rating):
                # Switch could work let's see if it's good
                if compMetric == 2: # loss metric
                    # assume full gate drive
                    C = switch_tech.gate_cap * switch_tech.gate_rating ** 2 + switch_tech.drain_cap * vr[i] ** 2 + \
                        switch_tech.body_cap * vrb[i] ** 2
                    M = switch_tech.conductance * vr[i] ** 2 / C
                else: #  area metric
                    C = switch_tech.area
                    M = switch_tech.conductance * vr[i] ** 2 / C

            if M > Msw:
                if Msw == 0:
                    switch_assign.append(switch_tech)
                else:
                    switch_assign[i] = switch_tech
                Msw = M
                Csw = C

        # check to make sure a suitable device exists
        if Msw == 0:
            raise ValueError("No switches meet the voltage requirement of: {}".format(vr[i]))

        # determine relative device size
        if ar[i] == 0:
            switch_rel_size.append(0)
        else:
            if compMetric == 2:
                switch_rel_size.append(ar[i]*vr[i] / (np.sqrt(Msw)*switch_assign[i].conductance))
            else:
                switch_rel_size.append(ar[i]*vr[i] / (np.sqrt(Msw)*switch_assign[i].area))

    # Scale Caps for unit area:

    cap_area = 0;
    for i in range(ac.shape[0]):
        cap_area = cap_area + cap_rel_size[i]*cap_assign[i].area

    cap_size = cap_rel_size/(cap_area+1e-30)


    # Scale Switches for unit area:
    sw_area = 0
    #print(switch_rel_size)
    for i in range(ar.shape[0]):
     #   print(i)
        sw_area = sw_area + switch_rel_size[i]*switch_assign[i].area

    #if sw_area > 0:
    #print(sw_area)
    switch_size = switch_rel_size / sw_area


    # Create implementation structure
    implementation = Implementation()
    implementation.topology = topology
    implementation.capacitors = cap_assign
    implementation.switches = switch_assign
    implementation.cap_size = cap_size
    implementation.switch_size = switch_size

    return implementation

if __name__ == '__main__':
    from pyseeman.generate_topology import generate_topology
    from pyseeman.techlib import ITRS16cap, ITRS16sw
    topo = generate_topology("series-parallel", 1, 3)

    implementation = implement_topology(topo, Vin=1, switchTechs=[ITRS16sw], capTechs=[ITRS16cap], compMetric=1)
    print(implementation)
