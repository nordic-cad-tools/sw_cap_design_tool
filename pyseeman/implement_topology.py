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
    for i in range(ac.shape[1]):
        print(i)
        Mc = 0
        Cc = 0 # cap cost
        for j, cap_tech in enumerate(capTechs):
            if vc(i) <= cap_tech:
                # Cap could work ... let's see if it's good
                # Use area-limited metric, which is usually applicable
                C = cap_tech.area
                M = cap_tech.capacitance*vc(i)**2 / C
                if (M > Mc):
                    if (Mc == 0):
                        cap_assign = [cap_assign, cap_tech]
                    else:
                        cap_assign.append(cap_tech)

                    Mc = M
                    Cc = C
        # check to make sure a suitable device exists
        if Mc == 0:
            raise ValueError("No capacitors meet the voltage requirement of: {}".format(vc(i)))


if __name__ == '__main__':
    pass
