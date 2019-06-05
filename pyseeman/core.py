import numpy as np

def fibfun(n):
    if n<0:
        raise ValueError("Invalid argument for fibfun")
    elif n == 1:
        return 0
    elif n == 2:
        return 1
    else:
        return fibfun(n-1) + fibfun(n-2)

class Topology:

    _valid_names = {"ladder", "dickson", "cockcroft-walton",
                    "doubler", "series-parallel", "fibonacci"}

    def __init__(self, name, num, den):
        if name not in self._valid_names:
                raise ValueError(f"{name} is not a valid name {self._valid_names}")
        self.name = name
        self.num = num
        self.den = den
        self._generate_topology()

    def _generate_topology(self):
        # If input ratio is step-down, flip to step up for ac computation
        # and then flip it again before returning result.
        flip = 0
        if self.num / self.den < 1:
            flip = 1
            den, num = self.num, self.den
        else:
            num, den = self.num, self.den

        # *************************** SERIES - PARALLEL *******************************#
        if self.name.lower() == 'series-parallel':
            n = num
            m = den
            N = n

            # SSL values
            ac = np.ones(m * (n - m) // m)
            vc = np.ones(m * (n - m) // m)
            vcb = np.array([])

            for i in range(1, m + 1):
                for j in range(1, n - m + 1):
                    vcb = np.append(vcb, (i + j - 1) / m)

            # FSL values
            vr = np.array([])
            vrb = np.array([])
            for i in range(1, m + 1):
                for j in range(1, n - m + 2):
                    if j == 1:
                        vr = np.append(vr, i / m)
                        vrb = np.append(vrb, (i + j - 1) / m)
                    elif j == n - m + 1:
                        vr = np.append(vr, (n - m - 1 + i) / m)
                        vrb = np.append(vrb, (i + j - 2) / m)
                    else:
                        vr = np.append(vr, 1 / m)
                        vrb = np.append(vrb, (i + j - 1) / m)
            for i in range(1, m + 2):
                for j in range(1, n - m + 1):
                    if i == 1:
                        vr = np.append(vr, j / m)
                    elif i == m + 1:
                        vr = np.append(vr, (m - 1 + j) / m)
                    else:
                        vr = np.append(vr, 1 / m)
                    if i == 1 or i == m + 1:
                        vrb = np.append(vrb, 0)
                    else:
                        vrb = np.append(vrb, (i + j - 2) / m)
            ar = np.ones(len(vr)) / m

        elif self.name.lower() == 'ladder':
            n = num
            m = den
            N = n

            ac = np.array([])
            vc = np.array([])
            vcb = np.array([])

            for j in range(1, n - m):
                ac = np.append(ac, [j, j])
                vc = np.append(vc, [1 / m, 1 / m])
                vcb = np.append(vcb, [1 / m, 0])

            ac = np.append(ac, n - m)
            vc = np.append(vc, 1 / m)
            vcb = np.append(vcb, 1 / m)

            for j in range(m - 1, 0, -1):
                ac = np.append(ac, np.ones(2) * j * (n / m - 1))
                vc = np.append(vc, [1 / m, 1 / m])
                vcb = np.append(vcb, [1 / m, 0])

            ar = np.hstack([np.ones(2 * (n - m)), np.ones(2 * m) * (n / m - 1)])
            vr = np.ones(2 * n) / m
            vrb = np.mod(np.linspace(0, (2 * n - 1), 2 * n), 2) / m

        elif self.name.lower() == "dickson":

            if den != 1:
                raise ValueError('SWITCHCAP:nonIntegerRatio the Dickson topology supports integer ratios only')

            N = num

            # SSL values
            ac = np.ones(N - 1)
            vc = np.array([])
            vcb = np.ones(N - 1)

            for j in range(1, N):
                vc = np.append(vc, j)

            if N == 2:
                vr = np.ones(4)
                ar = np.ones(4)
                vrb = np.array([0, 1, 0, 1])
            else:
                vr = np.hstack([np.ones(5), 2 * np.ones(N - 2), 2])
                ar = np.hstack([np.ones(2) * np.floor((j + 1) / 2.), np.ones(2) * np.floor(j / 2.), np.ones(N)])
                vrb = np.hstack([np.array([0, 1, 0, 1]), np.ones(N)])

        elif self.name.lower() == "cockcroft-walton":

            if den != 1:
                raise ValueError('SWITCHCAP:nonIntegerRatio the Cockcroft-Walton topology supports integer ratios only')

            N = num
            # SSL values
            ac = np.array([])
            vc = np.hstack([1, np.ones(N - 2) * 2])
            vcb = np.ones(N - 1)

            for j in range(1, N):
                ac = np.append(np.floor((j + 1) / 2.), ac)

            # FSL values
            if N == 2:
                vr = np.ones(4)
                ar = np.ones(4)
                vrb = np.array([0, 1, 0, 1])
            else:
                vr = np.hstack([np.ones(5), np.ones(N - 2) * 2, 1])
                ar = np.hstack([np.floor((j + 1) / 2.) * np.ones(2), np.floor(j / 2.) * np.ones(2), np.ones(N)])
                vrb = np.hstack([np.array([0, 1, 0, 1]), np.ones(N)])

        elif self.name.lower() == "doubler":

            if den != 1:
                raise ValueError('SWITCHCAP:nonIntegerRatio the Doubler topology supports integer ratios only')

            n = np.ceil(np.log2(num)).astype(np.int)
            N = 2 ** n
            if N != num:
                raise ValueError('SWITCHCAP:badRatio the Doubler topology supports conversion ratios ~ 2^n')

            # SSL values
            ac = np.array([])
            vc = np.array([])
            vcb = np.array([])

            for j in range(1, 2 * n):
                ac = np.append(2 ** np.floor((j - 1) / 2.), ac)
                vc = np.append(vc, 2 ** np.floor(j / 2.))
                vcb = np.append(vcb, 2 ** np.floor(j / 2.) * np.mod(j, 2))

            # FSL values
            ar = np.array([])
            vr = np.array([])
            vrb = np.array([])

            for j in range(1, n + 1):
                ar = np.append(ar, np.ones(4) * 2 ** (j - 1))
                vr = np.append(vr, np.ones(4) * 2 ** (n - j))
                vrb = np.append(vrb, np.array([0, 1, 0, 1]) * 2 ** (n - j))

        elif self.name.lower() == "fibonacci":
            if den != 1:
                raise ValueError('SWITCHCAP:nonIntegerRatio the Fibonacci topology supports integer ratios only')

            i = 2

            while fibfun(i) < num:
                i = i + 1

            if fibfun(i) > num:
                raise ValueError('SWITCHCAP:badRatio the fibonacci topology supports ratios of F_n or 1/F_n only')
            N = fibfun(i)

            ac = np.array([])
            vc = np.array([])
            vcb = np.array([])

            for j in range(2, i):
                ac = np.append(fibfun(j - 1), ac)
                vc = np.append(vc, fibfun(j - 1))
                vcb = np.append(vcb, fibfun(j - 1))

            ar = np.array([1])
            vr = np.array([])
            vrb = np.array([0])

            for j in range(2, i):
                ar = np.append(np.array([fibfun(j), fibfun(j - 1), fibfun(j - 1)]), ar)
                vr = np.append(vr, np.array([fibfun(j), fibfun(j), fibfun(j - 1)]))
                vrb = np.append(vrb, np.array([fibfun(j - 1), 0, fibfun(j - 1)]))

            vr = np.append(vr, fibfun(i - 1))
        else:
            raise ValueError("Topology type not implemented yet")

        # TODO: Check if it makes sense that M values are claulcated before the flipping
        ratio = num / den
        Mssl = 2 * ratio ** 2 / np.sum(ac * vc) ** 2
        Mfsl = ratio ** 2 / (2 * np.sum(ar * vr) ** 2)

        if flip == 1:
            ac = ac / ratio
            vc = vc / ratio
            vcb = vcb / ratio
            ar = ar / ratio
            vr = vr / ratio
            vrb = vrb / ratio
            ratio = 1 / ratio

        self.ac = ac
        self.vc = vc
        self.vcb = vcb
        self.ar = ar
        self.vr = vr
        self.vrb = vrb
        self.Mssl = Mssl
        self.Mfsl = Mfsl
        self.ratio = ratio

    def implement(self, vin, switch_techs, cap_techs, comp_metric=1):
        return Implementation(self, vin, switch_techs, cap_techs, comp_metric)

    def permute(self, topology):
        """
        Returns a new topology consisting of every permutation of topologies. The specified topology will be cascaded
        with the current topology.
        :param topology:
        :return:
        """
        pass

    def __repr__(self):
        return f"<Topology(name={self.name}, ratio={self.ratio})>"


class Implementation:

    def __init__(self, topology, vin, switch_techs, cap_techs, comp_metric=1):
        """
        :param topology:  matches components with switches and components in the topology.
        :param vin: input voltage of converter
        :param switch_techs:  a list of technology structures available for switch use
        :param cap_techs: a list of technology structures available for cap use
        :param comp_metric: a metric (1=area, 2=loss) used for determining the best component (1=default)
        """
        self.topology = topology
        self.vin = vin
        self.switch_techs = switch_techs
        self.cap_techs = cap_techs
        self.comp_metric = comp_metric
        self._implement()

    def _implement(self):
        # Break out components of topology structure
        ratio = self.topology.ratio
        ac = self.topology.ac
        ar = self.topology.ar
        vc = self.topology.vc * self.vin
        vr = self.topology.vr * self.vin
        vcb = self.topology.vcb * self.vin
        vrb = self.topology.vrb * self.vin

        switch_assign = []
        cap_assign = []
        switch_rel_size = []
        cap_rel_size = []

        # Assign Capacitors
        for i in range(ac.shape[0]):
            Mc = 0
            Cc = 0  # cap cost
            for j, cap_tech in enumerate(self.cap_techs):
                if vc[i] <= cap_tech.rating:
                    # Cap could work ... let's see if it's good
                    # Use area-limited metric, which is usually applicable
                    C = cap_tech.area
                    M = cap_tech.capacitance * vc[i] ** 2 / C
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
                cap_rel_size.append((ac[i] * vc[i]) / (np.sqrt(Mc) * cap_assign[i].area))

        # Assign Switches
        for i in range(ar.shape[0]):

            Msw = 0
            Csw = 0  # switch cost;
            for j, switch_tech in enumerate(self.switch_techs):
                if (vr[i] <= switch_tech.drain_rating):
                    # Switch could work let's see if it's good
                    if self.comp_metric == 2:  # loss metric
                        # assume full gate drive
                        C = switch_tech.gate_cap * switch_tech.gate_rating ** 2 + switch_tech.drain_cap * vr[
                            i] ** 2 + \
                            switch_tech.body_cap * vrb[i] ** 2
                        M = switch_tech.conductance * vr[i] ** 2 / C
                    else:  # area metric
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
                if self.comp_metric == 2:
                    switch_rel_size.append(ar[i] * vr[i] / (np.sqrt(Msw) * switch_assign[i].conductance))
                else:
                    switch_rel_size.append(ar[i] * vr[i] / (np.sqrt(Msw) * switch_assign[i].area))

        # Scale Caps for unit area:

        cap_area = 0
        for i in range(ac.shape[0]):
            cap_area = cap_area + cap_rel_size[i] * cap_assign[i].area

        cap_size = cap_rel_size / (cap_area + 1e-30)

        # Scale Switches for unit area:
        sw_area = 0
        # print(switch_rel_size)
        for i in range(ar.shape[0]):
            #   print(i)
            sw_area = sw_area + switch_rel_size[i] * switch_assign[i].area

        # if sw_area > 0:
        # print(sw_area)
        switch_size = switch_rel_size / sw_area

        self.capacitors = cap_assign
        self.switches = switch_assign
        self.cap_size = cap_size
        self.switch_size = switch_size


    def evaluate_loss(self, vout, iout, fsw, asw, ac):
        """
        Evaluate_loss: evaluates the loss and other peformance metrics for a
        specific size and operating condition of a implemented SC converter

        :param vout: converter output voltage for this calc [V]
        :param iout: converter output current for this calc [A]
        :param fsw: switching frequency [Hz]
        :param asw: switch area [m^2]
        :param ac: capacitor area [m^2]
        :return:
        """
        pass

    def optimize_loss(self, iout, ac):
        """
        Finds the optimal design point for given conditions
        :param iout: converter output current for this calc [A]
        :param ac: capacitor area [m^2]
        :return:
        """
        pass


    def plot_opt_contour(self):
        pass

if __name__ == "__main__":
    from pyseeman.techlib import ITRS16cap, ITRS16sw
    my_topo = Topology("series-parallel", 1, 3)
    my_topo = Topology("ladder", 2, 3)
    my_topo = Topology("dickson", 1, 3)
    my_topo = Topology("cockcroft-walton", 1, 3)
    my_topo = Topology("doubler", 1, 2)
    #my_topo = Topology("fibonacci", 5, 1)

    print(my_topo.__dict__)
    my_imp = my_topo.implement(vin=1,  switch_techs=[ITRS16sw], cap_techs=[ITRS16cap], comp_metric=1)
    print(my_imp.__dict__)
    #my_imp.evaluate_loss(vout=0.6, iout=1, fsw=1e6, asw=1, ac=10)