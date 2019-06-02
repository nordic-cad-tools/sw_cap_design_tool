import numpy as np


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

    def _implement(self):
        switch_assign = []
        cap_assign = []
        switch_rel_size = []
        cap_rel_size = []

        self.capacitors = cap_assign
        self.switches = switch_assign
        self.cap_size = cap_size = None
        self.switch_size = switch_size = None

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
    print(my_topo.__dict__)
    my_imp = my_topo.implement(vin=2,  switch_techs=[ITRS16sw], cap_techs=[ITRS16cap], comp_metric=1)
    my_imp.evaluate_loss(vout=0.6, iout=1, fsw=1e6, asw=1, ac=10)