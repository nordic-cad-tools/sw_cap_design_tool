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
        self.ac = 1
        self.vc = 2
        self.vcb = 3
        self.ar = 4
        self.vr = 5
        self.vrb = 6
        self.Mssl = 7
        self.Mfsl = 8
        self.ratio = 9

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

    def plot_opt_contour(self):
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

    def evaluate_loss(self, vin, vout, iout, fsw, asw, ac):
        """
        Evaluate_loss: evaluates the loss and other peformance metrics for a
        specific size and operating condition of a implemented SC converter

        :param vin: converter input voltage for this calc [V]
        :param vout: converter output voltage for this calc [V]
        :param iout: converter output current for this calc [A]
        :param fsw: switching frequency [Hz]
        :param asw: switch area [m^2]
        :param ac: capacitor area [m^2]
        :return:
        """
        pass

    def optimize_loss(self, vin, iout, ac):
        """
        Finds the optimal design point for given conditions
        :param vin: converter input voltage for this calc [V]
        :param iout: converter output current for this calc [A]
        :param ac: capacitor area [m^2]
        :return:
        """
        pass

if __name__ == "__main__":
    topo = Topology("series-parallel", 1, 3)
    print(topo)