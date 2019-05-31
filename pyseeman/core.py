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

    def __repr__(self):
        return f"<Topology(name={self.name}, ratio={self.ratio})>"


if __name__ == "__main__":
    topo = Topology("series-parallel", 1, 3)
    print(topo)