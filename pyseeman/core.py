import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import copy


def fibfun(n):
    if n < 0:
        raise ValueError("Invalid argument for fibfun")
    elif n == 1:
        return 0
    elif n == 2:
        return 1
    else:
        return fibfun(n - 1) + fibfun(n - 2)


class Topology:

    _valid_names = {
        "ladder",
        "dickson",
        "cockcroft-walton",
        "doubler",
        "series-parallel",
        "fibonacci",
    }

    def __init__(self, name, num, den):
        if name.lower() not in self._valid_names:
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
        if self.name.lower() == "series-parallel":
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

        elif self.name.lower() == "ladder":
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
                raise ValueError(
                    "SWITCHCAP:nonIntegerRatio the Dickson topology supports integer ratios only"
                )

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
                ar = np.hstack(
                    [
                        np.ones(2) * np.floor((j + 1) / 2.0),
                        np.ones(2) * np.floor(j / 2.0),
                        np.ones(N),
                    ]
                )
                vrb = np.hstack([np.array([0, 1, 0, 1]), np.ones(N)])

        elif self.name.lower() == "cockcroft-walton":

            if den != 1:
                raise ValueError(
                    "SWITCHCAP:nonIntegerRatio the Cockcroft-Walton topology supports integer ratios only"
                )

            N = num
            # SSL values
            ac = np.array([])
            vc = np.hstack([1, np.ones(N - 2) * 2])
            vcb = np.ones(N - 1)

            for j in range(1, N):
                ac = np.append(np.floor((j + 1) / 2.0), ac)

            # FSL values
            if N == 2:
                vr = np.ones(4)
                ar = np.ones(4)
                vrb = np.array([0, 1, 0, 1])
            else:
                vr = np.hstack([np.ones(5), np.ones(N - 2) * 2, 1])
                ar = np.hstack(
                    [
                        np.floor((j + 1) / 2.0) * np.ones(2),
                        np.floor(j / 2.0) * np.ones(2),
                        np.ones(N),
                    ]
                )
                vrb = np.hstack([np.array([0, 1, 0, 1]), np.ones(N)])

        elif self.name.lower() == "doubler":

            if den != 1:
                raise ValueError(
                    "SWITCHCAP:nonIntegerRatio the Doubler topology supports integer ratios only"
                )

            n = np.ceil(np.log2(num)).astype(np.int)
            N = 2 ** n
            if N != num:
                raise ValueError(
                    "SWITCHCAP:badRatio the Doubler topology supports conversion ratios ~ 2^n"
                )

            # SSL values
            ac = np.array([])
            vc = np.array([])
            vcb = np.array([])

            for j in range(1, 2 * n):
                ac = np.append(2 ** np.floor((j - 1) / 2.0), ac)
                vc = np.append(vc, 2 ** np.floor(j / 2.0))
                vcb = np.append(vcb, 2 ** np.floor(j / 2.0) * np.mod(j, 2))

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
                raise ValueError(
                    "SWITCHCAP:nonIntegerRatio the Fibonacci topology supports integer ratios only"
                )

            i = 2

            while fibfun(i) < num:
                i = i + 1

            if fibfun(i) > num:
                raise ValueError(
                    "SWITCHCAP:badRatio the fibonacci topology supports ratios of F_n or 1/F_n only"
                )
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
        # return f"<Topology(name={self.name}, ratio={self.ratio})>"
        my_string = ""
        for k, v in self.__dict__.items():
            my_string += f"{k} = {v}\n"
        return my_string


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
                raise ValueError(
                    "No capacitors meet the voltage requirement of: {}".format(vc[i])
                )

            # determine relative device size
            if ac[i] == 0:
                cap_rel_size = cap_rel_size.append(0)  # avoid divide by 0
            else:
                cap_rel_size.append(
                    (ac[i] * vc[i]) / (np.sqrt(Mc) * cap_assign[i].area)
                )

        # Assign Switches
        for i in range(ar.shape[0]):

            Msw = 0
            Csw = 0  # switch cost;
            for j, switch_tech in enumerate(self.switch_techs):
                if vr[i] <= switch_tech.drain_rating:
                    # Switch could work let's see if it's good
                    if self.comp_metric == 2:  # loss metric
                        # assume full gate drive
                        C = (
                            switch_tech.gate_cap * switch_tech.gate_rating ** 2
                            + switch_tech.drain_cap * vr[i] ** 2
                            + switch_tech.body_cap * vrb[i] ** 2
                        )
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
                raise ValueError(
                    "No switches meet the voltage requirement of: {}".format(vr[i])
                )

            # determine relative device size
            if ar[i] == 0:
                switch_rel_size.append(0)
            else:
                if self.comp_metric == 2:
                    switch_rel_size.append(
                        ar[i] * vr[i] / (np.sqrt(Msw) * switch_assign[i].conductance)
                    )
                else:
                    switch_rel_size.append(
                        ar[i] * vr[i] / (np.sqrt(Msw) * switch_assign[i].area)
                    )

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

    def _expand_input(self, input, maxsize):
        # TODO: Finish this for cases with more dimensions
        if input.shape == (1, 1):
            # Scalar input
            return input * np.ones(maxsize)

        elif input.shape == (1, maxsize[1]):
            # Row vector input
            return input * np.ones((maxsize[0], 1))

        elif input.shape == (maxsize[0], 1):
            # Column vector input
            return input * np.ones((1, maxsize[1]))
        elif input.shape == (maxsize[0], maxsize[1]):
            # Input already a properly-sized array
            return input
        elif input.size == 0:
            # Input is empty
            return input
        elif input.shape == (0,):
            raise ValueError("Only fsw or Vout can be empty")
        else:
            raise ValueError(
                "All inputs must have the same number of rows and columns (if not 1)"
            )

    def evaluate_loss(self, vin, vout, iout, fsw, area_sw, area_cap):
        """
        Evaluate_loss: evaluates the loss and other peformance metrics for a
        specific size and operating condition of a implemented SC converter

        :param vout: converter output voltage for this calc [V]
        :param iout: converter output current for this calc [A]
        :param fsw: switching frequency [Hz]
        :param area_sw: switch area [m^2]
        :param ac: capacitor area [m^2]
        :return:
        """

        # Converter to numpy arrays
        vin = np.asarray(vin)
        vout = np.asarray(vout)
        iout = np.asarray(iout)
        fsw = np.asarray(fsw)
        area_sw = np.asarray(area_sw)
        area_cap = np.asarray(area_cap)

        # If scalar or 1D array convert to 2D row vector
        if vin.ndim == 0:
            vin = vin.reshape((1, 1))
        elif vin.ndim == 1:
            vin = vin.reshape((1, vin.shape[0]))

        if vout.ndim == 0:
            vout = vout.reshape((1, 1))
        elif vout.ndim == 1:
            vout = vout.reshape((1, vout.shape[0]))

        if iout.ndim == 0:
            iout = iout.reshape((1, 1))
        elif iout.ndim == 1:
            iout = iout.reshape((1, iout.shape[0]))

        if fsw.ndim == 0:
            fsw = fsw.reshape((1, 1))
        elif fsw.ndim == 1:
            fsw = fsw.reshape((1, fsw.shape[0]))

        if area_sw.ndim == 0:
            area_sw = area_sw.reshape((1, 1))
        elif area_sw.ndim == 1:
            area_sw = area_sw.reshape((1, area_sw.shape[0]))

        if area_cap.ndim == 0:
            area_cap = area_cap.reshape((1, 1))
        elif area_cap.ndim == 1:
            area_cap = area_cap.reshape((1, area_cap.shape[0]))

        # Break out components of topology structure
        ratio = self.topology.ratio
        ac = self.topology.ac
        ar = self.topology.ar
        vc = self.topology.vc  # * self.vin
        vr = self.topology.vr  # * self.vin
        vcb = self.topology.vcb  # * self.vin
        vrb = self.topology.vrb  # * self.vin

        caps = self.capacitors
        cap_size = self.cap_size
        switches = self.switches
        sw_size = self.switch_size

        # Expand input parameters:
        # If an input is given as a vector then expand all other inputs
        # to vectors.
        # If two inputs are given as a row and column vector, respectively,
        # then expand inputs to 2D arrays.

        eval_type = 0  # 0: undefined, 1: vout, 2: fsw

        paramdim = np.max(
            [
                vin.shape,
                vout.shape,
                iout.shape,
                fsw.shape,
                area_sw.shape,
                area_cap.shape,
            ],
            axis=0,
        )
        vin = self._expand_input(vin, paramdim)

        # If Vout is empty, then evaluate to find Vout
        if vout.size == 0:
            eval_type = 1
            vout = np.zeros(paramdim)
        else:
            vout = self._expand_input(vout, paramdim)

        iout = self._expand_input(iout, paramdim)

        # If fsw is empty then evaluate to find fsw
        if fsw.size == 0:
            if eval_type == 1:
                raise ValueError("Both fsw and Vout cannot be undefined")
            else:
                eval_type = 2
                fsw = np.zeros(paramdim)
        else:
            fsw = self._expand_input(fsw, paramdim)

        area_sw = self._expand_input(area_sw, paramdim)
        area_cap = self._expand_input(area_cap, paramdim)

        # ********************** Start analysis ***********************#
        # Calculate SSL output resistance
        Rssl_alpha = 0

        for i in range(len(caps)):
            if ac[i] > 0:
                Rssl_alpha += ac[i] ** 2 / (caps[i].capacitance * cap_size[i])

        # Calculate FSL output resistance
        Rfsl_alpha = 0

        for i in range(len(switches)):
            if ar[i] > 0:
                Rfsl_alpha += 2 * ar[i] ** 2 / (switches[i].conductance * sw_size[i])
        Rfsl = Rfsl_alpha / area_sw

        # Calculate ESR loss
        Resr_alpha = 0

        for i in range(len(caps)):
            if ac[i] > 0:
                Resr_alpha += 4 * ac[i] ** 2 * caps[i].esr / cap_size[i]
        Resr = Resr_alpha / area_cap
        if hasattr(self, "esr"):
            Resr += self.esr

        # Calculate the unknown variable
        if eval_type == 1:
            # Vout is unknown
            Rssl = Rssl_alpha / (fsw * area_cap)

            # Calculate total output resistance
            Rout = np.sqrt(Rssl ** 2 + (Rfsl + Resr) ** 2)
            vout = vin * ratio - Rout * iout
            Pout = vout * iout
            is_prac = np.ones(paramdim)

        elif eval_type == 2:
            # fsw is unknown
            #           # Calculate needed output resistance and switching frequency to match
            # output voltage
            # is_prac is 1 if a finite fsw that satisfies Iout, Vin, Vout exists

            Rreq = (vin * ratio - vout) / iout
            is_prac = ((Rreq > 0) & (Rfsl + Resr < Rreq)) * 1.0
            Rssl = np.real(np.sqrt(Rreq ** 2 - (Rfsl + Resr) ** 2))
            fsw = Rssl_alpha / (Rssl * area_cap)

            # Calculate total output resistance
            Rout = np.sqrt(Rssl ** 2 + (Rfsl + Resr) ** 2)
            Pout = vout * iout
        else:
            raise ValueError("Either Vout or fsw must be []")

        # Calculate resistance losses
        Pssl = Rssl * iout ** 2
        Pfsl = Rfsl * iout ** 2
        Pesr = Resr * iout ** 2
        Pres = Rout * iout ** 2

        # Calculate cap related parasitic loss
        Pc_alpha = 0
        for i in range(len(caps)):
            Pc_alpha += caps[i].bottom_cap * cap_size[i] * vcb[i] ** 2

        Pc = Pc_alpha * fsw * area_cap * vin ** 2

        # Calculate switch related parasitic loss
        Psw_alpha = 0
        Pg_alpha = 0
        for i in range(len(switches)):
            # Assume switch is driven at full gate_rating voltage
            Vgssw = switches[i].gate_rating
            Pg_alpha += switches[i].gate_cap * sw_size[i] * Vgssw ** 2
            Psw_alpha += (
                switches[i].drain_cap * sw_size[i] * vr[i] ** 2
                + switches[i].body_cap * sw_size[i] * vrb[i] ** 2
            )
        Psw = (Psw_alpha * vin ** 2 + Pg_alpha) * fsw * area_sw

        # Calculate total loss, efficiency, etc.
        Ploss = Pres + Pc + Psw
        eff = Pout / (Pout + Ploss)

        # Find dominant loss
        loss_arr = np.dstack((Pssl, Pfsl, Pesr, Pc, Psw))
        dominant_loss_args = np.argmax(loss_arr, axis=2)
        texts = ["SSL Loss", "FSL Loss", "ESR Loss", "Bottom-plate", "Switch Parasitic"]

        # Pack performance parameters
        # TODO: Decide whether we should use something else than a dict here

        performance = {}
        performance["Vout"] = vout
        performance["fsw"] = fsw
        performance["is_possible"] = is_prac
        performance["efficiency"] = eff
        performance["total_loss"] = Ploss
        performance["impedance"] = Rout
        performance["dominant_loss"] = dominant_loss_args
        performance["dominant_text"] = [
            [texts[x] for x in row] for row in dominant_loss_args
        ]

        return performance

    def optimize_loss(self, iout, area_cap):
        """
        Finds the optimal design point for given conditions
        :param iout: converter output current for this calc [A]
        :param area_cap: capacitor area [m^2]
        :return:
        """
        opt_func = lambda x: self.evaluate_loss(
            vin=self.vin,
            vout=[],
            iout=iout,
            fsw=np.exp(x[0]),
            area_sw=np.exp(x[1]),
            area_cap=area_cap,
        )["total_loss"]
        x0 = np.array([10, -10])
        result = minimize(opt_func, x0, bounds=((1, 100), (-100, 1)), tol=1e-12)
        if result["success"] is not True:
            raise RuntimeError("optimization is not sucessfull")
        performance = self.evaluate_loss(
            self.vin, [], iout, np.exp(result.x[0]), np.exp(result.x[1]), area_cap
        )
        return performance, np.exp(result.x[0]), np.exp(result.x[1])

    def __repr__(self):
        # return f"<Topology(name={self.name}, ratio={self.ratio})>"
        my_string = ""
        for k, v in self.__dict__.items():
            my_string += f"{k} = {v}\n"
        return my_string


def plot_opt_contour(imp, vin, iout, area_cap, plot_points=100, plot_axes=None):
    """
    Create efficiency contour plot of efficiency
    :param imp: implementation object
    :param vin: input voltage [V]
    :param iout: output current [A]
    :param area_cap: Capacitor area [mm^2]
    :param plot_points: Number of points to plot (default: 100)
    :param plot_axes: Array of plot limits 10-exponents (optional)
    :return:
    """
    # Find optimal operating point
    [opt_perf, fsw_opt, asw_opt] = imp.optimize_loss(iout, area_cap)

    # Define plot region around the optimum
    if plot_axes is None:
        fsw_min = np.floor(np.log10(fsw_opt) - 1)
        fsw_max = np.ceil(np.log10(fsw_opt) + 1)
        asw_min = np.floor(np.log10(asw_opt) - 1)
        asw_max = np.ceil(np.log10(asw_opt) + 1)

    # Generate plot mesh and evaluate performance
    fsw, asw = np.meshgrid(
        np.logspace(fsw_min, fsw_max, plot_points),
        np.logspace(asw_min, asw_max, plot_points),
    )
    p = imp.evaluate_loss(vin, [], iout, fsw, asw, area_cap)

    # Find indices of optimum point
    idx_max = np.where(p["efficiency"] == p["efficiency"].max())

    # Plot contours, maximum point, and dominant loss regions
    fig, ax = plt.subplots()
    ax.contour(fsw, asw * 1e6, p["dominant_loss"], [0.5, 1.5, 2.5, 3.5], colors="grey")
    cs_eff = ax.contour(
        fsw,
        asw * 1e6,
        p["efficiency"],
        [0.05, 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1],
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.loglog(
        fsw[idx_max], asw[idx_max] * 1e6, marker="*", color="black", markersize=20
    )
    ax.clabel(cs_eff, inline=True, fontsize=10, inline_spacing=1)
    ax.set_xlabel("Switching frequency [Hz]")
    ax.set_ylabel("Switch area [mm^2]")
    ax.set_title(
        "Iout = %.2e A, Ac = %.4f mm^2, Eff = %.1f%%"
        % (iout, area_cap * 1e6, p["efficiency"].max() * 100)
    )
    return ax


def plot_regulation(
    topologies, vin, vout, iout, area_cap, switches, capacitors, esr=0, idesign=0
):
    """

    :param topologies: A matrix of topologies and ratios
    :param vin: Input voltage of converter (could be a vector)
    :param vout: Output voltage of converter (a vector if Vin is a scalar)
    :param iout: Matching vector of output currents [A]
    :param area_cap: Capacitor area constraint (in m^2).  fsw will be swept, Asw will be chosen automatically
    :param switches: a row vector of switch technology structures
    :param capacitors: a row vector of capacitor technology structures
    :param esr: the output-referred constant esr of requisite metal (ie, bondwires).  Default = 0
    :param idesign: a vector (size = number of topologies) containing the nominal design current for each topology
    :return:
    """
    mode = 0  # 1 = Vout, 2 = Vin swept
    if isinstance(vout, (int, float)):
        mode = 2
        xlabel = "Input Voltage [V]"
        xlim = (vin.min(), vin.max())
        title = f"Regulation @Vout={vout}"

    if isinstance(vin, (int, float)):
        if mode == 2:
            raise RuntimeError("Cannot sweep both Vin and Vout")
        mode = 1
        vin_nom = vin
        xlabel = "Output Voltage [V]"
        xlim = (vout.min(), vout.max())
        title = f"Regulation @Vin={vin}"

    numtops = len(topologies)

    eff = [0] * numtops

    fig, ax = plt.subplots()

    for i, t in enumerate(topologies):
        if mode == 2:
            vin_nom = vout / t.ratio
        imp = t.implement(vin_nom, switches, capacitors)
        [opt_perf, fsw_opt, asw_opt] = imp.optimize_loss(iout, area_cap)
        p = imp.evaluate_loss(vin, vout, iout, [], asw_opt, area_cap)
        eff = 100 * p["efficiency"] * p["is_possible"]
        eff[eff == 0] = "nan"  # or use np.nan
        eff_max = np.nanmax(eff)
        x_eff_max = np.nanargmax(eff)
        if mode == 1:
            vout_1 = p["Vout"] * p["is_possible"]

            eff_trace, = ax.plot(
                np.squeeze(vout_1), np.squeeze(eff), label=f"{t.num}:{t.den}"
            )
            ax.plot(
                np.squeeze(vout_1)[x_eff_max],
                eff_max,
                marker="o",
                color=eff_trace.get_color(),
            )
        if mode == 2:
            eff_trace, = ax.plot(
                np.squeeze(vin), np.squeeze(eff), label=f"{t.num}:{t.den}"
            )
            ax.plot(
                np.squeeze(vin)[x_eff_max],
                eff_max,
                marker="o",
                color=eff_trace.get_color(),
            )

    ax.set(xlabel=xlabel, ylabel="Efficiency [%]", title=title, xlim=xlim)
    ax.grid()
    ax.legend()
    return ax


def cascade_topologies(topology1, topology2):
    """
    Returns a new topology consisting of series connection of the two input topologies

    :param topology1: First topology (connected to input)
    :param topology2: Second topology (connected to output)
    :return: New object of the cascaded topology
    """
    casc_top = copy.deepcopy(topology1)

    ratio1 = topology1.ratio
    ratio2 = topology2.ratio

    casc_top.name = topology1.name + "->" + topology2.name
    casc_top.ac = np.hstack([topology1.ac * ratio2, topology2.ac])
    casc_top.ar = np.hstack([topology1.ar * ratio2, topology2.ar])
    casc_top.vc = np.hstack([topology1.vc, topology2.vc * ratio1])
    casc_top.vcb = np.hstack([topology1.vcb, topology2.vcb * ratio1])
    casc_top.vr = np.hstack([topology1.vr, topology2.vr * ratio1])
    casc_top.vrb = np.hstack([topology1.vrb, topology2.vrb * ratio1])
    casc_top.ratio = ratio1 * ratio2
    casc_top.Mssl = 2 * casc_top.ratio ** 2 / (np.sum(casc_top.ac * casc_top.vc)) ** 2
    casc_top.Mfsl = casc_top.ratio ** 2 / (2 * np.sum(casc_top.ar * casc_top.vr)) ** 2

    return casc_top


def permute_topologies(topologies1, topologies2):
    """
    For two lists of topology objects, returns a list of every permutation
    of topologies in the first list connected in series (cascaded) with
    the topologies of the second list,

    :param topologies1: First list of topologies
    :param topologies2: Second list of topologies
    :return: list of topology permutations
    """
    newtops = []

    for m in range(len(topologies1)):
        for n in range(len(topologies2)):
            top1 = topologies1[m]
            top2 = topologies2[n]

            newtops.append(cascade_topologies(top1, top2))

    return newtops


if __name__ == "__main__":
    plt.close("all")
    from pyseeman.techlib import ITRS16cap, HVcap, ITRS16sw

    top_stepdown_1 = Topology("Ladder", 2, 3)
    top_stepdown_2 = Topology("Dickson", 1, 3)
    top_stepdown_3 = Topology("Cockcroft-Walton", 1, 4)
    top_stepdown_4 = Topology("Doubler", 1, 2)

    top_stepup_1 = Topology("Ladder", 3, 2)
    top_stepup_2 = Topology("Dickson", 3, 1)
    top_stepup_3 = Topology("Cockcroft-Walton", 4, 1)
    top_stepup_4 = Topology("Doubler", 2, 1)

    imp_stepdown_1 = top_stepdown_1.implement(1.5, [ITRS16sw], [ITRS16cap, HVcap], 1)
    imp_stepdown_2 = top_stepdown_2.implement(1.5, [ITRS16sw], [ITRS16cap, HVcap], 1)
    imp_stepdown_3 = top_stepdown_3.implement(1.5, [ITRS16sw], [ITRS16cap, HVcap], 1)
    imp_stepdown_4 = top_stepdown_4.implement(1.5, [ITRS16sw], [ITRS16cap, HVcap], 1)

    imp_stepup_1 = top_stepup_1.implement(0.25, [ITRS16sw], [ITRS16cap, HVcap], 1)
    imp_stepup_2 = top_stepup_2.implement(0.25, [ITRS16sw], [ITRS16cap, HVcap], 1)
    imp_stepup_3 = top_stepup_3.implement(0.25, [ITRS16sw], [ITRS16cap, HVcap], 1)
    imp_stepup_4 = top_stepup_4.implement(0.25, [ITRS16sw], [ITRS16cap, HVcap], 1)

    ax1 = plot_opt_contour(imp_stepdown_1, 2.0, 10e-3, 1e-6)
    ax2 = plot_opt_contour(imp_stepdown_2, 2.0, 10e-3, 1e-6)
    ax3 = plot_opt_contour(imp_stepdown_3, 2.0, 10e-3, 1e-6)
    ax4 = plot_opt_contour(imp_stepdown_4, 2.0, 10e-3, 1e-6)

    ax5 = plot_opt_contour(imp_stepup_1, 0.25, 10e-3, 1e-6)
    ax6 = plot_opt_contour(imp_stepup_2, 0.25, 10e-3, 1e-6)
    ax7 = plot_opt_contour(imp_stepup_3, 0.25, 10e-3, 1e-6)
    ax8 = plot_opt_contour(imp_stepup_4, 0.25, 10e-3, 1e-6)

    # test optimization
    print(imp_stepdown_1.optimize_loss(iout=1e-3, area_cap=1e-6))

    topologies = [
        Topology("series-parallel", 1, 3),
        Topology("series-parallel", 1, 2),
        Topology("series-parallel", 2, 3),
        Topology("series-parallel", 4, 5),
    ]

    # Test Vout sweep
    Vout = np.linspace(0, 1.2, 1000).reshape((1, 1000))
    plot_regulation(topologies, 1.2, Vout, 500e-3, 1e-7, [ITRS16sw], [ITRS16cap])

    # Test Vin sweep
    ITRS16cap.rating = 3
    ITRS16sw.drain_rating = 3
    ITRS16sw.gate_rating = 3

    Vin = np.linspace(1, 5, 1000).reshape((1, 1000))
    plot_regulation(topologies, Vin, 1.2, 100e-3, 1e-6, [ITRS16sw], [ITRS16cap])

    plt.show()
