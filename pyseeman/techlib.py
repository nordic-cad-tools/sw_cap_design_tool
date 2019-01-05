
class ITRS90cap:
    """
    90 nm Oxide Capacitors
    """
    tech_name = 'ITRS 90nm'
    dev_name = 'ITRS 90nm Oxide Capacitor'
    capacitance = (40 ** 2) * 1e-15 * .7 / .09
    area = 40e-6 * 40e-6
    bottom_cap = 135.3e-15
    esr = 0
    rating = 1.2

class ITRS65cap:
    """
    65 nm Oxide Capacitors
    """
    tech_name = 'ITRS 65nm'
    dev_name = 'ITRS 65nm Oxide Capacitor'
    capacitance = (40 ** 2) * 0.9e-15 * .6 / .065
    area = 40e-6 * 40e-6
    bottom_cap = 135.3e-15
    esr = 0
    rating = 1.2

class ITRS45cap:
    pass

class ITRS32cap:
    pass

class ITRS22cap:
    pass

class ITRS16cap:
    pass

class HVcap:
    pass

class ITRS90sw:
    pass

class ITRS65sw:
    pass

class ITRS45sw:
    pass

class ITRS32sw:
    pass

class ITRS22sw:
    pass

class ITRS16sw:
    pass