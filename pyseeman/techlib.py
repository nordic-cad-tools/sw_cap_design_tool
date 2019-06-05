
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
    """
    45 nm Oxide Capacitors
    """
    tech_name = 'ITRS 45nm'
    dev_name = 'ITRS 45nm Oxide Capacitor'
    capacitance = (40 ** 2) * .84e-15 * .5 / .045
    area = 40e-6 * 40e-6
    bottom_cap = 135.3e-15
    esr = 0
    rating = 1.2

class ITRS32cap:
    """
    32 nm Oxide Capacitors
    """
    tech_name = 'ITRS 32nm'
    dev_name = 'ITRS 32nm Oxide Capacitor'
    capacitance = (40 ** 2) * .8e-15 * .45 / .032
    area = 40e-6 * 40e-6
    bottom_cap = 135.3e-15
    esr = 0
    rating = 1.2

class ITRS22cap:
    """
    22 nm Oxide Capacitors
    """
    tech_name = 'ITRS 22nm'
    dev_name = 'ITRS 22nm Oxide Capacitor'
    capacitance = (40 ** 2) * .58e-15 * .4 / .022
    area = 40e-6 * 40e-6
    bottom_cap = 135.3e-15
    esr = 0
    rating = 1.2


class ITRS16cap:
    """
    16 nm Oxide Capacitors
    """
    tech_name = 'ITRS 16nm'
    dev_name = 'ITRS 16nm Oxide Capacitor'
    capacitance = (40 ** 2) * .48e-15 * .35 / .016
    area = 40e-6 * 40e-6
    bottom_cap = 135.3e-15
    esr = 0
    rating = 1.2

#TODO check the name and type of that cap
class HVcap:
    """
    16 nm Oxide Capacitors
    """
    tech_name = 'ITRS 16nm'
    dev_name = 'ITRS 16nm Oxide Capacitor'
    capacitance = (40 ** 2) * .48e-15 * .35 / .016 * 0.01
    area = 40e-6 * 40e-6
    bottom_cap = 135.3e-15
    esr = 0
    rating = 1.2

class ITRS90sw:
    """
    90 nm switch
    """
    tech_name = 'ITRS 90 nm'
    dev_name = 'ITRS 90 nm NMOS native transistor'
    area = 1e-6 * 90e-9 / .4
    conductance = 1.11e-3
    gate_rating = 1.25
    drain_rating = 1.25
    gate_cap = 1e-15
    drain_cap = .33 * gate_cap
    body_cap = .2 * gate_cap

class ITRS65sw:
    """
    65 nm switch
    """
    tech_name = 'ITRS 65 nm'
    dev_name = 'ITRS 65 nm NMOS native transistor'
    area = 1e-6 * 65e-9 / .35
    conductance = 1.3e-3
    gate_rating = 1.25
    drain_rating = 1.25
    gate_cap = 0.9e-15
    drain_cap = .33 * gate_cap
    body_cap = .2 * gate_cap

class ITRS45sw:
    """
    45 nm switch
    """
    tech_name = 'ITRS 45 nm'
    dev_name = 'ITRS 45 nm NMOS native transistor'
    area = 1e-6 * 45e-9 / .35
    conductance = 1.51e-3
    gate_rating = 1.25
    drain_rating = 1.25
    gate_cap = 0.84e-15
    drain_cap = .33 * gate_cap
    body_cap = .2 * gate_cap

class ITRS32sw:
    """
    32 nm switch
    """
    tech_name = 'ITRS 32 nm'
    dev_name = 'ITRS 32 nm NMOS native transistor'
    area = 1e-6 * 32e-9 / .35
    conductance = 1.82e-3
    gate_rating = 1.25
    drain_rating = 1.25
    gate_cap = 0.8e-15
    drain_cap = .33 * gate_cap
    body_cap = .2 * gate_cap

class ITRS22sw:
    """
    22 nm switch
    """
    tech_name = 'ITRS 22 nm'
    dev_name = 'ITRS 22 nm NMOS native transistor'
    area = 1e-6 * 22e-9 / .3
    conductance = 2.245e-3
    gate_rating = 1.25
    drain_rating = 1.25
    gate_cap = 0.58e-15
    drain_cap = .33 * gate_cap
    body_cap = .2 * gate_cap


class ITRS16sw:
    """
    16 nm switch
    """
    tech_name = 'ITRS 16 nm'
    dev_name = 'ITRS 16 nm NMOS native transistor'
    area = 1e-6 * 16e-9 / .25
    conductance = 2.535e-3
    gate_rating = 1.25
    drain_rating = 1.25
    gate_cap = 0.48e-15
    drain_cap = .33 * gate_cap
    body_cap = .2 * gate_cap
