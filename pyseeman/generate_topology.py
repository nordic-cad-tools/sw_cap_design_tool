import numpy as np


def generate_topology(*args):
    """
    Generates a topology based on a choice of type of topology along
    with a specification of the desired ideal voltage conversion
    ratio.
    """
    if args is None or len(args) < 2:
        raise ValueError("generate_topology() takes at least two arguments")
    else:
        topology_name = args[0]

    # If only two arguments are given, convert float
    # num input to num and den
    if len(args) == 2:
        frac_val = args[1].as_interger_ratio()
        num = frac_val[0]
        den = frac_val[1]
    else:
        num = args[1]
        den = args[2]

    flip = 0
    if num/den < 1:
        # t = den
        # den = num
        # num = t
        # flip = 1
        den, num = num, den

    if topology_name.lower() is 'series-parallel':
        n = num
        m = den
        N = n

        ac = np.ones(m*(n-m)//m)
        vc = np.ones(m*(n-m)//m)
        vcb = []

        for i in range(1, m):
            for j in range(1, n-m+1):
                vcb.append((i+j)/m)

        vr = []
        vrb = []

        for i in range(1, m):
            for j in range(1, n-m+1):
                vr.append(i/m)
                vrb.append((i+j-1)/m)


    return 0
    #return result

