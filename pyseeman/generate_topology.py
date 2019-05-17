import numpy as np

class Topology:
    #def __init__(self, ac=None):
    #    self.ac = ac
    def __init__(self, *args, **kwargs):
        pass

    def __repr__(self):
        my_string = ""
        for k, v in self.__dict__.items():
            my_string += f"{k} = {v}\n"
        return my_string

def generate_topology(*args):
    """
    Generates a topology based on a choice of type of topology along
    with a specification of the desired ideal voltage conversion
    ratio.
    """
    if args is None or len(args) < 2:
        raise ValueError("generate_topology() takes at least two arguments")
    elif(len(args) > 3):
        raise ValueError("generate_topology() takes at most three arguments")
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

    # If input ratio is step-down, flip to step up for ac computation
    # and then flip it again before returning result.
    flip = 0
    if num/den < 1:
        flip = 1
        den, num = num, den
    #*************************** SERIES - PARALLEL *******************************#
    if topology_name.lower() == 'series-parallel':
        n = num
        m = den
        N = n

        # SSL values
        ac = np.ones(m*(n-m)//m)
        vc = np.ones(m*(n-m)//m)
        vcb = []

        for i in range(1, m + 1):
            for j in range(1, n - m + 2):
                vcb.append((i+j-1)/m)

        # FSL values
        vr = []
        vrb = []
        print(f"n = {n}")
        print(f"m = {m}")
        for i in range(1, m + 1):
            for j in range(1, n - m + 2):
                if j == 1:
                    vr.append(i/m)
                    vrb.append((i + j - 1)/m)
                    print("0")
                elif j == n - m + 1:
                    vr.append((n-m-1+i)/m)
                    vrb.append((i+j-2)/m)
                    print("1")
                else:
                    vr.append(1/m)
                    vrb.append((i+j-1)/m)
                    print("2")
        for i in range(1,m+2):
            for j in range(1, n-m+1):
                if i == 1:
                    vr.append(j/m)
                elif i == m+1:
                    vr.append((m-1+j)/m)
                else:
                    vr.append(1/m)
                if i == 1 or i == m+1:
                    vrb.append(0)
                    print("3")
                else:
                    vrb.append((i+j-2)/m)
                    print("4)")
        ar = np.ones(len(vr))/m
        print(vrb)
    else:
        raise ValueError("Topology type not implemented yet")

    ratio = num/den
    Mssl = 2*ratio**2/np.sum(ac*vc)**2
    Mfsl = ratio**2/(2*np.sum(ar*vr)**2)

    if flip == 1:
        ac = np.array(ac)/ratio
        vc = np.array(vc)/ratio
        vcb = np.array(vcb)/ratio
        ar = np.array(ar)/ratio
        vr = np.array(vr)/ratio
        vrb = np.array(vrb)/ratio
        ratio = 1/ratio

    result =  Topology()
    result.topName = topology_name
    result.ac = ac
    result.vc = vc
    result.vcb = vcb
    result.ar = ar
    result.vr = vr
    result.vrb = vrb
    result.Mssl = Mssl
    result.Mfsl = Mfsl
    result.ratio = ratio

    return result

if __name__ == "__main__":
    res1 = generate_topology("series-parallel",1,3)
    print(res1)
        #return result

