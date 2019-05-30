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
            for j in range(1, n - m + 1):
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

    elif topology_name.lower() == 'ladder':
        n = num
        m = den
        N = n

        ac = []
        vc = []
        vcb = []

        for j in range(1, n-m):
            ac.extend([j, j])
            vc.extend([1/m, 1/m])
            vcb.append([1/m, 0])

        ac.append(n-m)
        vc.append(1/m)
        vcb.append(1/m)

        for j in range(m-1,0,-1):
            ac.extend([np.ones(2)*j*(n/m-1)])
            vc.extend([1/m, 1/m])
            vcb.extend([1/m, 0])

        ar = np.hstack([np.ones(2*(n-m)), np.ones(2*m)*(n/m-1)])
        vr = np.ones(2*n)/m
        vrb = np.mod(np.linspace(0,(2*n-1),2*n),2)/m
    elif topology_name.lower() == "Dickson":
        raise ValueError('SWITCHCAP:nonIntegerRatio \nthe Dickson topology supports integer ratios only')

        if den != 1:
            raise ValueError("")

        N = num

        # SSL values
        ac = np.ones(N-1)
        vc = []
        vcb = np.ones(N-1)

        for j in range(1,N):
            vc.append(j)

        if N == 2:
            vr = np.ones(4)
            ar = np.ones(4)
            vrb = np.array([0,1,0,1])
        else:
            vr = np.hstack([np.ones(5),2*np.ones(N-2),2])
            ar = np.hstack([np.ones(2)*np.floor((j+1)/2.), np.ones(2)*np.floor(j/2.), np.ones(N)])
            vrb = np.hstack([np.array([0,1,0,1]),np.ones(N)])

    else:
        raise ValueError("Topology type not implemented yet")

    ac = np.array(ac)
    vc = np.array(vc)
    vcb = np.array(vcb)
    ar = np.array(ar)
    vr = np.array(vr)
    vrb = np.array(vrb)

    # TODO: Check if it makes sense that M values are claulcated before the flipping
    ratio = num/den
    Mssl = 2*ratio**2/np.sum(ac*vc)**2
    Mfsl = ratio**2/(2*np.sum(ar*vr)**2)

    if flip == 1:
        ac = ac/ratio
        vc = vc/ratio
        vcb = vcb/ratio
        ar = ar/ratio
        vr = vr/ratio
        vrb = vrb/ratio
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
    res1 = generate_topology("series-parallel",3,1)
    print(res1)

    res2 = generate_topology("ladder",2,3)

    print(res2)

    res3 = generate_topology("Dickson",3,1)
    print(res2)
        #return resulet

