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


def fibfun(n):
    if n<0:
        raise ValueError("Invalid argument for fibfun")
    elif n == 1:
        return 0
    elif n == 2:
        return 1
    else:
        return fibfun(n-1) + fibfun(n-2)

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
        vcb = np.array([])

        for i in range(1, m + 1):
            for j in range(1, n - m + 1):
                vcb = np.append(vcb,(i+j-1)/m)

        # FSL values
        vr = np.array([])
        vrb = np.array([])
        for i in range(1, m + 1):
            for j in range(1, n - m + 2):
                if j == 1:
                    vr = np.append(vr, i/m)
                    vrb = np.append(vrb, (i + j - 1)/m)
                elif j == n - m + 1:
                    vr = np.append(vr, (n-m-1+i)/m)
                    vrb = np.append(vrb, (i+j-2)/m)
                else:
                    vr = np.append(vr, 1/m)
                    vrb = np.append(vrb, (i+j-1)/m)
        for i in range(1,m+2):
            for j in range(1, n-m+1):
                if i == 1:
                    vr = np.append(vr, j/m)
                elif i == m+1:
                    vr = np.append(vr, (m-1+j)/m)
                else:
                    vr = np.append(vr, 1/m)
                if i == 1 or i == m+1:
                    vrb = np.append(vrb, 0)
                else:
                    vrb = np.append(vrb, (i+j-2)/m)
        ar = np.ones(len(vr))/m

    elif topology_name.lower() == 'ladder':
        n = num
        m = den
        N = n

        ac = np.array([])
        vc = np.array([])
        vcb = np.array([])

        for j in range(1, n-m):
            ac = np.append(ac, [j, j])
            vc = np.append(vc, [1/m, 1/m])
            vcb = np.append(vcb, [1/m, 0])

        ac = np.append(ac, n-m)
        vc = np.append(vc, 1/m)
        vcb = np.append(vcb, 1/m)

        for j in range(m-1,0,-1):
            ac = np.append(ac, np.ones(2)*j*(n/m-1))
            vc = np.append(vc, [1/m, 1/m])
            vcb = np.append(vcb, [1/m, 0])

        ar = np.hstack([np.ones(2*(n-m)), np.ones(2*m)*(n/m-1)])
        vr = np.ones(2*n)/m
        vrb = np.mod(np.linspace(0,(2*n-1),2*n),2)/m
    elif topology_name.lower() == "dickson":

        if den != 1:
            raise ValueError('SWITCHCAP:nonIntegerRatio the Dickson topology supports integer ratios only')

        N = num

        # SSL values
        ac = np.ones(N-1)
        vc = np.array([])
        vcb = np.ones(N-1)

        for j in range(1,N):
            vc = np.append(vc, j)

        if N == 2:
            vr = np.ones(4)
            ar = np.ones(4)
            vrb = np.array([0,1,0,1])
        else:
            vr = np.hstack([np.ones(5),2*np.ones(N-2),2])
            ar = np.hstack([np.ones(2)*np.floor((j+1)/2.), np.ones(2)*np.floor(j/2.), np.ones(N)])
            vrb = np.hstack([np.array([0,1,0,1]),np.ones(N)])

    elif topology_name.lower() == "cockcroft-walton":

        if den != 1:
            raise ValueError('SWITCHCAP:nonIntegerRatio the Cockcroft-Walton topology supports integer ratios only')

        N = num
        # SSL values
        ac = np.array([])
        vc = np.hstack([1, np.ones(N-2)*2])
        vcb = np.ones(N-1)

        for j in range(1,N):
            ac = np.append(np.floor((j+1)/2.), ac)

        # FSL values
        if N == 2:
            vr = np.ones(4)
            ar = np.ones(4)
            vrb = np.array([0, 1, 0, 1])
        else:
            vr = np.hstack([np.ones(5), np.ones(N-2)*2, 1])
            ar = np.hstack([np.floor((j+1)/2.)*np.ones(2), np.floor(j/2.)*np.ones(2), np.ones(N)])
            vrb = np.hstack([np.array([0, 1, 0, 1]), np.ones(N)])

    elif topology_name.lower() == "doubler":

        if den != 1:
            raise ValueError('SWITCHCAP:nonIntegerRatio the Doubler topology supports integer ratios only')

        n = np.ceil(np.log2(num)).astype(np.int)
        N = 2**n
        if N != num:
            raise ValueError('SWITCHCAP:badRatio the Doubler topology supports conversion ratios ~ 2^n')

        # SSL values
        ac = np.array([])
        vc = np.array([])
        vcb = np.array([])

        for j in range(1,2*n):
            ac = np.append(2**np.floor((j-1)/2.), ac)
            vc = np.append(vc, 2**np.floor(j/2.))
            vcb = np.append(vcb, 2**np.floor(j/2.)*np.mod(j,2))

        # FSL values
        ar = np.array([])
        vr = np.array([])
        vrb = np.array([])

        for j in range(1,n+1) :
            ar = np.append(ar, np.ones(4)*2**(j-1))
            vr = np.append(vr, np.ones(4)*2**(n-j))
            vrb = np.append(vrb, np.array([0, 1, 0, 1])*2**(n-j))

    elif topology_name.lower() == "fibonacci":
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


        for j in range(2,i):
            ac = np.append(fibfun(j-1), ac)
            vc = np.append(vc, fibfun(j-1))
            vcb = np.append(vcb, fibfun(j-1))

        ar = np.array([1])
        vr = np.array([])
        vrb = np.array([0])

        for j in range(2,i):
            ar = np.append(np.array([fibfun(j),fibfun(j-1),fibfun(j-1)]),ar)
            vr = np.append(vr, np.array([fibfun(j),fibfun(j),fibfun(j-1)]))
            vrb = np.append(vrb, np.array([fibfun(j-1),0,fibfun(j-1)]))

        vr = np.append(vr, fibfun(i-1))
    else:
        raise ValueError("Topology type not implemented yet")

    #ac = np.array(ac)
    #vc = np.array(vc)
    #vcb = np.array(vcb)
    #ar = np.array(ar)
    #vr = np.array(vr)
    #vrb = np.array(vrb)

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
    print(res3)

    res4 = generate_topology("Cockcroft-walton",3,1)
    print(res4)

    res5 = generate_topology("Doubler",2,1)
    print(res5)

    res6 = generate_topology("Fibonacci",5,1)
    print(res6)
