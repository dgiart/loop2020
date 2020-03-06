from math import sqrt


def closest(array,value):
    '''Given an ``array`` , and given a ``value`` , returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively.'''
    n = len(array)
    # if (value < array[0]):
    #     return -1
    # elif (value > array[n-1]):
    #     return n
    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return array[jl]



def e(r,z,R,k,ee,ek,P):
    c11=1.554; c12=0.672; c13=0.646; c33=1.725; c44=0.363; c66=(c11-c22)/2
    c11=1.554; c12=0.672; c13=0.646; c33=1.725; c44=0.363; c66=(c11-c12)/2
    v1=(-(c13*c13+2*c13*c44-c11*c33)+sqrt((c13*c13+2*c13*c44-c11*c33)*(c13*c13+2*c13*c44-c11*c33)-4*(c11*c44*c44*c33)))/(2*c11*c44)
    v2=(-(c13*c13+2*c13*c44-c11*c33)-sqrt((c13*c13+2*c13*c44-c11*c33)*(c13*c13+2*c13*c44-c11*c33)-4*(c11*c44*c44*c33)))/(2*c11*c44)
