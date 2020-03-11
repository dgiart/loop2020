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

def I01(r, z, R,K,ee,ek):
    from math import sqrt
    _k = float( sqrt( 4. * r * R / float((r + R) * (r + R) + z * z) ) )
    k=closest(K,_k)
    # print 'k=',_k
    out = 1. / ( sqrt( z * z + (r + R) * (r + R) ))
    # print 'out=',out
    ins = (R * R - r * r - z * z) / float(z * z + (R - r) * (R - r))
    # print 'ins=',ins
    EE =ee[k]
    EK =ek[k]
    return float( out * (ins * EE + EK) )

def e_an(r,z,R,eps,k,ee,ek,P):
    pi=3.1416
    c11=1.554; c12=0.672; c13=0.646; c33=1.725; c44=0.363; c66=(c11-c12)/2
    v1=(-(c13*c13+2*c13*c44-c11*c33)+sqrt((c13*c13+2*c13*c44-c11*c33)*(c13*c13+2*c13*c44-c11*c33)-4*(c11*c44*c44*c33)))/(2*c11*c44)
    v2=(-(c13*c13+2*c13*c44-c11*c33)-sqrt((c13*c13+2*c13*c44-c11*c33)*(c13*c13+2*c13*c44-c11*c33)-4*(c11*c44*c44*c33)))/(2*c11*c44)
    v3=c44/c66
    d=(eps*(c11+c12)-2*c13)/(c33-eps*c13)
    B1=(c13+v2*c11)/(2*c11*(c13+c44)*(v1-v2)*sqrt(v1))
    B2=(c13+v1*c11)/(2*c11*(c13+c44)*(v2-v1)*sqrt(v2))
    B=P*(c11+c12+d*c13)/(pi*(2+d))
    First=B*B1*(v1*(c13+c44)+eps*(c44-v1*c11))
    Second=B*B2*(v2*(c13+c44)+eps*(c44-v2*c11))
    return First*I01(r, z/(sqrt(v1)), R,k,ee,ek)+Second*I01(r, z/(sqrt(v2)), R,k,ee,ek)

def I_z(r,z,R,K,ee,ek,dee,dek):
    m = float((R - r) * (R - r) + z * z)
    p = float(((R + r) * (R + r) + z * z) ** 0.5)
    t = float(R * R - r * r - z * z)
    _k = (4 * R * r / p ** 2)**0.5
    k2z=-(8.*z*r*R)/float((z*z+(r+R)*(r+R))*(z*z+(r+R)*(r+R)))
    k = closest(K, _k)
    _ee = ee[k]
    _ek = ek[k]
    _dee = dee[k]
    _dek=dek[k]
    I_z=-I01(r, z, R,K,ee,ek)*z/(p*p)-(4*z*R*(R-r)*_ee)/(p*m*m)+k2z*(t*_dee/m+_dek)/p
    return I_z

def E_an_z(r,z,R,eps,k,ee,ek,dee,dek,P):
    pi=3.1416
    c11=1.554; c12=0.672; c13=0.646; c33=1.725; c44=0.363; c66=(c11-c12)/2
    v1=(-(c13*c13+2*c13*c44-c11*c33)+sqrt((c13*c13+2*c13*c44-c11*c33)*(c13*c13+2*c13*c44-c11*c33)-4*(c11*c44*c44*c33)))/(2*c11*c44)
    v2=(-(c13*c13+2*c13*c44-c11*c33)-sqrt((c13*c13+2*c13*c44-c11*c33)*(c13*c13+2*c13*c44-c11*c33)-4*(c11*c44*c44*c33)))/(2*c11*c44)
    v3=c44/c66
    d=(eps*(c11+c12)-2*c13)/(c33-eps*c13)
    B1=(c13+v2*c11)/(2*c11*(c13+c44)*(v1-v2)*sqrt(v1))
    B2=(c13+v1*c11)/(2*c11*(c13+c44)*(v2-v1)*sqrt(v2))
    B=P*(c11+c12+d*c13)/(pi*(2+d))
    First=B*B1*(v1*(c13+c44)+eps*(c44-v1*c11))
    Second=B*B2*(v2*(c13+c44)+eps*(c44-v2*c11))
    E_z =First*I_z(r, z/(sqrt(v1)), R,k,ee,ek,dee,dek)+Second*I_z(r, z/(sqrt(v2)), R,k,ee,ek,dee,dek)
    return E_z
