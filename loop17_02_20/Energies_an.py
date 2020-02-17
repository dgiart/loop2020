
# path="D:\Babich\Science\Ostapchuk\LOOP_New_Try_05_2018\Python\\"
# import sys
# sys.path.append(path)
# from closest import closest

#from Load_Energies import k, ee, ek, dee, dek


'''m= (R - r) * (R - r) + r * r
p = ((R + r) * (R + r) + z * z)  0.5
t = R * R - r * r - z * z

k_2 = 4 * R * r / p ^ 2'''


def closest(L,x):
    dist=abs(L[0]-x)
    closest=L[0]
    for el in L:
        if dist>abs(el-x):
            dist=abs(el-x)
            closest=el
    return closest#, dist
def I01(r, z, R,K,ee,ek):
    # from math import sqrt
    _k = float(( 4. * r * R / float((r + R) * (r + R) + z * z) )**0.5 )
    # _k = float( sqrt( 4. * r * R / float((r + R) * (r + R) + z * z) ) )
    k=closest(K,_k)
    # print 'k=',_k
    out = 1 / ( ( z * z + (r + R) * (r + R) )**0.5)
    # out = 1 / ( sqrt( z * z + (r + R) * (r + R) ))
    # print 'out=',out
    ins = (R * R - r * r - z * z) / float(z * z + (R - r) * (R - r))
    # print 'ins=',ins
    EE =ee[k]
    EK =ek[k]
    return float( out * (ins * EE + EK) )
def E_an(r,z,R,k,ee,ek,P):

    a1 = 0.234339408368
    a2 = -0.305925877794
    E=-2*P*(a1*I01(r, z/(1.6679), R,k,ee,ek)-a2*I01(r, z/(0.63166), R,k,ee,ek))
    return E

def E_an_r(r,z,R,K,ee,ek,dee,dek, P):
    m = float((R - r) * (R - r) + z * z)
    p = float(((R + r) * (R + r) + z * z) ** 0.5)
    t = float(R * R - r * r - z * z)
    _k = (4 * R * r / p ** 2)**0.5
    # print 'm={}, p={},t={}'.format(m,p,t)
    k2r=4*R*(1-2*r*(r+R)/(p*p))/(p*p)
    # print 'k2z=',k2z
    # print 'up=',(8*z*r*R)
    # print 'down=',((z*z+(r+R)*(r+R))*(z*z+(r+R)*(r+R)))
    k = closest(K, _k)
    _ee = ee[k]
    _ek = ek[k]
    _dee = dee[k]
    _dek=dek[k]
    # print 'I01(r, z, R,K,ee,ek)*z/(p*p)=',I01(r, z, R,K,ee,ek)*z/(p*p)
    # print '(4*z*R*(R-r)=',(4*z*R*(R-r))
    # print 'I01(r, z, R,K,ee,ek)=',I01(r, z, R,K,ee,ek)
    # print '(p*m*m)=',(p*m*m)
    # print '(4*z*R*(R-r)*I01(r, z, R,K,ee,ek))/(p*m*m)=',(4*z*R*(R-r)*I01(r, z, R,K,ee,ek))/(p*m*m)
    t1=-(R+r)*(t*_ee/m+_ek)/(p*p*p)
    t2=-(2*r*m-2*(R-r)*t)*_ee/(p*m*m)
    t3=k2r*(t*_dee/m+_dek)/p
    I_r=t1+t2+t3
    return -2*P*I_r
def I01_z(r,z,R,K,ee,ek,dee,dek):
    m = float((R - r) * (R - r) + z * z)
    p = float(((R + r) * (R + r) + z * z) ** 0.5)
    t = float(R * R - r * r - z * z)
    _k = (4 * R * r / p ** 2)**0.5
    # print 'm={}, p={},t={}'.format(m,p,t)
    k2z=-(8.*z*r*R)/float((z*z+(r+R)*(r+R))*(z*z+(r+R)*(r+R)))
    # print 'k2z=',k2z
    # print 'up=',(8*z*r*R)
    # print 'down=',((z*z+(r+R)*(r+R))*(z*z+(r+R)*(r+R)))
    k = closest(K, _k)
    _ee = ee[k]
    _ek = ek[k]
    _dee = dee[k]
    _dek=dek[k]
    # print 'I01(r, z, R,K,ee,ek)*z/(p*p)=',I01(r, z, R,K,ee,ek)*z/(p*p)
    # print '(4*z*R*(R-r)=',(4*z*R*(R-r))
    # print 'I01(r, z, R,K,ee,ek)=',I01(r, z, R,K,ee,ek)
    # print '(p*m*m)=',(p*m*m)
    # print '(4*z*R*(R-r)*I01(r, z, R,K,ee,ek))/(p*m*m)=',(4*z*R*(R-r)*I01(r, z, R,K,ee,ek))/(p*m*m)
    I_z=-I01(r, z, R,K,ee,ek)*z/(p*p)-(4*z*R*(R-r)*_ee)/(p*m*m)+k2z*(t*_dee/m+_dek)/p
    return I_z







def E_an_z(r,z,R,k,ee,ek,dee,dek,P):
    a1 = 0.234339408368
    a2 = -0.305925877794
    E_z = -2*P * (a1 * I01_z( r, z/1.6679, R,k,ee,ek,dee,dek )/1.6679 - a2 * I01_z( r, z/0.63166, R,k,ee,ek,dee,dek) /0.63166)
    return E_z
