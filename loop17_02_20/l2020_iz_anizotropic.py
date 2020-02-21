from numba import jit
import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import warnings
from math import exp,sqrt
# warnings.filterwarnings('ignore')

#SINGLE ITERATON
@jit
def iter(M,a1,a2,a3,a4,a5,imax,jbeg,jend,i1,i2,w):
    # print ('!!!!!')
    for i in range(0,i1):
        im1 = abs( i - 1 )
        for j in range(jbeg[i],jend[i]):
            jm1 = abs(j - 1)
            M[i][j] = M[i][j]-w *(a1[i][j]*M[i+1][j]+a2[i][j]*M[im1][j]+a3[i]*M[i][j]+a4[i][j]*M[i][j+1]+a5[i][j]*M[i][jm1])/a3[i]

    for i in range(i1,i2+1):
        im1 = abs( i - 1 )
        for j in range(jbeg[i]+1,jend[i]):
            jm1 = abs(j - 1)
            M[i][j] = M[i][j]-w *(a1[i][j]*M[i+1][j]+a2[i][j]*M[im1][j]+a3[i]*M[i][j]+a4[i][j]*M[i][j+1]+a5[i][j]*M[i][jm1])/a3[i]

    for i in range(i2+1,imax):
        im1 = abs( i - 1 )
        for j in range(jbeg[i],jend[i]):
            jm1 = abs(j - 1)
            M[i][j] = M[i][j]-w *(a1[i][j]*M[i+1][j]+a2[i][j]*M[im1][j]+a3[i]*M[i][j]+a4[i][j]*M[i][j+1]+a5[i][j]*M[i][jm1])/a3[i]
    return M

#FLOWS
    #Small
        #F corresponds to M
        #Ff corresponds to E
def flow_small(R, F,Ff,Lr,Lz,dr,dz):
    # from Energies import e as E
    # from Load_Energies import k, ee, ek
    from math import exp
    Pi=3.14

    ie = int(Lr/dr)
    je = int(Lz/dz)

    sum = 0
    for j in range(je+1):
        sum += exp(-Ff[ie][j])*(F[ie+1][j]-F[ie-1][j])
    sum = sum - 0.5*exp(-Ff[ie][0])*(F[ie+1][0]-F[ie-1][0])-0.5*exp(-Ff[ie][je])*(F[ie+1][je]-F[ie-1][je])
    Ie = Pi*Lr*dz/dr*sum
    sum = 0
    for i in range(ie+1):
        sum += i*exp(-Ff[i][je])*(F[i][je+1] - F[i][je-1])
    sum = sum - 0.5*ie*exp(-Ff[ie][je]) *(F[ie][je+1]  - F[ie][je-1])
    Ie = Ie + Pi*dr*dr/dz*sum
    Ie = Ie/(Pi*R)
    return Ie


    #Big
        #Y corresponds to M
        #Ff corresponds to E
def flow_big(R,L,Y,Ff,ie0, ie, je,dr,dz):
    from math import exp
    Pi=3.14

    sum = 0  # Integration along outer surface
    for j in range(je+1):
         sum += exp(-Ff[ie][j])*(Y[ie+1][j] - Y[ie-1][j])
    sum = sum-0.5*exp(-Ff[ie][0])*(Y[ie+1][0]-Y[ie-1][0])-0.5*exp(-Ff[ie][je])*(Y[ie+1][je]-Y[ie-1][je])
    Ie = Pi*(R-L+ie*dr)*dz/dr*sum

    sum = 0  # Integration along inner surface
    for j in range(je+1):
         sum += exp(-Ff[ie0][j])*(Y[ie0+1][j] - Y[ie0-1][j])
    sum = sum - 0.5*exp(-Ff[ie0][0]) *(Y[ie0+1][0]-Y[ie0-1][0])- 0.5*exp(-Ff[ie0][je])*(Y[ie0+1][je]-Y[ie0-1][je])
    Ie = Ie - Pi*(R-L+ie0*dr)*dz/dr*sum

    sum = 0  #Integral S_top
    for i in range(ie0,ie+1):
        sum += (R-L+i*dr)*exp(-Ff[i][je])*(Y[i][je+1] - Y[i][je-1])
    sum = sum-0.5*(R-L+ie0*dr)*exp(-Ff[ie0][je])*(Y[ie0][je+1]-Y[ie0][je-1])-0.5*(R-L+ie*dr)*exp(-Ff[ie][je])*(Y[ie][je+1]-Y[ie][je-1])

    Ie = Ie + Pi*dr/dz*sum
    Ie = Ie/(Pi*R)
    return(Ie)
def calc_coeff(Re,E,E_r,E_z,e, er , ez ,R, k,ee,ek,dee,dek,P,i1,i2,imax,jmax,jbeg,dr,dz):
    Ff = np.zeros( [imax+1,jmax+1] )
    dFr = np.zeros( [imax+1,jmax+1] )
    dFz = np.zeros( [imax+1,jmax+1] )
    Ffiz = np.zeros( [imax+1,jmax+1] )
    dFriz = np.zeros( [imax+1,jmax+1] )
    dFziz = np.zeros( [imax+1,jmax+1] )
    for i in range(imax+1):
        r=max(0,R-Re)+i*dr
        for j in range(jmax+1):
            z=j*dz
            if i1<=i<=i2 and j<jbeg[i]:
                Ff[i][j]=0
                dFr[i][j]=0
                dFz[i][j]=0

                Ffiz[i][j]=0
                dFriz[i][j]=0
                dFziz[i][j]=0
            else:
                Ff[i][j]=E(r,z,R,k,ee,ek,P)
                dFr[i][j]=E_r( r, z, R, k,ee,ek,dee,dek,P)
                dFz[i][j]=E_z( r, z, R, k,ee,ek,dee,dek,P)

                Ffiz[i][j]=e(r,z,R,k,ee,ek,P)
                dFriz[i][j]=er( r, z, R, k,ee,ek,dee,dek,P)
                dFziz[i][j]=ez( r, z, R, k,ee,ek,dee,dek,P)
    return {'Ff':Ff,'dFr':dFr,'dFz':dFz,'Ffiz':Ffiz,'dFriz':dFriz,'dFziz':dFziz}


def LiebmanAniz(_type,Ri,Re,R, P=1., _w=1.95, dr=1.,dz=1.):

    start = time.time()
    print ('start with jit and CALC!: ',time.asctime())
    from Energies_an import E_an as E, E_an_r as E_r, E_an_z as E_z
    from Energies import e, er , ez
    from Load_Energies import k, ee, ek, dee, dek

    k = k()
    ee = ee()
    ek = ek()
    dee = dee()
    dek = dek()
    L=Re
    rc=Ri
    imax=int((min(Re,R)+Re)/dr)
    jmax=int(Re/dz)

    jend=[0]*(imax+1)
    jbeg=[0]*(imax+1)
    # jend={}
    # jbeg={}
    for i in range(imax+1):
        jend[i] = int((sqrt(Re*Re - (min(R,Re)-i*dr)*(min(R,Re)-i*dr)))/dr)


    #Initialization of matrixes

    M = np.zeros( [imax+1,jmax+1] )
    #anizotropic functions
    a1 = np.zeros( [imax+1,jmax+1] )
    a2 = np.zeros( [imax+1,jmax+1] )
    a4 = np.zeros( [imax+1,jmax+1] )
    a5 = np.zeros( [imax+1,jmax+1] )
    # Ff = np.zeros( [imax+1,jmax+1] )
    # dFr = np.zeros( [imax+1,jmax+1] )
    # dFz = np.zeros( [imax+1,jmax+1] )

    a3 = np.zeros(imax+1)
    #izotropic functions
    a1iz = np.zeros( [imax+1,jmax+1] )
    a2iz = np.zeros( [imax+1,jmax+1] )
    a4iz = np.zeros( [imax+1,jmax+1] )
    a5iz = np.zeros( [imax+1,jmax+1] )
    # Ffiz = np.zeros( [imax+1,jmax+1] )
    # dFriz = np.zeros( [imax+1,jmax+1] )
    # dFziz = np.zeros( [imax+1,jmax+1] )

    a3iz = np.zeros(imax+1)


#**********************************************
#TODO
#Find left and right border of capture contour
    i1 = int((min(R,Re)-Ri)/dr)+1
    i2 = int((min(R,Re)+Ri)/dr)
    for i in range(imax+1):
        if 0<=i<=i1:
            jbeg[i]=0
        if i1<i<=i2:
            jbeg[i] = int(sqrt((Ri/dr)*(Ri/dr)-(i-min(R,Re)/dr)*(i-min(R,Re)/dr)))
        if i2<i:
            jbeg[i]=0
    jbeg=tuple(jbeg)
    jend=tuple(jend)
    # print 'jbeg '
    # print jbeg
    N_points=0
    for i in range(imax+1):
        for j in range(jbeg[i],jend[i]):
            N_points+=1




#Initial matrix

    for i in range(imax+1):
        r = max(0,R-Re)+i*dr
        # print jend[i]
        for j in range(jmax+1):
            z=j*dz
            M[i][j]=0.
            if j>=jend[i]:
                if _type=='exp':
                    M[i][j]= exp(e(r,z,R,k,ee,ek,P))
                if _type=='one':
                    M[i][j]= 1
# Definition of the coefficients of interaction
    print ('Start fill the matrixes')
    startfill=time.time()
    coefs=calc_coeff(Re,E,E_r,E_z,e, er , ez ,R, k,ee,ek,dee,dek,P,i1,i2,imax,jmax,jbeg,dr,dz)
    Ff=coefs['Ff'];dFr=coefs['dFr'];dFz=coefs['dFz']
    Ffiz=coefs['Ffiz'];dFriz=coefs['dFriz'];dFziz=coefs['dFziz']

#Find coefficients for yDot
    print(f'Time for coefficient calculation: {time.time()-startfill}')
    if R<=Re:
        a3[0]=-6.
        a3iz[0]=-6.
        for i in range(1,imax+1):
            a3[i]=-4.
            a3iz[i]=-4.
            for j in range(jmax+1):
                a1[i][j] = 1. + (1/(i*dr) - dFr[i][j])*dr/2
                a2[i][j] = 1. - (1/(i*dr) - dFr[i][j])*dr/2
                a4[i][j] = 1. - dFz[i][j]*dr/2
                a5[i][j] = 1. + dFz[i][j]*dr/2

                a1iz[i][j] = 1. + (1/(i*dr) - dFriz[i][j])*dr/2
                a2iz[i][j] = 1. - (1/(i*dr) - dFriz[i][j])*dr/2
                a4iz[i][j] = 1. - dFziz[i][j]*dr/2
                a5iz[i][j] = 1. + dFziz[i][j]*dr/2

        for j in range(jmax+1):
            a1[0][j] = 2. - dFr[0][j]*dr/2
            a2[0][j] = 2. + dFr[0][j]*dr/2
            a4[0][j] = 1. - dFz[0][j]*dr/2
            a5[0][j] = 1. + dFz[0][j]*dr/2

            a1iz[0][j] = 2. - dFriz[0][j]*dr/2
            a2iz[0][j] = 2. + dFriz[0][j]*dr/2
            a4iz[0][j] = 1. - dFziz[0][j]*dr/2
            a5iz[0][j] = 1. + dFziz[0][j]*dr/2
    else:
        for i in range(0,imax+1):
            a3[i]=-4.
            a3iz[i]=-4.
            for j in range(jmax+1):
                a1[i][j] = 1. + (1/(R-Re+i*dr) - dFr[i][j])*dr/2
                a2[i][j] = 1. - (1/(R-Re+i*dr) - dFr[i][j])*dr/2
                a4[i][j] = 1. - dFz[i][j]*dr/2
                a5[i][j] = 1. + dFz[i][j]*dr/2

                a1iz[i][j] = 1. + (1/(R-Re+i*dr) - dFriz[i][j])*dr/2
                a2iz[i][j] = 1. - (1/(R-Re+i*dr) - dFriz[i][j])*dr/2
                a4iz[i][j] = 1. - dFziz[i][j]*dr/2
                a5iz[i][j] = 1. + dFziz[i][j]*dr/2
    print (f'Finish fill the matrixes, T={time.time()-startfill}')
    # Define contours
    Lz1 = max(3.*rc,L/10.)
    Lr1 = sqrt(L*L-Lz1*Lz1)-Lz1+min(R,L)
    Lr01= min(R,L)-sqrt(L*L - Lz1*Lz1)+Lz1 #For R>L
    Lr2 = min(R,L)+Lz1
    Lr02 = L-4.*rc   #For R>L


    if(R<L):
         Lz2 = sqrt(L*L - R*R)-Lz1
    else:
        Lz2 = min(sqrt(L*L-(Lr02-L)*(Lr02-L)), sqrt(L*L-(Lr2-L)*(Lr2-L)))-Lz1

    ie1 =  int(Lr1/dr)
    ie01 = int(Lr01/dr)
    je1 =  int(Lz1/dz)
    ie02 = int(Lr02/dr)
    ie2 =  int(Lr2/dr)
    je2 =  int(Lz2/dz)

    dev=100
    e=1
    dev_old=0
    count=0
    dev=10
    fl=0
    print ('start EMPTY iterating, T=',time.time()-start)
    startempty=time.time()
    iter(M,a1iz,a2iz,a3iz,a4iz,a5iz,0,jbeg,jend,0,0,_w)
    print (f'Time for empty={time.time()-startempty}')
    print ('start IZotropic iterating, T=',time.time()-start)

    startisotropic=time.time()
    while  (dev>0.0001 or e>0.000001) and count<10000:
        dev_old=dev
        startiter=time.time()
        M=iter(M,a1iz,a2iz,a3iz,a4iz,a5iz,imax,jbeg,jend,i1,i2,_w)
        if count<5:
            print (f'Time for iter: {time.time()-startiter}')
        count+=1
        if R<=L:
            flow1=flow_small(R, M,Ffiz,Lr1,Lz1,dr,dz)
            flow2=flow_small(R, M,Ffiz,Lr2,Lz2,dr,dz)
        else:
            flow1=flow_big(R,L, M,Ffiz,ie01, ie1, je1,dr,dz)
            flow2=flow_big(R,L, M,Ffiz,ie02, ie2, je2,dr,dz)


        if count>100:# and count%100==0:

            fl=(flow1+flow2)/2

            dev=abs(flow1-flow2)/fl
            e=abs(dev-dev_old)
    # except:
    #     print ('!count={}, flow={},e={}, dev={},w={},T={}'.format(count,fl,e,dev, _w,time.time()-start))
    print ('FINISH IZotropic iterating with jit, T={}, flow={}'.format(time.time()-start, fl))
    print ('!count={}, flow={},e={}, dev={},w={},T={}'.format(count,fl,e,dev, _w,time.time()-start))
    print (f'Time for isotropic={time.time()-startisotropic}')

    print ('****************************************\n')
    for i in range(imax+1):
        r = max(0,R-Re)+i*dr
        for j in range(jmax+1):
            z=j*dz
            if j>=jend[i]:
                if _type=='exp':
                    M[i][j]= exp(E(r,z,R,k,ee,ek,P))
                if _type=='one':
                    M[i][j]= 1
    _w=1.95
    e=1
    dev_old=0
    count=0
    dev=10
    fl=0
    # try:
    startaniz=time.time()
    while  ( e>0.0000001) and count<10000:#DELETED dev>0.001 or
        dev_old=dev
        tt=time.time()
        M=iter(M,a1,a2,a3,a4,a5,imax,jbeg,jend,i1,i2,_w)
        # print (f'Time for iter: {time.time()-tt}')
        count+=1
        if R<=L:
            flow1=flow_small(R, M,Ff,Lr1,Lz1,dr,dz)
            flow2=flow_small(R, M,Ff,Lr2,Lz2,dr,dz)
        else:
            flow1=flow_big(R,L, M,Ff,ie01, ie1, je1,dr,dz)
            flow2=flow_big(R,L, M,Ff,ie02, ie2, je2,dr,dz)


        if count>100:# and count%100==0:

            fl=(flow1+flow2)/2

            dev=abs(flow1-flow2)/fl
            e=abs(dev-dev_old)
        if count%100==0:
            fl=(flow1+flow2)/2
            print ('ANIZ!count={}, flow={},e={}, dev={},w={},T={}'.format(count,fl,e,dev, _w,time.time()-start))
            # p=M[:,0]
            # z=np.array(range(len(p)))
            # plt.plot(z,p)
            # plt.show()
    # except:
    print ('Finish !count={}, flow={},e={}, dev={},w={},T={}'.format(count,fl,e,dev, _w,time.time()-start))
    print (f'Time for anizotropic={time.time()-startaniz}')





    if R<=L:
        flow1=flow_small(R, M,Ff,Lr1,Lz1,dr,dz)
        flow2=flow_small(R, M,Ff,Lr2,Lz2,dr,dz)
    else:
        flow1=flow_big(R,L, M,Ff,ie01, ie1, je1,dr,dz)
        flow2=flow_big(R,L, M,Ff,ie02, ie2, je2,dr,dz)

    fl=(flow1+flow2)/2
    dev=abs(flow1-flow2)/fl
    t=time.time()-start
    res={'sol':M,'flow':fl,'dev':dev,'N_points':N_points,'count':count, 'e':e,'i_max':imax,'j_max':jmax,'t':t}
    print ('DONE with ',time.time()-start)
    p=M[:,0]
    z=np.array(range(len(p)))
    plt.plot(z,p)
    plt.show()
    return res

def main():
    LiebmanAniz('one',5,55,3000, P=-12., _w=1.95, dr=1.,dz=1.)

if __name__ == '__main__':
    main()
