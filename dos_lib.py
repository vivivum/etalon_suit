# Function definition is here
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import time

# Version 1.0

print('###################### FJB&DOS python library #########################')
print('## David Orozco (oroco@iaa.es) and Francisco Bailén (fbailen@iaa.es) ##')
print('## Instituto de Astrofísica de Granada (IAA-CSIC)                    ##')
print('## Current: Version 1.0 released on 6th April, 2018                  ##')
print('#######################################################################')

    #Niobato de Litio parametros por defecto (no birefringencia)
no = 2.2926446 ##Ordinary refractive index
n3 = 2.2 #n3
h = 251.811e-6 #Etalon h [Meters]
Reflectance = 0.925 #Reflectance
Absorptance = 0.0 #Absorptance

def array(*args, **kwargs):
    kwargs.setdefault("dtype", np.float32)
    return np.zeros(*args, **kwargs)

# def fineza(R):
#     F = 4.0*R/(1.0-R)**2.0
#     return F
# "This calculates the finesse"
FN = lambda R: 4.0*R/(1.0-R)**2.0

#Transmittance
tau = lambda A,R: (1-A/(1-R))**2

#Calculates de cos(theta_o) angle from incidence angle theta_i"
# See eqn.48 in Paper I
# output is in radians
cos_theta_o = lambda theta_i:np.sqrt( 1.0 - (np.sin(theta_i)/no)**2 )

# "This calculates the phase shift between rays for the ordinary axis"
# lambda and h in SI
# output is in rad
delta = lambda theta_i, lda , theta = 2*pi: 4.0*pi*no*h*cos_theta_o(theta_i)/lda + 2*np.cos(theta)

def phi(theta_i, lda , theta_3):
    # "This calculates the retardance between the ord. and ext. rays"
    # lambda in nm and h in mm
    # output is in rad
    n = (no + n3) / 2.0
    theta_prime_rad = np.arcsin( np.sin( theta_i )/n )
    phi_val = ( 4.0 * pi * h * (n3 - no)
        * np.sin( theta_prime_rad - theta_3 )**2
        / ( lda * np.cos(theta_prime_rad) ) )
    return phi_val

def FPM(theta_i, lda , theta_3, r_angle, doprint='false'):
    #etalon mueller matrix
    #All angles in radians
    try:
        elsize = lda.size
        files, columns, length = 4, 4, elsize; #(files,columns)
        Matrix = array([files,columns,length])
    except AttributeError as err:
        print('Single wavelength')
        elsize = 1
        files, columns = 4, 4 #(files,columns)
        Matrix = array([files,columns])

    #print('Before: ',no,n3,h,Reflectance,Absorptance,theta_i,theta_3,r_angle)

    #start_time = time.time()

    phi_angle = phi(theta_i, lda , theta_3)
    delta_angle = delta(theta_i, lda)
    M1 = np.zeros((4, 4, elsize)) #retarder
    M2 = np.zeros((4, 4, elsize)) #mirror

    psi = ( (1+FN(Reflectance)*np.sin( delta_angle/2.0 )**2.0)
        * (1+FN(Reflectance)*np.sin( (delta_angle + phi_angle )/2.0 )**2.0 ) )
    alpha = ( FN(Reflectance)/2.0 * (np.sin( (delta_angle + phi_angle )/2.0 )**2.0 +
            np.sin( delta_angle/2.0 )**2.0) )
    beta =  ( FN(Reflectance)/2.0 * (np.sin( (delta_angle + phi_angle )/2.0 )**2.0 -
            np.sin( delta_angle/2.0 )**2.0) )
    gamma = ( FN(Reflectance)/4.0 * ( Reflectance*np.cos( phi_angle/2.0 ) -
            2.0*np.cos( delta_angle + phi_angle/2.0 ) ) )
    sigma = ( -FN(Reflectance)/4.0 * Reflectance * np.sin( phi_angle/2.0 ) )

    M1[0,0,:] = (1.0 - Reflectance)**2.0
    M1[1,1,:] = (1.0 - Reflectance)**2.0
    M1[2,2,:] = np.cos(phi_angle/2)
    M1[2,3,:] = -np.sin(phi_angle/2)
    M1[3,2,:] = np.sin(phi_angle/2)
    M1[3,3,:] = np.cos(phi_angle/2)
    M1 = M1/(1.0 - Reflectance)**2.0
##
    M2[0,0,:] = alpha
    M2[0,1,:] = beta
    M2[1,0,:] = beta
    M2[1,1,:] = alpha
    M2[2,2,:] = gamma
    M2[2,3,:] = sigma
    M2[3,2,:] = -sigma
    M2[3,3,:] = gamma

    MM = (M1 + M2) * tau(Absorptance,Reflectance)/psi

    MB = MM
    MM[:,:,:] = 0.0

    a,b,c,d = abcd(lda,theta_i,theta_3)
    #plt.close()
    MM[0,0,:] = a.real
    MM[1,1,:] = a.real
    MM[0,1,:] = b.real
    MM[1,0,:] = b.real
    MM[2,2,:] = c.real
    MM[3,3,:] = c.real
    MM[2,3,:] = d.real
    MM[3,2,:] = -d.real
    MR1 = ROT(r_angle)
    MR2 = ROT(-r_angle)
    if r_angle != 0:
        for k in range(elsize):
            MM[:,:,k] = np.matmul(MR2,np.matmul(MM[:,:,k],MR1))

    if doprint.upper() == 'TRUE':
        lda2 = (lda-617.33356e-9)*1e9
        plt.subplot(3, 5, 1)
        plt.plot(lda2,MM[0,0,:])
        plt.subplot(3, 5, 2)
        plt.plot(lda2,MM[0,1,:])
        plt.subplot(3, 5, 3)
        plt.plot(lda2,MM[2,2,:])
        plt.subplot(3, 5, 4)
        plt.plot(lda2,MM[2,3,:])
        plt.subplot(3, 5, 5)
        plt.plot(lda2,MM[1,3,:])
        plt.subplot(3, 5, 6)
        plt.plot(lda2,MB[0,0,:])
        plt.subplot(3, 5, 7)
        plt.plot(lda2,MB[0,1,:])
        plt.subplot(3, 5, 8)
        plt.plot(lda2,MB[2,2,:])
        plt.subplot(3, 5, 9)
        plt.plot(lda2,MB[2,3,:])
        plt.subplot(3, 5, 10)
        plt.plot(lda2,MB[1,3,:])
        plt.subplot(3, 5, 11)
        plt.plot(MM[0,0,:]-MB[0,0,:])
        plt.subplot(3, 5, 12)
        plt.plot(MM[0,1,:]-MB[0,1,:])
        plt.subplot(3, 5, 13)
        plt.plot(MM[2,2,:]-MB[2,2,:])
        plt.subplot(3, 5, 14)
        plt.plot(MM[2,3,:]-MB[2,3,:])
        plt.subplot(3, 5, 15)
        plt.plot(MM[1,3,:]-MB[1,3,:])
        plt.show()

    #print('After: ',no,n3,h,Reflectance,Absorptance,theta_i,theta_3,r_angle)
#    print("--- %s seconds ---" % (time.time() - start_time))

    return MB

#FJB version (numerically)

def deltao(wave,theta):
    deltao=(4*pi*h/wave)*np.sqrt(no**2-1+(np.cos(theta))**2)
    return deltao

def deltae(wave,theta,theta3):
    n=(no+n3)/2
    thetat=np.arcsin(np.sin(theta)/n)
    phi=(4*pi*h*n)*(n3-no)*\
    (np.sin(thetat-theta3))**2/(wave*np.sqrt(n**2-(np.sin(theta))**2))
    deltae=phi+deltao(wave,theta)
    return deltae

def H11(wave,theta):
    dto=deltao(wave,theta)
    H11=(np.sqrt(tau(Absorptance,Reflectance))/(1-Reflectance))*(1-Reflectance*np.exp(-1j*dto))*\
    np.exp(1j*dto/2)/(1+FN(Reflectance)*(np.sin(dto/2))**2)
    return H11

def H22(wave,theta,theta3):
    dte=deltae(wave,theta,theta3)
    H22=(np.sqrt(tau(Absorptance,Reflectance))/(1-Reflectance))*(1-Reflectance*np.exp(-1j*dte))*\
    np.exp(1j*dte/2)/(1+FN(Reflectance)*(np.sin(dte/2))**2)
    return H22

def abcd(wave,theta,theta3):  #R
    if np.size(theta)==1:
        h11=H11(wave,theta)
        h22=H22(wave,theta,theta3)
        a=(1/2)*(h11*np.conj(h11)+h22*np.conj(h22))
        b=(1/2)*(h11*np.conj(h11)-h22*np.conj(h22))
        c=(1/2)*(h22*np.conj(h11)+h11*np.conj(h22))
        d=(1j/2)*(h22*np.conj(h11)-h11*np.conj(h22))
    else:
        a=np.zeros((np.size(wave),np.size(theta)))
        b=np.copy(a)
        c=np.copy(a)
        d=np.copy(a)
        i=-1
        for theta_i in theta:
            i+=1
            h11=H11(wave,theta_i)
            h22=H22(wave,theta_i,theta3)
            a[i,:]=(1/2)*np.real(h11*np.conj(h11)+h22*np.conj(h22))
            b[i,:]=(1/2)*np.real(h11*np.conj(h11)-h22*np.conj(h22))
            c[i,:]=(1/2)*np.real(h22*np.conj(h11)+h11*np.conj(h22))
            d[i,:]=np.real((1j/2)*(h22*np.conj(h11)-h11*np.conj(h22)))
    return a,b,c,d

def abcd2(wave,theta,theta3):
    if np.size(theta)==1:
        h11=H11(wave,theta)
        h22=H22(wave,theta,theta3)
        a=(1/2)*(h11*np.conj(h11)+h22*np.conj(h22))
        b=(1/2)*(h11*np.conj(h11)-h22*np.conj(h22))
        c=(1/2)*(h22*np.conj(h11)+h11*np.conj(h22))
        d=(1j/2)*(h22*np.conj(h11)-h11*np.conj(h22))
    else:
        WVL,THETA=np.meshgrid(wave,theta)
        h11=H11(WVL,THETA)
        h22=H22(WVL,THETA,theta3)
        a=(1/2)*(h11*np.conj(h11)+h22*np.conj(h22))
        b=(1/2)*(h11*np.conj(h11)-h22*np.conj(h22))
        c=(1/2)*(h22*np.conj(h11)+h11*np.conj(h22))
        d=(1j/2)*(h22*np.conj(h11)-h11*np.conj(h22))
    return a,b,c,d



def LCVR(theta, delta):
    #LCRV mueller matrix
    # ;-------------------------------------------------------------------------------
    # ;quarter-wave plate, fast axis horizontal  THETA = 0, DELTA = 90
    # ;quarter-wave plate, fast axis vertical  THETA = 90, DELTA = 90
    # ;angles are measured positive counterclockwise
    # ;z axis positive direction is along the direction of light propagation
    # ;this is a right handed system
    # ;
    # ;   y
    # ;    |
    # ;    |
    # ;    |
    # ;    ---------->  X
    # ;

    deltar = np.radians(delta)
    thethar = np.radians(theta)
    filas = 4
    columnas = 4
    LCVRM = np.zeros([filas,columnas])
    DLCVR = np.zeros([filas,columnas,2])

    c2 = np.cos(2*thethar)
    s2 = np.sin(2*thethar)

    LCVRM[0,0] = 1.0
    LCVRM[1,1] = c2**2.0+s2**2.0*np.cos(deltar)
    LCVRM[2,1] = c2*s2*(1.0-np.cos(deltar))
    LCVRM[3,1] = s2*np.sin(deltar)
    LCVRM[1,2] = c2*s2*(1.0-np.cos(deltar))
    LCVRM[2,2] = s2**2.0+c2**2.0*np.cos(deltar)
    LCVRM[3,2] = -c2*np.sin(deltar)
    LCVRM[1,3] = -s2*np.sin(deltar)
    LCVRM[2,3] = c2*np.sin(deltar)
    LCVRM[3,3] = np.cos(deltar)

    DLCVR[1,1,1] = -s2**2.0*np.sin(deltar)
    DLCVR[2,1,1] = c2*s2*np.sin(deltar)
    DLCVR[3,1,1] = s2*np.cos(deltar)
    DLCVR[1,2,1] = c2*s2*np.sin(deltar)
    DLCVR[2,2,1] = -c2**2.0*np.sin(deltar)
    DLCVR[3,2,1] = -c2*np.cos(deltar)
    DLCVR[1,3,1] = -s2*np.cos(deltar)
    DLCVR[2,3,1] = c2*np.cos(deltar)
    DLCVR[3,3,1] = -np.sin(deltar)

    DLCVR[1,1,0] = -4.*c2*s2+4.*s2*c2*np.cos(deltar)
    DLCVR[2,1,0] = 2.*(c2*c2-s2*s2)*(1.0-np.cos(deltar))
    DLCVR[3,1,0] = 2.*c2*np.sin(deltar)
    DLCVR[1,2,0] = 2.*(c2*c2-s2*s2)*(1.0-np.cos(deltar))
    DLCVR[2,2,0] = 4.*s2*c2-4.*c2*s2*np.cos(deltar)
    DLCVR[3,2,0] = +2.*s2*np.sin(deltar)
    DLCVR[1,3,0] = -2.*c2*np.sin(deltar)
    DLCVR[2,3,0] = -2.*s2*np.sin(deltar)
    DLCVR[3,3,0] = 0.0

    DLCVR = DLCVR * pi/180.0

    return LCVRM

def POL_LIN(alpha):
    # ;OK
    # ; ALPHA IS = 0 WHEN THE POLARIZER IS HORIZONTAL!!!! Vertical => THETA = 90
    # ;
    # ;   y
    # ;    |
    # ;    |
    # ;    |
    # ;    ---------->  X
    # ;

    angle = np.radians(alpha)
    filas = 4
    columnas = 4
    #PL = np.zeros([filas,columnas],dtype=float64)
    PL = array([filas,columnas])
    c2 = np.cos(2*angle,dtype=np.float32)
    s2 = np.sin(2*angle,dtype=np.float32)

    PL[0,0] = 1.0
    PL[0,1] = c2
    PL[0,2] = s2
    PL[1,0] = c2
    PL[1,1] = c2**2
    PL[1,2] = s2*c2
    PL[2,0] = s2
    PL[2,1] = c2*s2
    PL[2,2] = s2**2

    # ;calcula la derivada
    # c4 = cos(4d0*angle)
    # s4 = sin(4d0*angle)
    #
    # DPOL_LIN=[[0d0 , -s2 , c2  , 0d0],$
    #           [-s2 , -s4 , c4  , 0d0],$
    #           [c2  , c4  , s4  , 0d0],$
    #           [0d0 , 0d0 , 0d0 , 0d0]] * !dpi/180d0

    return 0.5*PL

def ROT(alpha):
    # ROTATION MATRIX. Input in radians
    filas = 4
    columnas = 4
    ROT = array([filas,columnas])
    c2 = np.cos(2*alpha)
    s2 = np.sin(2*alpha)

    ROT[0,0] = 1.
    ROT[1,1] = c2
    ROT[1,2] = s2
    ROT[2,1] = -s2
    ROT[2,2] = c2
    ROT[3,3] = 1.

    return ROT


# # Function definition is here
# def printinfo( arg1, *vartuple ):
#    "This prints a variable passed arguments"
#    print ("Output is: ")
#    print (arg1)
#    for var in vartuple:
#       print (var)
#    return
# # Function definition is here
# sum = lambda arg1, arg2: arg1 + arg2
#
# # Now you can call sum as a function
# print ("Value of total : ", sum( 10, 20 ))
