import numpy as np
import astropy.units as u


def r(x, y):
    return np.sqrt(x ** 2 + y ** 2)


def phi(x, y):
    if y >= 0:
        return np.arccos(x / r(x, y))
    return 2. * np.pi - np.arccos(x / r(x, y))


# X-field component

BX0 = (4.6e-6 * u.uG).value
thetaX0 = (49 * u.degree.to(u.rad))
rXc = (4.8 * u.kpc).value
rX = (2.9 * u.kpc).value
h = rXc * np.tan(thetaX0)


def H(x, y, z):
    return np.heaviside(np.sqrt(x ** 2 + y ** 2 + z ** 2) - 1, 0.)


def rp0(x, y, z):
    return r(x, y) * rXc / (rXc + np.absolute(z) / np.tan(thetaX0))


def rp1(x, y, z):
    return r(x, y) - np.absolute(z) / np.tan(thetaX0)


def BX(x, y, z):
    if r(x, y) > rXc + np.absolute(z) / np.tan(thetaX0):
        return H(x, y, z) * BX0 * np.exp(-rp1(x, y, z) / rX) * (rp1(x, y, z) / r(x, y))
    return H(x, y, z) * BX0 * np.exp(-rp0(x, y, z) / rX) * ((rp0(x, y, z) / r(x, y)) ** 2)


def thetaX(x, y, z):
    if r(x, y) > rXc + np.absolute(z) / np.tan(thetaX0):
        return thetaX0
    return np.pi / 2. - np.arctan(r(x, y) / (np.absolute(z) + h))


def BXr(x, y, z):
    if z >= 0:
        return BX(x, y, z) * np.cos(thetaX(x, y, z))
    return -BX(x, y, z) * np.cos(thetaX(x, y, z))


def BXx(x, y, z):
    return BXr(x, y, z) * np.cos(phi(x, y))


def BXy(x, y, z):
    return BXr(x, y, z) * np.sin(phi(x, y))


def BXz(x, y, z):
    return BX(x, y, z) * np.sin(thetaX(x, y, z))


# Halo component

hdisk = (0.4 * u.kpc).value
wdisk = (0.27 * u.kpc).value
Bn = (1.4e-6 * u.uG).value
Bs = (-1.1e-6 * u.uG).value
rn = (9.22 * u.kpc).value
rs = (16.7 * u.kpc).value
wh = (0.2 * u.kpc).value
z0 = (5.3 * u.kpc).value


def L(z, h, w):
    return (1 + np.exp(-2 * (np.absolute(z) - h) / w)) ** (-1)


def BH(x, y ,z):
    if z >= 0:
        return np.exp(-(np.absolute(z) / z0)) * L(z, hdisk, wdisk) * Bn * (1 - L(r(x, y), rn, wh))
    return np.exp(-(np.absolute(z) / z0)) * L(z, hdisk, wdisk) * Bs * (1 - L(r(x, y), rs, wh))


def BHx(x, y, z):
    return -BH(x, y, z) * np.sin(phi(x, y))


def BHy(x, y, z):
    return BH(x, y, z) * np.cos(phi(x, y))


def BHz(x, y, z):
    return 0.


#Disk component


bb = 1 / np.tan(np.pi / 2. - 11.5 * u.degree.to(u.rad))
Bi = [0.1 * 1e-6, 3.0 * 1e-6, -0.9 * 1e-6, -0.8 * 1e-6, -2.0 * 1e-6,
                -4.2 * 1e-6, 0.0 * 1e-6, 2.7 * 1e-6, 0.1 * 1e-6]  # uG mag field
ai = [-5.1, -6.3, -7.1, -8.3, -9.8, -11.4, -12.7, -15.5]  # kpc
aim2p = [i * np.exp(bb * -2 * np.pi) for i in ai]
ai2p = [i * np.exp(bb * 2 * np.pi) for i in ai]
a7m4p = ai[7] * np.exp(bb * -4 * np.pi)


def H0(x, y):
    return np.heaviside(np.sqrt(x ** 2 + y ** 2) - 5, 0)


def H20(x, y):
    return np.heaviside(20 - np.sqrt(x ** 2 + y ** 2), 0)


def Ai(x, y):
    if x == 0.:
        return (x + 1e-5) / (np.cos(phi(x + 1e-5, y) - np.pi) * np.exp(bb * (phi(x + 1e-5, y) - np.pi)))
    return x / (np.cos(phi(x, y) - np.pi) * np.exp(bb * (phi(x, y) - np.pi)))


def B0(x, y):
    if ai[0] < Ai(x, y) < aim2p[7]:
        return Bi[0] * 5 / r(x, y)
    return 0.
def B0m2p(x, y):
    if aim2p[0] < Ai(x, y) < a7m4p:
        return Bi[0] * 5 / r(x, y)
    return 0.
def B02p(x, y):
    if ai2p[0] < Ai(x, y) < ai[7]:
        return Bi[0] * 5 / r(x, y)
    return 0.
def B0tot(x, y):
    return H0(x, y) * (B0(x, y) + B0m2p(x, y) + B02p(x, y))


def B1(x, y):
    if ai[1] < Ai(x, y) < ai[0]:
        return Bi[1] * 5 / r(x, y)
    return 0.
def B1m2p(x, y):
    if aim2p[1] < Ai(x, y) < aim2p[0]:
        return Bi[1] * 5 / r(x, y)
    return 0.
def B12p(x, y):
    if ai2p[1] < Ai(x, y) < ai2p[0]:
        return Bi[1] * 5 / r(x, y)
    return 0.
def B1tot(x, y):
    return H0(x, y) * (B1(x, y) + B1m2p(x, y) + B12p(x, y))


def B2(x, y):
    if ai[2] < Ai(x, y) < ai[1]:
        return Bi[2] * 5 / r(x, y)
    return 0.
def B2m2p(x, y):
    if aim2p[2] < Ai(x, y) < aim2p[1]:
        return Bi[2] * 5 / r(x, y)
    return 0.
def B22p(x, y):
    if ai2p[2] < Ai(x, y) < ai2p[1]:
        return Bi[2] * 5 / r(x, y)
    return 0.
def B2tot(x, y):
    return H0(x, y) * (B2(x, y) + B2m2p(x, y) + B22p(x, y))


def B3(x, y):
    if ai[3] < Ai(x, y) < ai[2]:
        return Bi[3] * 5 / r(x, y)
    return 0.
def B3m2p(x, y):
    if aim2p[3] < Ai(x, y) < aim2p[2]:
        return Bi[3] * 5 / r(x, y)
    return 0.
def B32p(x, y):
    if ai2p[3] < Ai(x, y) < ai2p[2]:
        return Bi[3] * 5 / r(x, y)
    return 0.
def B3tot(x, y):
    return H0(x, y) * (B3(x, y) + B3m2p(x, y) + B32p(x, y))


def B4(x, y):
    if ai[4] < Ai(x, y) < ai[3]:
        return Bi[4] * 5 / r(x, y)
    return 0.
def B4m2p(x, y):
    if aim2p[4] < Ai(x, y) < aim2p[3]:
        return Bi[4] * 5 / r(x, y)
    return 0.
def B42p(x, y):
    if ai2p[4] < Ai(x, y) < ai2p[3]:
        return Bi[4] * 5 / r(x, y)
    return 0.
def B4tot(x, y):
    return H0(x, y) * (B4(x, y) + B4m2p(x, y) + B42p(x, y))


def B5(x, y):
    if ai[5] < Ai(x, y) < ai[4]:
        return Bi[5] * 5 / r(x, y)
    return 0.
def B5m2p(x, y):
    if aim2p[5] < Ai(x, y) < aim2p[4]:
        return Bi[5] * 5 / r(x, y)
    return 0.
def B52p(x, y):
    if ai2p[5] < Ai(x, y) < ai2p[4]:
        return Bi[5] * 5 / r(x, y)
    return 0.
def B5tot(x, y):
    return H0(x, y) * (B5(x, y) + B5m2p(x, y) + B52p(x, y))


def B6(x, y):
    if ai[6] < Ai(x, y) < ai[5]:
        return Bi[6] * 5 / r(x, y)
    return 0.
def B6m2p(x, y):
    if aim2p[6] < Ai(x, y) < aim2p[5]:
        return Bi[6] * 5 / r(x, y)
    return 0.
def B62p(x, y):
    if ai2p[6] < Ai(x, y) < ai2p[5]:
        return Bi[6] * 5 / r(x, y)
    return 0.
def B6tot(x, y):
    return H0(x, y) * (B6(x, y) + B6m2p(x, y) + B62p(x, y))


def B7(x, y):
    if ai[7] < Ai(x, y) < ai[6]:
        return Bi[7] * 5 / r(x, y)
    return 0.
def B7m2p(x, y):
    if aim2p[7] < Ai(x, y) < aim2p[6]:
        return Bi[7] * 5 / r(x, y)
    return 0.
def B72p(x, y):
    if ai2p[7] < Ai(x, y) < ai2p[6]:
        return Bi[7] * 5 / r(x, y)
    return 0.
def B7tot(x, y):
    return H0(x, y) * (B7(x, y) + B7m2p(x, y) + B72p(x, y))


def BSptot(x, y, z):
    return (1 - L(z, hdisk, wdisk)) * H20(x, y) * (B0tot(x, y) + B1tot(x, y) + B2tot(x, y) + B3tot(x, y)
                                                   + B4tot(x, y) + B5tot(x, y) + B6tot(x, y) + B7tot(x, y))
def BSptotx(x, y, z):
    return BSptot(x, y, z) * np.cos(np.pi / 2. + phi(x, y) - 11.5 * u.degree.to(u.rad))
def BSptoty(x, y, z):
    return BSptot(x, y, z) * np.sin(np.pi / 2. + phi(x, y) - 11.5 * u.degree.to(u.rad))


def Bring(x, y, z):
    if 3. < r(x, y) < 5.:
        return 0.1e-6 * (1 - L(z, hdisk, wdisk))
    return 0.
def Bringx(x, y, z):
    return Bring(x, y ,z) * np.cos(phi(x, y) + np.pi / 2.)
def Bringy(x, y, z):
    return Bring(x, y ,z) * np.sin(phi(x, y) + np.pi / 2.)


def Ballx(x, y, z):
    return BXx(x, y, z) + BHx(x, y, z) + BSptotx(x, y, z) + Bringx(x, y, z)
def Bally(x, y, z):
    return BXy(x, y, z) + BHy(x, y, z) + BSptoty(x, y, z) + Bringy(x, y, z)
def Ballz(x, y, z):
    return BXz(x, y, z) + BHz(x, y, z) + 0. + 0.


#
# RANDOM MAGNETIC FIELD COMPONENT
#

# random Halo component

B0hR = 4.68e-6
r0hR = 10.97
z0haloR = 2.84


def BhaloR(x, y, z):
    return B0hR * np.exp(-r(x, y) / r0hR) * np.exp(-z ** 2 / (2 * z0haloR ** 2)) * H20(x, y)


# random DISK

BiR = [10.81e-6, 6.96e-6, 9.59e-6, 6.96e-6, 1.96e-6, 16.34e-6, 37.29e-6, 10.35e-6, 7.63e-6]  # uG mag field
z0diskR = 0.61


def B0R(x, y):
    if ai[0] < Ai(x, y) < aim2p[7]:
        return BiR[0] * 5 / r(x, y)
    return 0.
def B0m2pR(x, y):
    if aim2p[0] < Ai(x, y) < a7m4p:
        return BiR[0] * 5 / r(x, y)
    return 0.
def B02pR(x, y):
    if ai2p[0] < Ai(x, y) < ai[7]:
        return BiR[0] * 5 / r(x, y)
    return 0.
def B0totR(x, y):
    return H0(x, y) * (B0R(x, y) + B0m2pR(x, y) + B02pR(x, y))


def B1R(x, y):
    if ai[1] < Ai(x, y) < ai[0]:
        return BiR[1] * 5 / r(x, y)
    return 0.
def B1m2pR(x, y):
    if aim2p[1] < Ai(x, y) < aim2p[0]:
        return BiR[1] * 5 / r(x, y)
    return 0.
def B12pR(x, y):
    if ai2p[1] < Ai(x, y) < ai2p[0]:
        return BiR[1] * 5 / r(x, y)
    return 0.
def B1totR(x, y):
    return H0(x, y) * (B1R(x, y) + B1m2pR(x, y) + B12pR(x, y))


def B2R(x, y):
    if ai[2] < Ai(x, y) < ai[1]:
        return BiR[2] * 5 / r(x, y)
    return 0.
def B2m2pR(x, y):
    if aim2p[2] < Ai(x, y) < aim2p[1]:
        return BiR[2] * 5 / r(x, y)
    return 0.
def B22pR(x, y):
    if ai2p[2] < Ai(x, y) < ai2p[1]:
        return BiR[2] * 5 / r(x, y)
    return 0.
def B2totR(x, y):
    return H0(x, y) * (B2R(x, y) + B2m2pR(x, y) + B22pR(x, y))


def B3R(x, y):
    if ai[3] < Ai(x, y) < ai[2]:
        return BiR[3] * 5 / r(x, y)
    return 0.
def B3m2pR(x, y):
    if aim2p[3] < Ai(x, y) < aim2p[2]:
        return BiR[3] * 5 / r(x, y)
    return 0.
def B32pR(x, y):
    if ai2p[3] < Ai(x, y) < ai2p[2]:
        return BiR[3] * 5 / r(x, y)
    return 0.
def B3totR(x, y):
    return H0(x, y) * (B3R(x, y) + B3m2pR(x, y) + B32pR(x, y))


def B4R(x, y):
    if ai[4] < Ai(x, y) < ai[3]:
        return BiR[4] * 5 / r(x, y)
    return 0.
def B4m2pR(x, y):
    if aim2p[4] < Ai(x, y) < aim2p[3]:
        return BiR[4] * 5 / r(x, y)
    return 0.
def B42pR(x, y):
    if ai2p[4] < Ai(x, y) < ai2p[3]:
        return BiR[4] * 5 / r(x, y)
    return 0.
def B4totR(x, y):
    return H0(x, y) * (B4R(x, y) + B4m2pR(x, y) + B42pR(x, y))


def B5R(x, y):
    if ai[5] < Ai(x, y) < ai[4]:
        return BiR[5] * 5 / r(x, y)
    return 0.
def B5m2pR(x, y):
    if aim2p[5] < Ai(x, y) < aim2p[4]:
        return BiR[5] * 5 / r(x, y)
    return 0.
def B52pR(x, y):
    if ai2p[5] < Ai(x, y) < ai2p[4]:
        return BiR[5] * 5 / r(x, y)
    return 0.
def B5totR(x, y):
    return H0(x, y) * (B5R(x, y) + B5m2pR(x, y) + B52pR(x, y))


def B6R(x, y):
    if ai[6] < Ai(x, y) < ai[5]:
        return BiR[6] * 5 / r(x, y)
    return 0.
def B6m2pR(x, y):
    if aim2p[6] < Ai(x, y) < aim2p[5]:
        return BiR[6] * 5 / r(x, y)
    return 0.
def B62pR(x, y):
    if ai2p[6] < Ai(x, y) < ai2p[5]:
        return BiR[6] * 5 / r(x, y)
    return 0.
def B6totR(x, y):
    return H0(x, y) * (B6R(x, y) + B6m2pR(x, y) + B62pR(x, y))


def B7R(x, y):
    if ai[7] < Ai(x, y) < ai[6]:
        return BiR[7] * 5 / r(x, y)
    return 0.
def B7m2pR(x, y):
    if aim2p[7] < Ai(x, y) < aim2p[6]:
        return BiR[7] * 5 / r(x, y)
    return 0.
def B72pR(x, y):
    if ai2p[7] < Ai(x, y) < ai2p[6]:
        return BiR[7] * 5 / r(x, y)
    return 0.
def B7totR(x, y):
    return H0(x, y) * (B7R(x, y) + B7m2pR(x, y) + B72pR(x, y))



def BringR(x, y, z):
    if 0. < r(x, y) < 5.:
        return BiR[8]
    return 0.

def BsptotR(x, y, z):
    return np.exp(-z ** 2 / (2 * z0diskR ** 2)) * H20(x, y) * (B0totR(x, y) + B1totR(x, y) + B2totR(x, y) + B3totR(x, y)
                                                               + B4totR(x, y) + B5totR(x, y) + B6totR(x, y) + B7totR(x, y))


def BtotR(x, y, z):
    return np.sqrt(BhaloR(x, y, z) ** 2 + (BsptotR(x, y, z) + BringR(x, y, z)) ** 2)
