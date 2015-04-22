#
#   Copyright (c) 2015, Scott J Maddox
#
#   This file is part of Open Band Parameters Device Simulator (OBPDS).
#
#   OBPDS is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as published
#   by the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   OBPDS is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with OBPDS.  If not, see <http://www.gnu.org/licenses/>.
#
#############################################################################
'''
Fermi-Dirac integrals
    
[1] T. Fukushima, "Precise and fast computation of Fermi-Dirac integral
    of integer and half integer order by piecewise minimax rational
    approximation," Applied Mathematics and Computation, vol. 259,
    pp. 708-729, May 2015.
'''
import numpy
from numpy import exp, sqrt, log
from scipy.optimize import newton


__all__ = ['fdm1h', 'fd1h', 'dfd1h', 'ifd1h']


def _fdm1h(x):
    '''
    Double-precision rational minimax approximation of Fermi-Dirac
    integral of order -1/2.
    
    [1] T. Fukushima, "Precise and fast computation of Fermi-Dirac integral
        of integer and half integer order by piecewise minimax rational
        approximation," Applied Mathematics and Computation, vol. 259,
        pp. 708-729, May 2015.
    '''
    if(x < -2e0):
        ex=exp(x)
        t=ex*7.38905609893065023e0
        fd=ex*(1.77245385090551603e0 \
        -ex*(40641.4537510284430e0 \
        +t*(9395.7080940846442e0 \
        +t*(649.96168315267301e0 \
        +t*(12.7972295804758967e0 \
        +t*0.00153864350767585460e0 \
        ))))/(32427.1884765292940e0 \
        +t*(11079.9205661274782e0 \
        +t*(1322.96627001478859e0 \
        +t*(63.738361029333467e0 \
        +t)))))
    elif(x < 0e0):
        s=-0.5e0*x
        t=1e0-s
        fd=(272.770092131932696e0 \
        +t*(30.8845653844682850e0 \
        +t*(-6.43537632380366113e0 \
        +t*(14.8747473098217879e0 \
        +t*(4.86928862842142635e0 \
        +t*(-1.53265834550673654e0 \
        +t*(-1.02698898315597491e0 \
        +t*(-0.177686820928605932e0 \
        -t*0.00377141325509246441e0 \
        ))))))))/(293.075378187667857e0 \
        +s*(305.818162686270816e0 \
        +s*(299.962395449297620e0 \
        +s*(207.640834087494249e0 \
        +s*(92.0384803181851755e0 \
        +s*(37.0164914112791209e0 \
        +s*(7.88500950271420583e0 \
        +s)))))))
    elif(x < 2e0):
        t=0.5e0*x
        fd=(3531.50360568243046e0 \
        +t*(6077.5339658420037e0 \
        +t*(6199.7700433981326e0 \
        +t*(4412.78701919567594e0 \
        +t*(2252.27343092810898e0 \
        +t*(811.84098649224085e0 \
        +t*(191.836401053637121e0 \
        +t*23.2881838959183802e0 \
        )))))))/(3293.83702584796268e0 \
        +t*(1528.97474029789098e0 \
        +t*(2568.48562814986046e0 \
        +t*(925.64264653555825e0 \
        +t*(574.23248354035988e0 \
        +t*(132.803859320667262e0 \
        +t*(29.8447166552102115e0 \
        +t)))))))
    elif(x < 5e0):
        t=0.3333333333333333333e0*(x-2e0)
        fd=(4060.70753404118265e0 \
        +t*(10812.7291333052766e0 \
        +t*(13897.5649482242583e0 \
        +t*(10628.4749852740029e0 \
        +t*(5107.70670190679021e0 \
        +t*(1540.84330126003381e0 \
        +t*(284.452720112970331e0 \
        +t*29.5214417358484151e0 \
        )))))))/(1564.58195612633534e0 \
        +t*(2825.75172277850406e0 \
        +t*(3189.16066169981562e0 \
        +t*(1955.03979069032571e0 \
        +t*(828.000333691814748e0 \
        +t*(181.498111089518376e0 \
        +t*(32.0352857794803750e0 \
        +t)))))))
    elif(x < 10e0):
        t=0.2e0*x-1e0
        fd=(1198.41719029557508e0 \
        +t*(3263.51454554908654e0 \
        +t*(3874.97588471376487e0 \
        +t*(2623.13060317199813e0 \
        +t*(1100.41355637121217e0 \
        +t*(267.469532490503605e0 \
        +t*(25.4207671812718340e0 \
        +t*0.389887754234555773e0 \
        )))))))/(273.407957792556998e0 \
        +t*(595.918318952058643e0 \
        +t*(605.202452261660849e0 \
        +t*(343.183302735619981e0 \
        +t*(122.187622015695729e0 \
        +t*(20.9016359079855933e0 \
        +t))))))
    elif(x < 20e0):
        t=0.1e0*x-1e0
        fd=(9446.00169435237637e0 \
        +t*(36843.4448474028632e0 \
        +t*(63710.1115419926191e0 \
        +t*(62985.2197361074768e0 \
        +t*(37634.5231395700921e0 \
        +t*(12810.9898627807754e0 \
        +t*(1981.56896138920963e0 \
        +t*81.4930171897667580e0 \
        )))))))/(1500.04697810133666e0 \
        +t*(5086.91381052794059e0 \
        +t*(7730.01593747621895e0 \
        +t*(6640.83376239360596e0 \
        +t*(3338.99590300826393e0 \
        +t*(860.499043886802984e0 \
        +t*(78.8565824186926692e0 \
        +t)))))))
    elif(x < 40e0):
        t=0.05e0*x-1e0
        fd=(22977.9657855367223e0 \
        +t*(123416.616813887781e0 \
        +t*(261153.765172355107e0 \
        +t*(274618.894514095795e0 \
        +t*(149710.718389924860e0 \
        +t*(40129.3371700184546e0 \
        +t*(4470.46495881415076e0 \
        +t*132.684346831002976e0 \
        )))))))/(2571.68842525335676e0 \
        +t*(12521.4982290775358e0 \
        +t*(23268.1574325055341e0 \
        +t*(20477.2320119758141e0 \
        +t*(8726.52577962268114e0 \
        +t*(1647.42896896769909e0 \
        +t*(106.475275142076623e0 \
        +t)))))))
    else:
        w=1e0/(x*x)
        t=1600e0*w
        fd=sqrt(x)*2.*(1e0 \
        -w*(0.411233516712009968e0 \
        +t*(0.00110980410034088951e0 \
        +t*(0.0000113689298990173683e0 \
        +t*(2.56931790679436797e-7 \
        +t*(9.97897786755446178e-9 \
        +t*8.67667698791108582e-10))))))
    return fd
fdm1h = numpy.vectorize(_fdm1h)

def _fd1h(x):
    '''
    Double-precision rational minimax approximation of Fermi-Dirac
    integral of order 1/2.
    
    [1] T. Fukushima, "Precise and fast computation of Fermi-Dirac integral
        of integer and half integer order by piecewise minimax rational
        approximation," Applied Mathematics and Computation, vol. 259,
        pp. 708-729, May 2015.
    '''
    if(x < -2.e0):
        ex=exp(x)
        t=ex*7.38905609893065023e0
        fd=ex*(0.886226925452758014e0 \
        -ex*(19894.4553386951666e0 \
        +t*(4509.64329955948557e0 \
        +t*(303.461789035142376e0 \
        +t*(5.7574879114754736e0 \
        +t*0.00275088986849762610e0 \
        ))))/(63493.915041308052e0 \
        +t*(19070.1178243603945e0 \
        +t*(1962.19362141235102e0 \
        +t*(79.250704958640158e0 \
        +t)))))
    elif(x < 0.e0):
        s=-0.5e0*x
        t=1.e0-s
        fd=(149.462587768865243e0 \
        +t*(22.8125889885050154e0 \
        +t*(-0.629256395534285422e0 \
        +t*(9.08120441515995244e0 \
        +t*(3.35357478401835299e0 \
        +t*(-0.473677696915555805e0 \
        +t*(-0.467190913556185953e0 \
        +t*(-0.0880610317272330793e0 \
        -t*0.00262208080491572673e0 \
        ))))))))/(269.94660938022644e0 \
        +s*(343.6419926336247e0 \
        +s*(323.9049470901941e0 \
        +s*(218.89170769294024e0 \
        +s*(102.31331350098315e0 \
        +s*(36.319337289702664e0 \
        +s*(8.3317401231389461e0 \
        +s)))))))
    elif(x < 2.e0):
        t=0.5e0*x
        fd=(71652.717119215557e0 \
        +t*(134954.734070223743e0 \
        +t*(153693.833350315645e0 \
        +t*(123247.280745703400e0 \
        +t*(72886.293647930726e0 \
        +t*(32081.2499422362952e0 \
        +t*(10210.9967337762918e0 \
        +t*(2152.71110381320778e0 \
        +t*232.906588165205042e0 \
        ))))))))/(105667.839854298798e0 \
        +t*(31946.0752989314444e0 \
        +t*(71158.788776422211e0 \
        +t*(15650.8990138187414e0 \
        +t*(13521.8033657783433e0 \
        +t*(1646.98258283527892e0 \
        +t*(618.90691969249409e0 \
        +t*(-3.36319591755394735e0 \
        +t))))))))
    elif(x < 5.e0):
        t=0.3333333333333333333e0*(x-2.e0)
        fd=(23744.8706993314289e0 \
        +t*(68257.8589855623002e0 \
        +t*(89327.4467683334597e0 \
        +t*(62766.3415600442563e0 \
        +t*(20093.6622609901994e0 \
        +t*(-2213.89084119777949e0 \
        +t*(-3901.66057267577389e0 \
        -t*948.642895944858861e0 \
        )))))))/(9488.61972919565851e0 \
        +t*(12514.8125526953073e0 \
        +t*(9903.44088207450946e0 \
        +t*(2138.15420910334305e0 \
        +t*(-528.394863730838233e0 \
        +t*(-661.033633995449691e0 \
        +t*(-51.4481470250962337e0 \
        +t)))))))
    elif(x < 10.e0):
        t=0.2e0*x-1.e0
        fd=(311337.452661582536e0 \
        +t*(1.11267074416648198e6 \
        +t*(1.75638628895671735e6 \
        +t*(1.59630855803772449e6 \
        +t*(910818.935456183774e0 \
        +t*(326492.733550701245e0 \
        +t*(65507.2624972852908e0 \
        +t*4809.45649527286889e0 \
        )))))))/(39721.6641625089685e0 \
        +t*(86424.7529107662431e0 \
        +t*(88163.7255252151780e0 \
        +t*(50615.7363511157353e0 \
        +t*(17334.9774805008209e0 \
        +t*(2712.13170809042550e0 \
        +t*(82.2205828354629102e0 \
        -t)))))))*0.999999999999999877e0
    elif(x < 20.e0):
        t=0.1e0*x-1.e0
        fd=(7.26870063003059784e6 \
        +t*(2.79049734854776025e7 \
        +t*(4.42791767759742390e7 \
        +t*(3.63735017512363365e7 \
        +t*(1.55766342463679795e7 \
        +t*(2.97469357085299505e6 \
        +t*154516.447031598403e0 \
        ))))))/(340542.544360209743e0 \
        +t*(805021.468647620047e0 \
        +t*(759088.235455002605e0 \
        +t*(304686.671371640343e0 \
        +t*(39289.4061400542309e0 \
        +t*(582.426138126398363e0 \
        +t*(11.2728194581586028e0 \
        -t)))))))
    elif(x < 40.e0):
        t=0.05e0*x-1.e0
        fd=(4.81449797541963104e6 \
        +t*(1.85162850713127602e7 \
        +t*(2.77630967522574435e7 \
        +t*(2.03275937688070624e7 \
        +t*(7.41578871589369361e6 \
        +t*(1.21193113596189034e6 \
        +t*63211.9545144644852e0 \
        ))))))/(80492.7765975237449e0 \
        +t*(189328.678152654840e0 \
        +t*(151155.890651482570e0 \
        +t*(48146.3242253837259e0 \
        +t*(5407.08878394180588e0 \
        +t*(112.195044410775577e0 \
        -t))))))
    else:
        w=1.e0/(x*x)
        s=1.e0-1600.e0*w
        fd=x*sqrt(x)*0.666666666666666667e0*(1.e0+w \
        *(8109.79390744477921e0 \
        +s*(342.069867454704106e0 \
        +s*1.07141702293504595e0)) \
        /(6569.98472532829094e0 \
        +s*(280.706465851683809e0 \
        +s)))
    return fd
fd1h = numpy.vectorize(_fd1h)


def _dfd1h(x):
    '''
    Double-precision rational minimax approximation of the derivative of the
    Fermi-Dirac integral of order 1/2.
    
    [1] T. Fukushima, "Precise and fast computation of Fermi-Dirac integral
        of integer and half integer order by piecewise minimax rational
        approximation," Applied Mathematics and Computation, vol. 259,
        pp. 708-729, May 2015.
    '''
    return _fdm1h(x)/2.
dfd1h = numpy.vectorize(_dfd1h)


def _ifd1h(nu):
    '''
    Inverse Fermi-Dirac integral of order 1/2.

    Parameters
    ----------
    nu : float
        normalized carrier concentration, n/Nc.
    
    Returns
    -------
    eta : float
        normalized Fermi energy, (phi_n-Ec)/Vt
    '''
    f = lambda eta: _fd1h(eta) - nu
    fprime = lambda eta: _dfd1h(eta)
    if nu < 10:
        guess = log(nu)
    else:
        guess = nu**1.5
    return newton(f, guess, fprime=fprime)
ifd1h = numpy.vectorize(_ifd1h)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    _, ax = plt.subplots()
    phi = numpy.linspace(-30, 50, 10000)
    ax.semilogy(phi, fd1h(phi))
    ax.semilogy(phi, dfd1h(phi))
#     plt.show()
    _, ax = plt.subplots()
    nu = numpy.logspace(-10, 3, 10000)
    ax.semilogx(nu, ifd1h(nu))
    plt.show()