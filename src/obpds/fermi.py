from numpy import exp, log

def fermi_p(psi, Vp, phi_p, nieffref, Vt):
    return nieffref*exp((-psi+Vp+phi_p)/Vt)

def fermi_n(psi, Vn, phi_n, nieffref, Vt):
    return nieffref*exp((psi+Vn-phi_n)/Vt)

def inv_fermi_p(psi, Vp, p, nieffref, Vt):
    return psi - Vp + log(p / nieffref)*Vt

def inv_fermi_n(psi, Vn, n, nieffref, Vt):
    return psi + Vn - log(n / nieffref)*Vt