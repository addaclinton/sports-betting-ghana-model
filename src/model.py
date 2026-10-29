"""Compartmental model of sports betting among Ghanaian adults.

Compartments: S (susceptible), E (experimenters), R (regular bettors),
P (problem bettors), Q (recovered). Time is in years.

    dS/dt = Lam + delta*Q - beta*S*(E+R+P)/N - mu*S
    dE/dt = beta*S*(E+R+P)/N - (a1+g1+mu)*E
    dR/dt = a1*E - (a2+g2+mu)*R
    dP/dt = a2*R - (g3+mu)*P
    dQ/dt = g1*E + g2*R + g3*P - (delta+mu)*Q

Parameter classes (see README and Table 1 of the manuscript):
  A  Ghanaian national statistics (N0, mu, Lam)
  B  calibrated to a Ghanaian survey (a1, a2 from the 16:44:40 severity composition)
  C  no Ghanaian estimate (beta via R0, g1, g2, g3, delta): scenario values and ranges,
     NOT estimates.
"""
import numpy as np

# ---------- Class A: Ghanaian statistics ----------
N0 = 17_931_673            # adults aged 18+, 2021 Population and Housing Census (GSS)
LE = 65.69                 # life expectancy at birth, 2024 (World Bank)
mu = 1.0 / (LE - 18.0)     # exit rate: 1 / mean remaining adult lifetime
Lam = mu * N0              # recruitment, so that the adult population is stationary

# ---------- Class B: severity composition of 403 bettors, Volta Region (Manu et al. 2024) ----------
share = dict(E=0.16, R=0.44, P=0.40)

# Past-12-month sports betting, ages 18-44, online survey (TGM Research 2022). Used only as a rough bound.
TGM_PREV = 0.4174

# ---------- Class C: baseline scenario values (illustrative, not estimates) ----------
g1, g2, g3, delta = 0.50, 0.10, 0.10, 0.15
R0_target = 1.30


def derive(g2=g2, g3=g3, mu=mu):
    """alpha1, alpha2 such that the equilibrium composition E:R:P equals `share`."""
    c = g3 + mu
    a2 = (share['P'] / share['R']) * c
    b = a2 + g2 + mu
    a1 = (share['R'] / share['E']) * b
    return a1, a2


def R0f(beta, a1, a2, g1, g2, g3, mu):
    """Basic reproduction number (next-generation matrix method)."""
    a = a1 + g1 + mu
    b = a2 + g2 + mu
    c = g3 + mu
    return beta * (1 / a + a1 / (a * b) + a1 * a2 / (a * b * c))


a1, a2 = derive()
beta = R0_target / R0f(1.0, a1, a2, g1, g2, g3, mu)
P = dict(Lam=Lam, mu=mu, beta=beta, a1=a1, a2=a2, g1=g1, g2=g2, g3=g3, delta=delta, N0=N0)


def rhs(y, p):
    S, E, R, Pp, Q = y
    N = S + E + R + Pp + Q
    lam = p['beta'] * (E + R + Pp) / N
    dS = p['Lam'] + p['delta'] * Q - lam * S - p['mu'] * S
    dE = lam * S - (p['a1'] + p['g1'] + p['mu']) * E
    dR = p['a1'] * E - (p['a2'] + p['g2'] + p['mu']) * R
    dP = p['a2'] * R - (p['g3'] + p['mu']) * Pp
    dQ = p['g1'] * E + p['g2'] * R + p['g3'] * Pp - (p['delta'] + p['mu']) * Q
    return np.array([dS, dE, dR, dP, dQ])


def rk4(y0, p, T, h=0.005, out_every=0.25):
    """Classical fourth-order Runge-Kutta; returns (times, states) sampled every `out_every` years."""
    n = int(round(T / h))
    y = np.array(y0, float)
    ts = [0.0]
    ys = [y.copy()]
    k = int(round(out_every / h))
    for i in range(1, n + 1):
        k1 = rhs(y, p)
        k2 = rhs(y + h / 2 * k1, p)
        k3 = rhs(y + h / 2 * k2, p)
        k4 = rhs(y + h * k3, p)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i % k == 0:
            ts.append(i * h)
            ys.append(y.copy())
    return np.array(ts), np.array(ys)


# Initial state: seeding scenario (1.0% E, 0.3% R, 0.1% P of N0), horizon 30 years
E0, R0i, P0i = 0.010 * N0, 0.003 * N0, 0.001 * N0
y0 = [N0 - E0 - R0i - P0i, E0, R0i, P0i, 0.0]
T = 30.0


def endemic(p):
    """Closed-form endemic equilibrium (Theorem 5); returns (R0, [S*, E*, R*, P*, Q*])."""
    a = p['a1'] + p['g1'] + p['mu']
    b = p['a2'] + p['g2'] + p['mu']
    c = p['g3'] + p['mu']
    R0 = R0f(p['beta'], p['a1'], p['a2'], p['g1'], p['g2'], p['g3'], p['mu'])
    pp = p['a1'] * p['a2'] / (b * c)
    m = 1 + p['a1'] / b + pp
    q = (p['g1'] + p['g2'] * p['a1'] / b + p['g3'] * pp) / (p['delta'] + p['mu'])
    Ns = p['Lam'] / p['mu']
    Es = Ns * (1 - 1 / R0) / (m + q)
    Rs = p['a1'] * Es / b
    Ps = p['a2'] * Rs / c
    Qs = q * Es
    Ss = Ns / R0
    return R0, np.array([Ss, Es, Rs, Ps, Qs])
