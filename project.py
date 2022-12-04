from math import exp

import matplotlib.pyplot as plt
import scipy.integrate as spi
import numpy as np


# y = [V, T, I, F, R, Cn, E, Bn, P, As, Al]

V0 = 1e4
T0 = 7e7


def f(y, t, pv=210, delta_v=5, kl=0.4, beta=5e-7, gt=0.8, T0=10e7, beta_prim=3e-8, rho=2.6, fi=0.33, delta_i=2, kn=2.5,
      ke=5e-5, pf=1e-5, delta_f=2, beta_cn=1, hc=1e4, teta_c=6, pc=1.2, delta_e=0.57, beta_bn=0.03, hb=1e4, teta_b=4,
      delta_p=0.5, pb=0.52, ps=12, delta_s=2, pl=4, delta_l=0.015):
    z = np.zeros(np.size(y))
    z[0] = pv * y[2] - delta_v * y[0] - kl * y[9] - beta * y[0] * y[1]
    z[1] = gt * (y[1] + y[4]) * (1 - ((y[1] + y[4] + y[2]) / T0)) - beta_prim * y[0] * y[1] + rho * y[4] - fi * y[3] * \
           y[1]
    z[2] = beta_prim * y[0] * y[1] - delta_i * y[2] - kn * y[2] * y[3] - ke * y[2] * y[6]
    z[3] = pf * y[2] - delta_f * y[3]
    z[4] = fi * y[3] * y[1] - rho * y[4]
    z[5] = - beta_cn * (y[0] / (y[0] + hc)) * y[5]
    z[6] = beta_cn * (y[0] * (t - teta_c) / (y[0] * (t - teta_c) + hc)) * y[5] * (t - teta_c) * exp(
        pc * teta_c) - delta_e * y[6]
    z[7] = - beta_bn * (y[0] / (y[0] + hb)) * y[7]
    z[8] = beta_bn * (y[0] * (t - teta_b) / (y[0] * (t - teta_b) + hb)) * y[7] * (t - teta_b) * exp(
        pb * teta_b) - delta_p * y[8]
    z[9] = ps * y[8] - delta_s * y[9]
    z[10] = pl * y[8] - delta_l * y[10]
    return z


def afficher():
    y0 = np.array([V0, T0, 0, 0, 0, 100, 0, 100, 0, 0, 0])
    t0 = 0
    T = 6
    pas = 0.01
    tps = np.arange(t0, T, pas)

    y = spi.odeint(f, y0, tps)

    for u in range(len(y0)):
        plt.plot(tps, y[:, u])
    plt.show()


afficher()
