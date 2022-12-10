from math import exp

import matplotlib.pyplot as plt
import numpy as np
from ddeint import ddeint

V0 = 1e4
T0 = 7e7


def model(Y, t, pv=210, delta_v=5, kl=0.4, beta=5e-7, gt=0.8, beta_prim=3e-8, rho=2.6, fi=0.33, delta_i=2,
          kn=2.5, ke=5e-5, pf=1e-5, delta_f=2, beta_cn=1, hc=1e4, tau_c=6, pc=1.2, delta_e=0.57, beta_bn=0.03, hb=1e4,
          tau_b=4, delta_p=0.5, pb=0.52, ps=12, delta_s=2, pl=4, delta_l=0.015, ks=0.8):

    viral_load, target_cells, infected_cells, ifn, virus_resistance, naive_CD8, effector, naive_B, plasma_cells, \
    short_antibodies, long_antibodies = Y(t)

    VIRAL_LOAD = pv * infected_cells - delta_v * viral_load - ks * viral_load * short_antibodies - kl * viral_load * \
    long_antibodies - beta * viral_load * target_cells

    TARGET_CELLS = gt * (target_cells + virus_resistance) * (1 - (target_cells + virus_resistance + infected_cells)/T0)\
    - beta_prim * viral_load * target_cells + rho * virus_resistance - fi * ifn * target_cells

    INFECTED_CELLS = beta_prim * viral_load * target_cells - delta_i * infected_cells - kn * infected_cells * ifn - ke \
    * infected_cells * effector

    IFN = pf * infected_cells - delta_f * ifn

    VIRUS_RESISTANCE = fi * ifn * target_cells - rho * virus_resistance

    NAIVE_CD8 = - beta_cn * (viral_load / (viral_load + hc)) * naive_CD8

    EFFECTOR = beta_cn * (viral_load*(t - tau_c) / (viral_load*(t - tau_c) + hc)) * naive_CD8*(t - tau_c) * \
    exp(pc * tau_c) - delta_e * effector

    NAIVE_B = - beta_bn * (viral_load / (viral_load + hb)) * naive_B

    PLASMA_CELLS = beta_bn * (viral_load*(t - tau_b) / (viral_load*(t - tau_b) + hb)) * naive_B*(t - tau_b) * exp(
        pb * tau_b) - delta_p * plasma_cells

    SHORT_ANTIBODIES = ps * plasma_cells - delta_s * short_antibodies

    LONG_ANTIBODIES = pl * plasma_cells - delta_l * long_antibodies

    return np.array([VIRAL_LOAD, TARGET_CELLS, INFECTED_CELLS, IFN, VIRUS_RESISTANCE, NAIVE_CD8, EFFECTOR, NAIVE_B, PLASMA_CELLS, SHORT_ANTIBODIES, LONG_ANTIBODIES])


legend = ['V', 'T', 'I', 'F', 'R', 'Cn', 'E', 'Bn', 'P', 'As', 'Al']


def run_model(y0=np.array([V0, T0, 0, 0, 0, 100, 0, 100, 0, 0, 0]), t0=0, interval=6, graph_step=10000):
    tt = np.linspace(t0, interval, graph_step)

    y = ddeint(model, lambda t: y0, tt)
    return {"y0": y0, "tt": tt, "y": y}


def plot_model(y, y0=np.array([V0, T0, 0, 0, 0, 100, 0, 100, 0, 0, 0]), tt=np.linspace(0, 6, 10000)):
    for u in range(len(y0)):
        plt.plot(tt, y[:, u], label=legend[u])
        plt.legend()
        plt.show()


my_model = run_model()
plot_model(my_model["y"])
print("Done")
