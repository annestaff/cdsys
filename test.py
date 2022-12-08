from ddeint import ddeint
import matplotlib.pyplot as plt
from numpy import linspace, array


def model(Y, t, d):
    x, y = Y(t)
    xd, yd = Y(t - d)
    return array([0.5 * x * (1 - yd), -0.5 * y * (1 - xd)])


g = lambda t: array([1, 2])
tt = linspace(2, 30, 20000)

for d in [0, 0.2]:
    yy = ddeint(model, g, tt, fargs=(d,))
    # WE PLOT X AGAINST Y
    plt.plot(yy[:, 0], yy[:, 1], lw=2, label='delay = %.01f' % d)

plt.legend()  # display the legend
plt.show()
