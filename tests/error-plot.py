import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

def plot_erros(h, error, error0, figname):
  x = np.log(h)
  y1 = np.log(error)
  y2 = np.log(error0)

  #fitting data into straight line:
  m1, b1 = np.polyfit(x, y1, 1)
  m2, b2 = np.polyfit(x, y2, 1) #error0

  xf = np.linspace(np.min(x), np.max(x), 200)
  y1f =  m1*xf + b1
  y2f =  m2*xf + b2  #error0

  # seting-up Figures parameters:
  rcParams.update({'figure.autolayout': True})
  rcParams["figure.figsize"] = [4.5, 5.5]
  rcParams['text.usetex'] = True

  plt.figure()
  plt.plot(x, y1, 'bo')
  plt.plot(xf, y1f, '--b', label='$||error||_E$')
  plt.plot(x, y2, 'ro')
  plt.plot(xf, y2f, '--r', label='$||error||_0$')
  plt.xlabel('log h')
  plt.ylabel('log error')
  plt.legend()
  plt.savefig(figname)
  #plt.show()

hq = np.array([0.5, 0.25, 0.125, 0.0625, 0.0625/2.])
ht = np.array([1.0, 0.5, 0.25, 0.125, 0.0625])

## Triangle element:
#error order1:
error = np.array([7.84384, 5.66874, 3.43309, 1.80281, 0.914325])
error0 = np.array([0.97281, 0.492347, 0.177818, 0.0491276, 0.0126559])
plot_erros(ht, error, error0, "triangle-order1.eps")

#error order 2:
error = np.array([5.91471, 2.77374, 0.777577, 0.204498, 0.0519473])
error0 = np.array([0.444053, 0.129974, 0.0187519, 0.00247243, 0.000314414])
plot_erros(ht, error, error0, "triangle-order2.eps")

## Quad element:
# error order 1:
error = np.array([7.72681, 4.29193, 2.13392, 1.06749, 0.533907 ])
error0 = np.array([0.872926, 0.278562, 0.0699183, 0.017506, 0.00437865])
plot_erros(hq, error, error0, "quad-order1.eps")

#error order 2:
error = np.array([4.29717, 2.33426, 1.92219, 1.81116, 2.86259])
error0 = np.array([0.310539, 0.138998, 0.0904695, 0.0712541, 0.344019])
plot_erros(hq, error, error0, "quad-order2.eps")

