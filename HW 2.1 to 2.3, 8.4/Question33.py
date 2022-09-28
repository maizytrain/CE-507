from string import digits
import unittest
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def PowerPlot(maxPower, resolution, labels):
    x = np.linspace(0,1,resolution)
    y = np.linspace(0,1,resolution)
    for j in range(maxPower+1):
        for i in range(len(x)):
            y[i] = x[i]**j

        plt.plot(x,y, label = "Power of " + str(j))
    if (labels == True):
        plt.legend()
    plt.show()
    return None


PowerPlot(maxPower = 10, resolution = 200, labels = True)
PowerPlot(maxPower = 1000, resolution = 200, labels = False)
