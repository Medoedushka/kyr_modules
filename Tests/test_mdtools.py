import sys
import unittest
import random
import numpy as np

sys.path.append(r'D:\Python\MD_UNI')
from kyr_modules.md_tools import calc_num_density3D # type: ignore


def genPoints(N, box):
    points = np.array([0,0,0], ndmin=2)
    for i in range(N):
        # x = random.gauss(mu=box[0]/2, sigma=5)
        # y = random.gauss(mu=box[1]/2, sigma=5)
        # z = random.gauss(mu=box[2]/2, sigma=5)
        x = random.uniform(0, box[0] + 1)
        y = random.uniform(0, box[1] + 1)
        z = random.uniform(0, box[2] + 1) 
        points = np.append(points, [[x, y, z]], axis=0)
    return np.delete(points, 0, axis=0)


class TestMDTools(unittest.TestCase):

    def test_numdens(self):
        box = (10, 10, 27)  # Box parameters
        N = 500           # number of points
        dens_soll = N / (box[0] * box[1] * box[2])
        a = 2
        dens_f_values = []
        for s in range(2500):
            pts = genPoints(N, box)
            dens_f = calc_num_density3D(box, pts, h=3, axis=a)
            dens_avr = 1 / box[a] * np.trapz(dens_f[:,1], dens_f[:,0])
            dens_f_values.append(dens_avr)
            print(s)
        
        print("Generated density: {}\nEvaluated density: {}".format(dens_soll, np.average(dens_f_values)))
        # self.assertLessEqual(abs(np.average(dens_f_values) - dens_soll), 0.5)
        self.assertAlmostEqual(dens_avr, dens_soll, places=1)


if __name__ == "__main__":
    unittest.main()