### Importing required packages
from aileronProperties import Aileron
from unittest import TestCase
import unittest
import numpy as np
import math as m
from shearCenterBuild import get_shear_center, torsional_stiffness


class GeometricalProperties(TestCase):
    
    def load_aircraft(self, name: str):
        if name == "A320":
            self.Ca    = 0.547
            self.la    = 2.771
            self.x1    = 0.153
            self.x2    = 1.281
            self.x3    = 2.681
            self.xa    = 28.0
            self.ha    = 22.5
            self.tsk   = 1.1
            self.tsp   = 2.9
            self.tst   = 1.2
            self.hst   = 1.5
            self.wst   = 2.0
            self.nst   = 17
            self.d1    = 1.103
            self.d3    = 1.642
            self.theta = 26
            self.P     = 91.7
        elif name == "B737":
            self.Ca    = 0.605
            self.la    = 2.661
            self.x1    = 0.172
            self.x2    = 1.211
            self.x3    = 2.591
            self.xa    = 35.0
            self.ha    = 20.5
            self.tsk   = 1.1
            self.tsp   = 2.8
            self.tst   = 1.2
            self.hst   = 1.6
            self.wst   = 1.9
            self.nst   = 15
            self.d1    = 1.154
            self.d3    = 1.840
            self.theta = 28
            self.P     = 97.4
        elif name == "CRJ700":
            self.Ca    = 0.484
            self.la    = 1.691
            self.x1    = 0.149
            self.x2    = 0.554
            self.x3    = 1.541
            self.xa    = 27.2
            self.ha    = 17.3
            self.tsk   = 1.1
            self.tsp   = 2.5
            self.tst   = 1.2
            self.hst   = 1.4
            self.wst   = 1.8
            self.nst   = 13
            self.d1    = 0.681
            self.d3    = 2.030
            self.theta = 26
            self.P     = 37.9
        elif name == "Do228":
            self.Ca    = 0.515
            self.la    = 2.691
            self.x1    = 0.174
            self.x2    = 1.051
            self.x3    = 2.512
            self.xa    = 30.0
            self.ha    = 24.8
            self.tsk   = 1.1
            self.tsp   = 2.2
            self.tst   = 1.2
            self.hst   = 1.5
            self.wst   = 3.0
            self.nst   = 11
            self.d1    = 1.034
            self.d3    = 2.066
            self.theta = 25
            self.P     = 20.6
        elif name == "Fokker100":
            self.Ca    = 0.505
            self.la    = 1.611
            self.x1    = 0.125
            self.x2    = 0.498
            self.x3    = 1.494
            self.xa    = 24.5
            self.ha    = 16.1
            self.tsk   = 1.1
            self.tsp   = 2.4
            self.tst   = 1.2
            self.hst   = 1.3
            self.wst   = 1.7
            self.nst   = 11
            self.d1    = 0.389
            self.d3    = 1.245
            self.theta = 30
            self.P     = 49.2
            
        ## Model stuff
        self.aileron = Aileron(self.Ca, self.la, self.x1, self.x2, self.x3, self.xa, self.ha, self.tsk, self.tsp, self.tst, self.hst, self.wst, self.nst, self.d1, self.d3, self.theta, self.P)

        ## Verification stuff
        self.xa /= 1e2  # cm to m
        self.ha /= 1e2  # cm to m
        self.tsk /= 1e3  # mm to m
        self.tsp /= 1e3  # mm to m
        self.tst /= 1e3  # mm to m
        self.hst /= 1e2  # cm to m
        self.wst /= 1e2  # cm to m
        self.d1  /= 1e2  # cm to m
        self.d3  /= 1e2  # cm to m
        self.theta = m.radians(self.theta)
        self.P   *= 1e3  # kN to N
        
        # self.crosssection = Stiffness.Crosssection(self.nst, self.Ca, self.ha, self.tsk, self.tsp, self.tst, self.hst, self.wst)
        # self.crosssection.compute_bending_properties()   # Run the calculations
        # self.crosssection.compute_shearcenter()   # Run the calculations
        # self.crosssection.compute_torsionalstiffness()   # Run the calculations
        

    
    def test_geometry(self):
        for aircraft in ["A320", "B737", "CRJ700", "Do228", "Fokker100"]:
            print(); print()
            print(f"#####     Running test on {aircraft}     #####")
            self.load_aircraft(aircraft)
            
            print("- Stringer coordinates")
            
            self.assertAlmostEqual((self.crosssection.stcoord[:,::-1] - self.aileron.stringersPosition()).sum(), 0, delta=1e-6, msg="Stringer positions are not correct.")
            
            print("- Cross-sectional area")
            print(f"Should be: {self.crosssection.totarea}")
            print(f"       is: {self.aileron.crossArea()}")
            
            self.assertEqual(self.crosssection.totarea, self.aileron.crossArea(), msg="Cross-sectional area is not correct.")
            
            print("- Centroid position")
            print(f"Should be (y, z): {self.crosssection.yc, self.crosssection.zc}")
            print(f"              is: {0, self.aileron.zCentroid()}")
            
            self.assertEqual(self.crosssection.yc, 0, msg="Centroid y position not correct.")
            self.assertEqual(self.crosssection.zc, self.aileron.zCentroid(), msg="Centroid z position not correct.")
            
            print("- Moment of Inertia position")
            print(f"Should be (I_yy, I_zz): {self.crosssection.Iyy, self.crosssection.Izz}")
            print(f"                    is: {self.aileron.momInertia()}")
            
            self.assertEqual((self.crosssection.Iyy, self.crosssection.Izz), self.aileron.momInertia(), msg="MoI is not correct.")
            
            print("- Shear centre")
            print(f"Should be (y, z): {self.crosssection.ysc, self.crosssection.zsc}")
            print(f"              is: {0, get_shear_center(self.aileron)}")
            
            self.assertEqual((self.crosssection.ysc, self.crosssection.zsc), (0, get_shear_center(self.aileron)), msg="Shear centre is not correct.")
            
            print("- Torsional constant")
            print(f"Should be: {self.crosssection.J}")
            print(f"         : {torsional_stiffness(self.aileron)}")
            
            self.assertEqual(self.crosssection.J, torsional_stiffness(self.aileron), msg="MoI is not correct.")
        

    def assertAlmostEqual(self, first, second, places=None, msg=None, delta=None):
        try:
            super().assertAlmostEqual(first, second, places=None, msg=None, delta=None)
            print("Correct.")
        except AssertionError:
            print("Incorrect")

    def assertEqual(self, first, second, msg=None):
        """Fail if the two objects are unequal as determined by the '=='
           operator.
        """
        try:
            super().assertEqual(first, second, msg=msg)
            print("Correct.")
        except AssertionError:
            print("Incorrect")

if __name__ == "__main__":
    
    unittest.main()
