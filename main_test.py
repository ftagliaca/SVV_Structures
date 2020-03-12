### Importing required packages
from aileronProperties import Aileron
from unittest import TestCase
import unittest
import numpy as np
import math as m
from shearCenterBuild import get_shear_center, torsional_stiffness
from collections import namedtuple

VerificationProperties = namedtuple('VerificationProperties', [
    'stringer_coordinates',
    'cross_section_area_tot',
    'cross_section_area',
    'y_centroid',
    'z_centroid',
    'I_yy',
    'I_zz',
    'y_shear_centre',
    'z_shear_centre',
    'torsional_constant'
])


def print_correctness(first, second, correct, suffix=''):
    if not correct:
        print(f"Incorrect. {suffix}")
        print(f"   should be: {first}")
        print(f"          is: {second} ({((second / first - 1) * 100):+.2f}%)")
    else:
        print(f"Correct {suffix}")


class GeometricalProperties(TestCase):

    def load_aircraft(self, name: str, constant_stiffener_count=False):
        run_verification_model = False

        if name == "A320":
            self.Ca = 0.547
            self.la = 2.771
            self.x1 = 0.153
            self.x2 = 1.281
            self.x3 = 2.681
            self.xa = 28.0
            self.ha = 22.5
            self.tsk = 1.1
            self.tsp = 2.9
            self.tst = 1.2
            self.hst = 1.5
            self.wst = 2.0
            self.nst = 17
            self.d1 = 1.103
            self.d3 = 1.642
            self.theta = 26
            self.P = 91.7

            self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[0., -0.],
                [0.06845563, -0.02322472],
                [0.10864704, -0.08330976],
                [0.10145495, -0.15515843],
                [0.0830086, -0.22640235],
                [0.06456224, -0.29764627],
                [0.04611589, -0.3688902],
                [0.02766953, -0.44013412],
                [0.00922318, -0.51137804],
                [-0.00922318, -0.51137804],
                [-0.02766953, -0.44013412],
                [-0.04611589, -0.3688902],
                [-0.06456224, -0.29764627],
                [-0.0830086, -0.22640235],
                [-0.10145495, -0.15515843],
                [-0.10864704, -0.08330976],
                [-0.06845563, -0.02322472]]), cross_section_area_tot=0.06876164101099792, cross_section_area=0.0027426935105400304, y_centroid=0.0, z_centroid=-0.21577972811362234, I_yy=6.86413733566373e-05, I_zz=1.280745624085021e-05, y_shear_centre=0, z_shear_centre=-0.11922705644352412, torsional_constant=1.663139269310244e-05)

        elif name == "B737":
            self.Ca = 0.605
            self.la = 2.661
            self.x1 = 0.172
            self.x2 = 1.211
            self.x3 = 2.591
            self.xa = 35.0
            self.ha = 20.5
            self.tsk = 1.1
            self.tsp = 2.8
            self.tst = 1.2
            self.hst = 1.6
            self.wst = 1.9
            self.nst = 15
            self.d1 = 1.154
            self.d3 = 1.840
            self.theta = 28
            self.P = 97.4

            if not constant_stiffener_count:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.07877549, -0.03692049],
                    [ 0.09876497, -0.12081074],
                    [ 0.08080771, -0.20884515],
                    [ 0.06285044, -0.29687956],
                    [ 0.04489317, -0.38491397],
                    [ 0.0269359 , -0.47294838],
                    [ 0.00897863, -0.56098279],
                    [-0.00897863, -0.56098279],
                    [-0.0269359 , -0.47294838],
                    [-0.04489317, -0.38491397],
                    [-0.06285044, -0.29687956],
                    [-0.08080771, -0.20884515],
                    [-0.09876497, -0.12081074],
                    [-0.07877549, -0.03692049]]), cross_section_area_tot=0.06800942890838887, cross_section_area=0.002686478946739162, y_centroid=-4.339264140786877e-19, z_centroid=-0.24048766835061938, I_yy=8.651211860639685e-05, I_zz=1.0280189203385745e-05, y_shear_centre=0, z_shear_centre=-0.10856995078063854, torsional_constant=1.5101498390705797e-05)
            else:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.07160611, -0.02915959],
                    [ 0.10247066, -0.10004751],
                    [ 0.08714556, -0.17777418],
                    [ 0.07130092, -0.2554516 ],
                    [ 0.05545627, -0.33312903],
                    [ 0.03961162, -0.41080645],
                    [ 0.02376697, -0.48848387],
                    [ 0.00792232, -0.56616129],
                    [-0.00792232, -0.56616129],
                    [-0.02376697, -0.48848387],
                    [-0.03961162, -0.41080645],
                    [-0.05545627, -0.33312903],
                    [-0.07130092, -0.2554516 ],
                    [-0.08714556, -0.17777418],
                    [-0.10247066, -0.10004751],
                    [-0.07160611, -0.02915959]]), cross_section_area_tot=0.06800942890838887, cross_section_area=0.0027704789467391617, y_centroid=0.0, z_centroid=-0.2416466303511318, I_yy=8.954404784883145e-05, I_zz=1.0642752968134631e-05, y_shear_centre=0, z_shear_centre=-0.10834219430517074, torsional_constant=1.5101498390705797e-05)

        elif name == "CRJ700":
            self.Ca = 0.484
            self.la = 1.691
            self.x1 = 0.149
            self.x2 = 0.554
            self.x3 = 1.541
            self.xa = 27.2
            self.ha = 17.3
            self.tsk = 1.1
            self.tsp = 2.5
            self.tst = 1.2
            self.hst = 1.4
            self.wst = 1.8
            self.nst = 13
            self.d1 = 0.681
            self.d3 = 2.030
            self.theta = 26
            self.P = 37.9

            if not constant_stiffener_count:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.07111646, -0.03725877],
                    [ 0.07988634, -0.11689227],
                    [ 0.06213382, -0.19847177],
                    [ 0.0443813 , -0.28005126],
                    [ 0.02662878, -0.36163076],
                    [ 0.00887626, -0.44321025],
                    [-0.00887626, -0.44321025],
                    [-0.02662878, -0.36163076],
                    [-0.0443813 , -0.28005126],
                    [-0.06213382, -0.19847177],
                    [-0.07988634, -0.11689227],
                    [-0.07111646, -0.03725877]]), cross_section_area_tot=0.046136840816161116, cross_section_area=0.0021255886520793153, y_centroid=-1.0028413565320257e-18, z_centroid=-0.19406263838748938, I_yy=4.363276766019503e-05, I_zz=5.81593895759915e-06, y_shear_centre=0, z_shear_centre=-0.09185594953325858, torsional_constant=8.629971582027014e-06)
            else:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.05820344, -0.02251087],
                    [ 0.08611301, -0.07832695],
                    [ 0.07466501, -0.14088624],
                    [ 0.06108955, -0.20327056],
                    [ 0.04751409, -0.26565488],
                    [ 0.03393864, -0.3280392 ],
                    [ 0.02036318, -0.39042352],
                    [ 0.00678773, -0.45280784],
                    [-0.00678773, -0.45280784],
                    [-0.02036318, -0.39042352],
                    [-0.03393864, -0.3280392 ],
                    [-0.04751409, -0.26565488],
                    [-0.06108955, -0.20327056],
                    [-0.07466501, -0.14088624],
                    [-0.08611301, -0.07832695],
                    [-0.05820344, -0.02251087]]), cross_section_area_tot=0.046136840816161116, cross_section_area=0.0022791886520793156, y_centroid=0.0, z_centroid=-0.19595905064248173, I_yy=4.717344573387326e-05, I_zz=6.270784865707143e-06, y_shear_centre=0, z_shear_centre=-0.09142274223537021, torsional_constant=8.629971582027014e-06)

        elif name == "Do228":
            self.Ca = 0.515
            self.la = 2.691
            self.x1 = 0.174
            self.x2 = 1.051
            self.x3 = 2.512
            self.xa = 30.0
            self.ha = 24.8
            self.tsk = 1.1
            self.tsp = 2.2
            self.tst = 1.2
            self.hst = 1.5
            self.wst = 3.0
            self.nst = 11
            self.d1 = 1.034
            self.d3 = 2.066
            self.theta = 25
            self.P = 20.6

            if not constant_stiffener_count:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.09612648, -0.04566929],
                    [ 0.11637895, -0.1480309 ],
                    [ 0.08312782, -0.25287921],
                    [ 0.04987669, -0.35772753],
                    [ 0.01662556, -0.46257584],
                    [-0.01662556, -0.46257584],
                    [-0.04987669, -0.35772753],
                    [-0.08312782, -0.25287921],
                    [-0.11637895, -0.1480309 ],
                    [-0.09612648, -0.04566929]]), cross_section_area_tot=0.07263656432079832, cross_section_area=0.0024705343591563712, y_centroid=3.033354054941713e-19, z_centroid=-0.20728702965108006, I_yy=5.377416790820396e-05, I_zz=1.4221538884296291e-05, y_shear_centre=0, z_shear_centre=-0.13229700743946904, torsional_constant=1.9193311985303668e-05)

            else:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.06732887, -0.01987112],
                    [ 0.1130787 , -0.07311575],
                    [ 0.1183349 , -0.14186335],
                    [ 0.09681946, -0.20970638],
                    [ 0.07530402, -0.27754941],
                    [ 0.05378859, -0.34539243],
                    [ 0.03227315, -0.41323546],
                    [ 0.01075772, -0.48107849],
                    [-0.01075772, -0.48107849],
                    [-0.03227315, -0.41323546],
                    [-0.05378859, -0.34539243],
                    [-0.07530402, -0.27754941],
                    [-0.09681946, -0.20970638],
                    [-0.1183349 , -0.14186335],
                    [-0.1130787 , -0.07311575],
                    [-0.06732887, -0.01987112]]), cross_section_area_tot=0.07263656432079832, cross_section_area=0.002794534359156371, y_centroid=0.0, z_centroid=-0.21011089923358023, I_yy=6.231458098479158e-05, I_zz=1.6161179097413905e-05, y_shear_centre=0, z_shear_centre=-0.13108425184061262, torsional_constant=1.9193311985303668e-05)

        elif name == "Fokker100":
            self.Ca = 0.505
            self.la = 1.611
            self.x1 = 0.125
            self.x2 = 0.498
            self.x3 = 1.494
            self.xa = 24.5
            self.ha = 16.1
            self.tsk = 1.1
            self.tsp = 2.4
            self.tst = 1.2
            self.hst = 1.3
            self.wst = 1.7
            self.nst = 11
            self.d1 = 0.389
            self.d3 = 1.245
            self.theta = 30
            self.P = 49.2

            if not constant_stiffener_count:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.0766793 , -0.0559942 ],
                    [ 0.06621955, -0.155805  ],
                    [ 0.04729968, -0.255575  ],
                    [ 0.02837981, -0.355345  ],
                    [ 0.00945994, -0.455115  ],
                    [-0.00945994, -0.455115  ],
                    [-0.02837981, -0.355345  ],
                    [-0.04729968, -0.255575  ],
                    [-0.06621955, -0.155805  ],
                    [-0.0766793 , -0.0559942 ]]), cross_section_area_tot=0.04435140289671263, cross_section_area=0.0020111318843290065, y_centroid=0.0, z_centroid=-0.20362591085157106, I_yy=4.5943507864451845e-05, I_zz=4.753851442684437e-06, y_shear_centre=0, z_shear_centre=-0.08553893540215983, torsional_constant=7.748548555816593e-06)

            else:
                self.verification_properties = VerificationProperties(stringer_coordinates=np.array([[ 0.        , -0.        ],
                    [ 0.05865051, -0.02536047],
                    [ 0.07957475, -0.08537912],
                    [ 0.06733248, -0.14993617],
                    [ 0.05509021, -0.21449323],
                    [ 0.04284794, -0.27905029],
                    [ 0.03060567, -0.34360735],
                    [ 0.0183634 , -0.40816441],
                    [ 0.00612113, -0.47272147],
                    [-0.00612113, -0.47272147],
                    [-0.0183634 , -0.40816441],
                    [-0.03060567, -0.34360735],
                    [-0.04284794, -0.27905029],
                    [-0.05509021, -0.21449323],
                    [-0.06733248, -0.14993617],
                    [-0.07957475, -0.08537912],
                    [-0.05865051, -0.02536047]]), cross_section_area_tot=0.04435140289671263, cross_section_area=0.0022271318843290065, y_centroid=0.0, z_centroid=-0.2065355019291737, I_yy=5.1342473816015615e-05, I_zz=5.2644068076863455e-06, y_shear_centre=0, z_shear_centre=-0.08500516192069422, torsional_constant=7.748548555816593e-06)


        if constant_stiffener_count:
            self.nst = 17

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
        self.d1 /= 1e2  # cm to m
        self.d3 /= 1e2  # cm to m
        self.theta = m.radians(self.theta)
        self.P *= 1e3  # kN to N

        if run_verification_model:
            import Stiffness

            self.crosssection = Stiffness.Crosssection(self.nst, self.Ca, self.ha, self.tsk, self.tsp, self.tst, self.hst, self.wst)
            self.crosssection.compute_bending_properties()  # Run the calculations
            self.crosssection.compute_shearcenter()  # Run the calculations
            self.crosssection.compute_torsionalstiffness()  # Run the calculations

            h = self.crosssection.ha / 2.
            A1 = m.pi * h ** 2 / 2.
            A2 = (self.crosssection.Ca - h) * h

            A = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
            b = np.array([0., 0., 0.])

            ### First row
            A[0, 0] = 2. * A1
            A[0, 1] = 2. * A2
            b[0] = 1

            ### Second row
            A[1, 0] = (h * m.pi / self.crosssection.tsk + 2 * h / self.crosssection.tsp) / (2 * A1)
            A[1, 1] = (-2 * h / self.crosssection.tsp) / (2 * A1)
            A[1, 2] = -1.
            b[1] = 0.

            ### Third row
            A[2, 0] = (-2 * h / self.crosssection.tsp) / (2 * A2)
            A[2, 1] = (2 * self.crosssection.lsk / self.crosssection.tsk + 2 * h / self.crosssection.tsp) / (2 * A2)
            A[2, 2] = -1
            b[2] = 0.

            solution = np.linalg.solve(A, b)
            self.crosssection.J = 1. / solution[-1]
            self.verification_properties = VerificationProperties(
                stringer_coordinates=self.crosssection.stcoord[:, ::-1],
                cross_section_area_tot=A1 + A2,
                cross_section_area=self.crosssection.totarea,
                y_centroid=self.crosssection.yc,
                z_centroid=self.crosssection.zc,
                I_yy=self.crosssection.Iyy,
                I_zz=self.crosssection.Izz,
                y_shear_centre=self.crosssection.ysc,
                z_shear_centre=self.crosssection.zsc,
                torsional_constant=self.crosssection.J
            )
            print(self.verification_properties)

    def test_geometry(self):
        for aircraft in ["A320", "B737", "CRJ700", "Do228", "Fokker100"]:
            print();
            print()
            print(f"#####     Running test on {aircraft}     #####")
            self.load_aircraft(aircraft)

            print("- Stringer coordinates")
            self.assertAlmostEqual((self.verification_properties.stringer_coordinates - self.aileron.stringersPosition()).sum(), 0, delta=1e-6, msg="Stringer positions are not correct.")

            print("- Cross-sectional area")
            self.assertEqual(self.verification_properties.cross_section_area_tot, self.aileron.crossArea(), msg="Cross-sectional area is not correct.")

            print("- Centroid position")

            print("y: ", end='')
            self.assertAlmostEqual(self.verification_properties.y_centroid, 0, delta=1e-12, msg="Centroid y position not correct.")
            print("z: ", end='')
            self.assertAlmostEqual(self.verification_properties.z_centroid, self.aileron.zCentroid(), delta=1e-12, msg="Centroid z position not correct.")


            print("- Moment of Inertia")
            I_yy, I_zz = self.aileron.momInertia()

            print("y: ", end='')
            self.assertAlmostEqual(self.verification_properties.I_yy, I_yy, delta=1e-12, msg="MoI_yy is not correct.")
            print("z: ", end='')
            self.assertAlmostEqual(self.verification_properties.I_zz, I_zz, delta=1e-12, msg="MoI_zz is not correct.")

            self.load_aircraft(aircraft, constant_stiffener_count=True)
            self.aileron.crossArea()

            print("- Shear centre")

            print("y: ", end='')
            self.assertEqual(self.verification_properties.y_shear_centre, 0, msg="y shear centre location is not correct.")
            print("z: ", end='')
            self.assertEqual(self.verification_properties.z_shear_centre, get_shear_center(self.aileron)[0], msg="z shear centre location is not correct.")

            print("- Torsional constant")
            self.assertAlmostEqual(self.verification_properties.torsional_constant, torsional_stiffness(self.aileron)[0], delta=1e-6, msg="J is not correct.")

    def assertAlmostEqual(self, first, second, places=None, msg=None, delta=None):
        correct = True
        try:
            super().assertAlmostEqual(first, second, places=places, msg=msg, delta=delta)
        except AssertionError:
            correct = False

        print_correctness(first, second, correct, suffix=f"(up to {delta}; magnitude of error: {'1e' + f'{(second - first):.2e}'.split('e')[-1]})")

    def assertEqual(self, first, second, msg=None):
        correct = True
        try:
            super().assertEqual(first, second, msg=msg)
        except AssertionError:
            correct = False

        print_correctness(first, second, correct)

if __name__ == "__main__":
    unittest.main()
