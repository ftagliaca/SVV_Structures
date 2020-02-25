from aileronProperties import Aileron
from internalLoadsStresses import solveInternal
from aero_loads import AerodynamicLoad


A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
aeroLoad = AerodynamicLoad(A320, "aerodynamicloada320.dat")
q = aeroLoad.get_value_at
_ = A320.crossArea()
print(_)
_ = A320.stringersPosition()
print(_)
_ = A320.zCentroid()
print(_)
_ = A320.momInertia()
print(_)
print(solveInternal(A320,q))
