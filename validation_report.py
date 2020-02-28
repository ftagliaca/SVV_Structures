from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

rpt_file = 'data/B737.rpt'

Bending_reg1 = genfromtxt(rpt_file,delimiter=None,skip_header=20,max_rows=5778,invalid_raise=False)

Bending_reg2 = genfromtxt(rpt_file,delimiter=None,skip_header=5816,max_rows=856,invalid_raise=False)


Jam_bent_reg1 = genfromtxt(rpt_file,delimiter=None,skip_header=6705,max_rows=5778,invalid_raise=False)

Jam_bent_reg2 = genfromtxt(rpt_file,delimiter=None,skip_header=12501,max_rows=856,invalid_raise=False)


Jam_str_reg1 = genfromtxt(rpt_file,delimiter=None,skip_header=13390,max_rows=5778,invalid_raise=False)

Jam_str_reg2 = genfromtxt(rpt_file,delimiter=None,skip_header=19186,max_rows=856,invalid_raise=False)


Bending_B737_1 = genfromtxt(rpt_file,delimiter=None,skip_header=20074,max_rows=6588,invalid_raise=False)

Bending_assembly = genfromtxt(rpt_file,delimiter=None,skip_header=26678,max_rows=16,invalid_raise=False)


Jam_bent_B737_1 = genfromtxt(rpt_file,delimiter=None,skip_header=26724,max_rows=6588,invalid_raise=False)

Jam_bent_assembly = genfromtxt(rpt_file,delimiter=None,skip_header=33328,max_rows=16,invalid_raise=False)


Jam_str_B737_1 = genfromtxt(rpt_file,delimiter=None,skip_header=33374,max_rows=6588,invalid_raise=False)

Jam_str_assembly = genfromtxt(rpt_file,delimiter=None,skip_header=39978,max_rows=16,invalid_raise=False)


Bending_RF_B737 = genfromtxt(rpt_file,delimiter=None,skip_header=40024,max_rows=6588,invalid_raise=False)

Bending_RF_assembly = genfromtxt(rpt_file,delimiter=None,skip_header=46628,max_rows=16,invalid_raise=False)


Jam_bent_RF_B737 = genfromtxt(rpt_file,delimiter=None,skip_header=46674,max_rows=6588,invalid_raise=False)

Jam_bent_RF_assembly = genfromtxt(rpt_file,delimiter=None,skip_header=53278,max_rows=16,invalid_raise=False)


Jam_str_RF_B737 = genfromtxt(rpt_file,delimiter=None,skip_header=53324,max_rows=6588,invalid_raise=False)

Jam_str_RF_assembly = genfromtxt(rpt_file,delimiter=None,skip_header=59928,max_rows=16,invalid_raise=False)


