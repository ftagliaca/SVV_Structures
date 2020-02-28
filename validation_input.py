from numpy import *

inp_file = 'data/B737.inp'

geom_data = genfromtxt(inp_file,delimiter=',',skip_header=9,max_rows = 6588)

surface_data = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=6598,max_rows = 6634) 

skin = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=13233, max_rows=364)
skin = append(skin,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=13597,max_rows = 1))

el_skin = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=13599, max_rows=361)
el_skin = append(el_skin,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=13960,max_rows = 1))

spar = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=13962,max_rows = 60)
spar = append(spar,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14022,max_rows = 1))

el_spar = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14024,max_rows = 53)
el_spar = append(el_spar,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14077,max_rows = 1))

RibA = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14079,max_rows = 3)
RibA = append(RibA,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14082,max_rows = 1))

RibB = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14084,max_rows = 3)
RibB = append(RibB,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14087,max_rows = 1))

RibC = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14089,max_rows = 3)
RibC = append(RibC,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14092,max_rows = 1))

RibD = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14094,max_rows = 3)
RibD = append(RibD,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14097,max_rows = 1))

LE = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14114,max_rows = 6)
LE = append(LE,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14120,max_rows = 1))

TE = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14122,max_rows = 6)
TE = append(TE,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14128,max_rows = 1))

Hinge = genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14178,max_rows = 3)
Hinge = append(Hinge,genfromtxt(inp_file,dtype=int,delimiter=',',skip_header=14181,max_rows = 1))


