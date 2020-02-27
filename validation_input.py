from numpy import *



geom_data = genfromtxt('B737.inp',delimiter=',',skip_header=9,max_rows = 6588)

surface_data = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=6598,max_rows = 6634) 

skin = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=13233, max_rows=364)
skin = append(skin,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=13597,max_rows = 1))

el_skin = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=13599, max_rows=361)
el_skin = append(el_skin,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=13960,max_rows = 1))

spar = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=13962,max_rows = 60)
spar = append(spar,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14022,max_rows = 1))

el_spar = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14024,max_rows = 53)
el_spar = append(el_spar,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14077,max_rows = 1))

RibA = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14079,max_rows = 3)
RibA = append(RibA,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14082,max_rows = 1))

RibB = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14084,max_rows = 3)
RibB = append(RibB,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14087,max_rows = 1))

RibC = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14089,max_rows = 3)
RibC = append(RibC,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14092,max_rows = 1))

RibD = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14094,max_rows = 3)
RibD = append(RibD,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14097,max_rows = 1))

LE = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14114,max_rows = 6)
LE = append(LE,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14120,max_rows = 1))

TE = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14122,max_rows = 6)
TE = append(TE,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14128,max_rows = 1))

Hinge = genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14178,max_rows = 3)
Hinge = append(Hinge,genfromtxt('B737.inp',dtype=int,delimiter=',',skip_header=14181,max_rows = 1))


