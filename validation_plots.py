from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from mpl_toolkits import mplot3d
from validation_report import *
from validation_input import *


##REPORT DATA:
##    1. Bending_reg1 (surf)
##    2. Bending_reg2 (nodes)
##    3. Jam_bent_reg1 (surf)
##    4. Jam_bent_reg2 (nodes)
##    5. Jam_str_reg1 (surf)
##    6. Jam_str_reg2 (nodes)
##    7. Bending_B737_1 (nodes)
##    8. Bending_assembly (nodes)
##    9. Jam_bent_B737_1 (nodes)
##    10. Jam_bent_assembly (nodes)*
##    11. Jam_str_B737_1 (nodes)*
##    12. Jam_str_assembly (nodes)
##    13. Bending_RF_B737 (nodes)
##    14. Bending_RF_assembly (nodes)
##    15. Jam_bent_RF_B737 (nodes)
##    16. Jam_bent_RF_assembly (nodes)
##    17. Jam_str_RF_B737 (nodes)
##    18. Jam_str_RF_assembly (nodes)
##
##DATA FORMATS:
##col 0 = label
##    
##for RPT.DATA 1, 2, 3, 4, 5, 6,
##col 1 = Int Pt. col 2 = S.Mises @ 1     col 3 = S.Mises @ 2     col 4 = S.S12 @ 1   col 5 = S.S12 @ 2
##
##for RPT.DATA 7, 8, 9, 10, 11, 12
##col 1 = U.Magn @ 1  col 1 = U.U1 @ 1  col 3 = U.U2 @ 1  col 4 = U.U3 @ 1
##
##for RPT.DATA 13, 14, 15, 16, 17, 18
##col 1 = RF.Magn @ 1  col 1 = RF.RF1 @ 1  col 3 = RF.RF2 @ 1  col 4 = RF.RF3 @ 1


def LC_meaning(loading_case):
    if loading_case is Bending_reg1:
        return "Bending_reg1", 1
    elif loading_case is Bending_reg2:
        return "Bending_reg2", 2
    elif loading_case is Jam_bent_reg1:
        return "Jam_bent_reg1", 3
    elif loading_case is Jam_bent_reg2:
        return "Jam_bent_reg2", 4
    elif loading_case is Jam_str_reg1:
        return "Jam_str_reg1", 5
    elif loading_case is Jam_str_reg2:
        return "Jam_str_reg2", 6
    elif loading_case is Bending_B737_1:
        return "Bending_B737_1", 7
    elif loading_case is Bending_assembly:
        return "Bending_assembly", 8
    elif loading_case is Jam_bent_B737_1:
        return "Jam_bent_B737_1", 9
    elif loading_case is Jam_bent_assembly:
        return "Jam_bent_assembly", 10
    elif loading_case is Jam_str_B737_1:
        return "Jam_str_B737_1", 11
    elif loading_case is Jam_str_assembly:
        return "Jam_str_assembly", 12
    elif loading_case is Bending_RF_B737:
        return "Bending_RF_B737", 13
    elif loading_case is Bending_RF_assembly:
        return "Bending_RF_assembly", 14
    elif loading_case is Jam_bent_RF_B737:
        return "Jam_bent_RF_B737", 15
    elif loading_case is Jam_bent_RF_assembly:
        return "Jam_bent_RF_assembly", 16
    elif loading_case is Jam_str_RF_B737:
        return "Jam_str_RF_B737", 17
    elif loading_case is Jam_str_RF_assembly:
        return "Jam_str_RF_assembly", 18

loading_cases_sorted = [
    "Bending_reg1",
    "Bending_reg2",
    "Jam_bent_reg1",
    "Jam_bent_reg2",
    "Jam_str_reg1",
    "Jam_str_reg2",
    "Bending_B737_1",
    "Bending_assembly",
    "Jam_bent_B737_1",
    "Jam_bent_assembly",
    "Jam_str_B737_1",
    "Jam_str_assembly",
    "Bending_RF_B737",
    "Bending_RF_assembly",
    "Jam_bent_RF_B737",
    "Jam_bent_RF_assembly",
    "Jam_str_RF_B737",
    "Jam_str_RF_assembly"
]

def column_meaning_and_units(loading_case: int):
    result = []
    result.append(("Index", "[-]"))
    if 1 <= loading_case <= 6:
        result.append(("Int Pt.",     "?"))     # col. 1
        result.append(("S.Mises @ 1", "GPa"))   # col. 2
        result.append(("S.Mises @ 2", "GPa"))   # col. 3
        result.append(("S.S12 @ 1",   "GPa"))   # col. 4
        result.append(("S.S12 @ 2",   "GPa"))   # col. 5
    elif 7 <= loading_case <= 12:
        result.append(("U.Magn @ 1",  "[mm]"))  # col. 1
        result.append(("U.U1 @ 1",    "[mm]"))    # col. 2
        result.append(("U.U2 @ 1",    "[mm]"))    # col. 3
        result.append(("U.U3 @ 1",    "[mm]"))    # col. 4
    elif 13 <= loading_case <= 18:
        result.append(("RF.Magn @ 1", "?"))  # col. 1
        result.append(("RF.RF1 @ 1",  "?"))   # col. 2
        result.append(("RF.RF2 @ 1",  "?"))   # col. 3
        result.append(("RF.RF3 @ 1",  "?"))   # col. 4
    
    return result


def Validation_plot(LC, value): #LC = loading case, value is type of stress to show
    global surface_data
    
    loading_case_name, loading_case_index = LC_meaning(LC)
    column_meaning, column_unit = column_meaning_and_units(loading_case_index)[value]

    fig = plt.figure(f"{loading_case_name} {column_meaning}")
    ax = plt.axes(projection='3d')
    ax.set_title(f"{loading_case_name} {column_meaning}")

    indices = LC[:,0].astype(int)-1

    color_dimension = LC[:,value]

    norm = matplotlib.colors.Normalize(color_dimension.min(), color_dimension.max())
    m = plt.cm.ScalarMappable(norm=norm, cmap='jet')

    print(color_dimension.min(), color_dimension.max())

    fcolors = m.to_rgba(color_dimension)
    fcolors[:, -1] = 1

    fig.colorbar(m, ax=ax, label=column_unit)
    

    if indices.max() > 6587:
        i = indices.min()

        for s in surface_data[indices,:]-1:
            surface = array([geom_data[s[1],1:4],geom_data[s[2],1:4],geom_data[s[3],1:4],geom_data[s[4],1:4]])
            X = surface[:,0]
            Y = surface[:,1]
            Z = surface[:,2]
            ax.plot_trisurf(X, Y, Z, color=fcolors[i])
            i += 1
        
    else:
        
        #cb1.set_label('Some Units')

        X = geom_data[indices,1]
        Y = geom_data[indices,2]
        Z = geom_data[indices,3]
        print(shape(X))
        ax.scatter(X,Y,Z,',',c=fcolors)

    plt.show()

Validation_plot(Bending_B737_1, 3)