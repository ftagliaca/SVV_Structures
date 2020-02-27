from numpy import *
import matplotlib.pyplot as plt
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


ax = plt.axes(projection='3d')



def Validation_plot(LC, value): #LC = loading case, value is type of stress to show
    index = LC[:,0].astype(int)-1

    if max(index > 6587):
        i = min(index)

        for s in surface_data[index,:]-1:
            surface = array([geom_data[s[1],1:4],geom_data[s[2],1:4],geom_data[s[3],1:4],geom_data[s[4],1:4]])
            X = surface[:,0]
            Y = surface[:,1]
            Z = surface[:,2]
            r = abs(LC[i,value]/max(LC[:,value]))
            g = 1-r
            c = [r,g,0.5]
            ax.plot_trisurf(X, Y, Z,color=c)
            i += 1
        plt.show()
        
    else:
        
        r = abs(LC[:,value]/(max(LC[:,value])+0.1))
        g = 1-r
        color = array([r,g,0.5*ones(len(r))]).T
        print(color)
        print(shape(color))
        X = geom_data[index,1]
        Y = geom_data[index,2]
        Z = geom_data[index,3]
        print(shape(X))
        ax.scatter(X,Y,Z,',',c=color)

        plt.show()

