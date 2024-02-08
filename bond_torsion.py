#!/usr/bin/python3.10

import argparse
from pandas import read_csv
from numpy import sqrt, inner, pi, mean, std, arctan2, cross, array, argsort

def PBC(dist, length):
    '''
    Correct a distance using Periodic Boundary Conditions (PBC).
    '''

    return dist - length * int(2*dist/length)

def main():

    name = args.input_file
    center1 = args.central_atom1
    center2 = args.central_atom2
    neigh1 = args.neighboring_atom1
    neigh2 = args.neighboring_atom2
    Rcut = args.Rcut
    Rcut2 = Rcut * Rcut
    convfact = 180 / pi
    Rmin = args.Rmin
    if Rmin == None:
        Rmin = 0.1
    Rmin2 = Rmin * Rmin

    xsf = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    xyz = xsf.iloc[8:,0:4].reset_index(drop=True)

    ID, BTA = [],[]
    aux1, aux2, aux3 = [], [], []

    # Mismos átomos centrales
    if center1 == center2:        
        if neigh1 == neigh2 == center1:
            nAt = xyz.iloc[:,0].value_counts()[center1]
            xyz = xyz[xyz['idAt'] == center1]
            rx = xyz.iloc[:, 1]
            ry = xyz.iloc[:, 2]
            rz = xyz.iloc[:, 3]

            for i in range(nAt):
                for j in range(i+1, nAt):
                    dx = rx.iloc[i] - rx.iloc[j]
                    dy = ry.iloc[i] - ry.iloc[j]
                    dz = rz.iloc[i] - rz.iloc[j]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux1.append((rx.index[i]+1,rx.index[j]+1,sqrt(d2), r2))

            L = len(aux1)
            for m in range(L):
                for n in range(m+1, L):
                    for k in range(L):
                        if (aux1[m][0] == aux1[n][0]) and (aux1[m][1] == aux1[k][1]) and (aux1[m][0] != aux1[k][0]):
                            b1 = aux1[n][3]
                            b2 = aux1[m][3]
                            b3 = -array(aux1[k][3])
                            angle = convfact * arctan2(aux1[m][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux1[n][1]}-{center1}{aux1[m][0]}-{center2}{aux1[m][1]}-{neigh2}{aux1[k][0]}-')
                            BTA.append(angle)

                        elif (aux1[m][0] == aux1[n][0]) and (aux1[n][1] == aux1[k][1]) and (aux1[n][0] != aux1[k][0]):
                            b1 = aux1[m][3]
                            b2 = aux1[n][3]
                            b3 = -array(aux1[k][3])
                            angle = convfact * arctan2(aux1[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux1[m][1]}-{center1}{aux1[n][0]}-{center2}{aux1[n][1]}-{neigh2}{aux1[k][0]}-')
                            BTA.append(angle)

        elif neigh1 == neigh2:
            nCenter = xyz.iloc[:,0].value_counts()[center1]
            nNeigh = xyz.iloc[:,0].value_counts()[neigh1]

            xyzCenter = xyz[xyz['idAt'] == center1]
            xyzNeigh = xyz[xyz['idAt'] == neigh1]

            rxCenter = xyzCenter.iloc[:, 1]
            ryCenter = xyzCenter.iloc[:, 2]
            rzCenter = xyzCenter.iloc[:, 3]

            rxNeigh = xyzNeigh.iloc[:, 1]
            ryNeigh = xyzNeigh.iloc[:, 2]
            rzNeigh = xyzNeigh.iloc[:, 3]

            for i in range(nCenter):
                for j in range(nNeigh):
                    dx = rxNeigh.iloc[j] - rxCenter.iloc[i]
                    dy = ryNeigh.iloc[j] - ryCenter.iloc[i]
                    dz = rzNeigh.iloc[j] - rzCenter.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux1.append((rxCenter.index[i]+1,rxNeigh.index[j]+1,sqrt(d2), r2))
            
            for i in range(nCenter):
                for j in range(i+1, nCenter):
                    dx = rxCenter.iloc[j] - rxCenter.iloc[i]
                    dy = ryCenter.iloc[j] - ryCenter.iloc[i]
                    dz = rzCenter.iloc[j] - rzCenter.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux2.append((rxCenter.index[i]+1,rxCenter.index[j]+1,sqrt(d2), r2))

            L1 = len(aux1)
            L2 = len(aux2)
            for m in range(L2):
                for n in range(L1):
                    for k in range(L1):
                        if (aux2[m][0] == aux1[n][0]) and (aux2[m][1] == aux1[k][0]):
                            b1 = aux1[n][3]
                            b2 = aux2[m][3]
                            b3 = aux1[k][3]
                            angle = convfact * arctan2(aux2[m][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux1[n][1]}-{center1}{aux2[m][0]}-{center2}{aux2[m][1]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)

        elif neigh1 == center1:
            nCenter = xyz.iloc[:,0].value_counts()[center1]
            nNeigh = xyz.iloc[:,0].value_counts()[neigh2]

            xyzCenter = xyz[xyz['idAt'] == center1]
            xyzNeigh = xyz[xyz['idAt'] == neigh2]

            rxCenter = xyzCenter.iloc[:, 1]
            ryCenter = xyzCenter.iloc[:, 2]
            rzCenter = xyzCenter.iloc[:, 3]

            rxNeigh = xyzNeigh.iloc[:, 1]
            ryNeigh = xyzNeigh.iloc[:, 2]
            rzNeigh = xyzNeigh.iloc[:, 3]

            for i in range(nCenter):
                for j in range(nNeigh):
                    dx = rxNeigh.iloc[j] - rxCenter.iloc[i]
                    dy = ryNeigh.iloc[j] - ryCenter.iloc[i]
                    dz = rzNeigh.iloc[j] - rzCenter.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux1.append((rxCenter.index[i]+1,rxNeigh.index[j]+1,sqrt(d2), r2))
            
            for i in range(nCenter):
                for j in range(i+1, nCenter):
                    dx = rxCenter.iloc[j] - rxCenter.iloc[i]
                    dy = ryCenter.iloc[j] - ryCenter.iloc[i]
                    dz = rzCenter.iloc[j] - rzCenter.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux2.append((rxCenter.index[i]+1,rxCenter.index[j]+1,sqrt(d2), r2))

            L1 = len(aux1)
            L2 = len(aux2)
            for m in range(L2):
                for n in range(m+1, L2):
                    for k in range(L1):
                        # Enlace m como enlace central
                        if (aux2[m][0] == aux2[n][0]) and (aux2[m][1] == aux1[k][0]):
                            print(1)
                            #print(f'-{neigh1}{aux2[n][1]}-{center1}{aux2[m][0]}-{center2}{aux2[m][1]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[n][3]
                            b2 = aux2[m][3]
                            b3 = aux1[k][3]
                            angle = convfact * arctan2(aux2[m][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[n][1]}-{center1}{aux2[m][0]}-{center2}{aux2[m][1]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)
                        
                        elif (aux2[m][1] == aux2[n][1]) and (aux2[m][0] == aux1[k][0]):
                            print(2)
                            #print(f'-{neigh1}{aux2[n][0]}-{center1}{aux2[m][1]}-{center2}{aux2[m][0]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[n][3]
                            b2 = aux2[m][3]
                            b3 = -array(aux1[k][3])
                            angle = convfact * arctan2(aux2[m][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[n][0]}-{center1}{aux2[m][1]}-{center2}{aux2[m][0]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)

                        elif (aux2[m][0] == aux2[n][1]) and (aux2[m][1] == aux1[k][0]):
                            print(3)
                            #print(f'-{neigh1}{aux2[n][0]}-{center1}{aux2[m][1]}-{center2}{aux2[m][1]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[n][3]
                            b2 = aux2[m][3]
                            b3 = -array(aux1[k][3])
                            angle = convfact * arctan2(aux2[m][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[n][0]}-{center1}{aux2[m][1]}-{center2}{aux2[m][1]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)
                        
                        elif (aux2[m][1] == aux2[n][0]) and (aux2[m][0] == aux1[k][0]):
                            print(4)
                            #print(f'-{neigh1}{aux2[n][1]}-{center1}{aux2[m][1]}-{center2}{aux2[m][0]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[n][3]
                            b2 = aux2[m][3]
                            b3 = aux1[k][3]
                            angle = convfact * arctan2(aux2[m][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[n][1]}-{center1}{aux2[m][1]}-{center2}{aux2[m][0]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)
                        
                        # Enlace n como enlace central
                        elif (aux2[m][0] == aux2[n][0]) and (aux2[n][1] == aux1[k][0]):
                            print(5)
                            #print(f'-{neigh1}{aux2[m][1]}-{center1}{aux2[n][0]}-{center2}{aux2[n][1]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[m][3]
                            b2 = aux2[n][3]
                            b3 = aux1[k][3]
                            angle = convfact * arctan2(aux2[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[m][1]}-{center1}{aux2[n][0]}-{center2}{aux2[n][1]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)
                        
                        elif (aux2[m][1] == aux2[n][1]) and (aux2[n][0] == aux1[k][0]):
                            print(6)
                            #print(f'-{neigh1}{aux2[m][0]}-{center1}{aux2[n][1]}-{center2}{aux2[n][0]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[m][3]
                            b2 = aux2[n][3]
                            b3 = -array(aux1[k][3])
                            angle = convfact * arctan2(aux2[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[m][0]}-{center1}{aux2[n][1]}-{center2}{aux2[n][0]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)

                        elif (aux2[m][0] == aux2[n][1]) and (aux2[n][0] == aux1[k][0]):
                            print(7)
                            #print(f'-{neigh1}{aux2[m][1]}-{center1}{aux2[n][1]}-{center2}{aux2[n][0]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[m][3]
                            b2 = aux2[n][3]
                            b3 = -array(aux1[k][3])
                            angle = convfact * arctan2(aux2[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[m][1]}-{center1}{aux2[n][1]}-{center2}{aux2[n][0]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)

                        elif (aux2[m][1] == aux2[n][0]) and (aux2[n][1] == aux1[k][0]):
                            print(8)
                            #print(f'-{neigh1}{aux2[m][0]}-{center1}{aux2[n][0]}-{center2}{aux2[n][1]}-{neigh2}{aux1[k][1]}-')
                            b1 = aux2[m][3]
                            b2 = aux2[n][3]
                            b3 = -array(aux1[k][3])
                            angle = convfact * arctan2(aux2[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                            angle = 180 - abs(angle)
                            ID.append(f'-{neigh1}{aux2[m][0]}-{center1}{aux2[n][0]}-{center2}{aux2[n][1]}-{neigh2}{aux1[k][1]}-')
                            BTA.append(angle)
                    
        elif neigh2 == center1:
            print('Invertir extremos')
            exit()

        else:
            nCenter = xyz.iloc[:,0].value_counts()[center1]
            xyzCenter = xyz[xyz['idAt'] == center1]
            rxCenter = xyzCenter.iloc[:, 1]
            ryCenter = xyzCenter.iloc[:, 2]
            rzCenter = xyzCenter.iloc[:, 3]

            nNeigh1 = xyz.iloc[:,0].value_counts()[neigh1]
            xyzNeigh1 = xyz[xyz['idAt'] == neigh1]
            rxNeigh1 = xyzNeigh1.iloc[:, 1]
            ryNeigh1 = xyzNeigh1.iloc[:, 2]
            rzNeigh1 = xyzNeigh1.iloc[:, 3]

            nNeigh2 = xyz.iloc[:,0].value_counts()[neigh2]
            xyzNeigh2 = xyz[xyz['idAt'] == neigh2]
            rxNeigh2 = xyzNeigh2.iloc[:, 1]
            ryNeigh2 = xyzNeigh2.iloc[:, 2]
            rzNeigh2 = xyzNeigh2.iloc[:, 3]

            for i in range(nCenter):
                for j in range(nNeigh1):
                    dx = rxNeigh1.iloc[j] - rxCenter.iloc[i]
                    dy = ryNeigh1.iloc[j] - ryCenter.iloc[i]
                    dz = rzNeigh1.iloc[j] - rzCenter.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux1.append((rxCenter.index[i]+1,rxNeigh1.index[j]+1,sqrt(d2), r2))

            for i in range(nCenter):
                for j in range(nNeigh2):
                    dx = rxNeigh2.iloc[j] - rxCenter.iloc[i]
                    dy = ryNeigh2.iloc[j] - ryCenter.iloc[i]
                    dz = rzNeigh2.iloc[j] - rzCenter.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux2.append((rxCenter.index[i]+1,rxNeigh2.index[j]+1,sqrt(d2), r2))

            for i in range(nCenter):
                for j in range(nCenter):
                    dx = rxCenter.iloc[j] - rxCenter.iloc[i]
                    dy = ryCenter.iloc[j] - ryCenter.iloc[i]
                    dz = rzCenter.iloc[j] - rzCenter.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux3.append((rxCenter.index[i]+1,rxCenter.index[j]+1,sqrt(d2), r2))

            L1 = len(aux1)
            L2 = len(aux2)
            L3 = len(aux3)
            for m in range(L1):
                for n in range(L3):
                    if aux1[m][0] == aux3[n][0]:
                        for k in range(L2):
                            if aux3[n][1] == aux2[k][0]:
                                b1 = aux1[m][3]
                                b2 = aux3[n][3]
                                b3 = aux2[k][3]
                                
                                angle = convfact * arctan2(aux3[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                                angle = 180 - abs(angle)
                                ID.append(f'-{neigh1}{aux1[m][1]}-{center1}{aux3[n][0]}-{center2}{aux3[n][1]}-{neigh2}{aux2[k][1]}-')
                                BTA.append(angle)

    # Distintos átomos centrales
    else:
        nCenter1 = xyz.iloc[:,0].value_counts()[center1]
        xyzCenter1 = xyz[xyz['idAt'] == center1]
        rxCenter1 = xyzCenter1.iloc[:, 1]
        ryCenter1 = xyzCenter1.iloc[:, 2]
        rzCenter1 = xyzCenter1.iloc[:, 3]

        nCenter2 = xyz.iloc[:,0].value_counts()[center2]
        xyzCenter2 = xyz[xyz['idAt'] == center2]
        rxCenter2 = xyzCenter2.iloc[:, 1]
        ryCenter2 = xyzCenter2.iloc[:, 2]
        rzCenter2 = xyzCenter2.iloc[:, 3]

        if neigh1 == neigh2:
            nNeigh = xyz.iloc[:,0].value_counts()[neigh1]
            xyzNeigh = xyz[xyz['idAt'] == neigh1]
            rxNeigh = xyzNeigh.iloc[:, 1]
            ryNeigh = xyzNeigh.iloc[:, 2]
            rzNeigh = xyzNeigh.iloc[:, 3]

            for i in range(nCenter1):
                for j in range(nNeigh):
                    dx = rxNeigh.iloc[j] - rxCenter1.iloc[i]
                    dy = ryNeigh.iloc[j] - ryCenter1.iloc[i]
                    dz = rzNeigh.iloc[j] - rzCenter1.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux1.append((rxCenter1.index[i]+1,rxNeigh.index[j]+1,sqrt(d2), r2))

            for i in range(nCenter2):
                for j in range(nNeigh):
                    dx = rxNeigh.iloc[j] - rxCenter2.iloc[i]
                    dy = ryNeigh.iloc[j] - ryCenter2.iloc[i]
                    dz = rzNeigh.iloc[j] - rzCenter2.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux2.append((rxCenter2.index[i]+1,rxNeigh.index[j]+1,sqrt(d2), r2))

            for i in range(nCenter1):
                for j in range(nCenter2):
                    dx = rxCenter2.iloc[j] - rxCenter1.iloc[i]
                    dy = ryCenter2.iloc[j] - ryCenter1.iloc[i]
                    dz = rzCenter2.iloc[j] - rzCenter1.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux3.append((rxCenter1.index[i]+1,rxCenter2.index[j]+1,sqrt(d2), r2))

            L1 = len(aux1)
            L2 = len(aux2)
            L3 = len(aux3)
            for m in range(L1):
                for n in range(L3):
                    if aux1[m][0] == aux3[n][0]:
                        for k in range(L2):
                            if aux3[n][1] == aux2[k][0]:
                                b1 = aux1[m][3]
                                b2 = aux3[n][3]
                                b3 = aux2[k][3]

                                angle = convfact * arctan2(aux3[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                                angle = 180 - abs(angle)
                                ID.append(f'-{neigh1}{aux1[m][1]}-{center1}{aux3[n][0]}-{center2}{aux3[n][1]}-{neigh2}{aux2[k][1]}-')
                                BTA.append(angle)
        else:
            nNeigh1 = xyz.iloc[:,0].value_counts()[neigh1]
            xyzNeigh1 = xyz[xyz['idAt'] == neigh1]
            rxNeigh1 = xyzNeigh1.iloc[:, 1]
            ryNeigh1 = xyzNeigh1.iloc[:, 2]
            rzNeigh1 = xyzNeigh1.iloc[:, 3]

            nNeigh2 = xyz.iloc[:,0].value_counts()[neigh2]
            xyzNeigh2 = xyz[xyz['idAt'] == neigh2]
            rxNeigh2 = xyzNeigh2.iloc[:, 1]
            ryNeigh2 = xyzNeigh2.iloc[:, 2]
            rzNeigh2 = xyzNeigh2.iloc[:, 3]

            for i in range(nCenter1):
                for j in range(nNeigh1):
                    dx = rxNeigh1.iloc[j] - rxCenter1.iloc[i]
                    dy = ryNeigh1.iloc[j] - ryCenter1.iloc[i]
                    dz = rzNeigh1.iloc[j] - rzCenter1.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux1.append((rxCenter1.index[i]+1,rxNeigh1.index[j]+1,sqrt(d2), r2))

            for i in range(nCenter2):
                for j in range(nNeigh2):
                    dx = rxNeigh2.iloc[j] - rxCenter2.iloc[i]
                    dy = ryNeigh2.iloc[j] - ryCenter2.iloc[i]
                    dz = rzNeigh2.iloc[j] - rzCenter2.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux2.append((rxCenter2.index[i]+1,rxNeigh2.index[j]+1,sqrt(d2), r2))

            for i in range(nCenter1):
                for j in range(nCenter2):
                    dx = rxCenter2.iloc[j] - rxCenter1.iloc[i]
                    dy = ryCenter2.iloc[j] - ryCenter1.iloc[i]
                    dz = rzCenter2.iloc[j] - rzCenter1.iloc[i]

                    dx = PBC(dx, Lx)
                    dy = PBC(dy, Ly)
                    dz = PBC(dz, Lz)
                    r2 = [dx, dy, dz]

                    d2 = dx**2 + dy**2 + dz**2
                    if ((d2>= Rmin2) and (d2<= Rcut2)):
                        aux3.append((rxCenter1.index[i]+1,rxCenter2.index[j]+1,sqrt(d2), r2))

            L1 = len(aux1)
            L2 = len(aux2)
            L3 = len(aux3)
            for m in range(L1):
                for n in range(L3):
                    if aux1[m][0] == aux3[n][0]:
                        for k in range(L2):
                            if aux3[n][1] == aux2[k][0]:
                                b1 = aux1[m][3]
                                b2 = aux3[n][3]
                                b3 = aux2[k][3]
                                
                                angle = convfact * arctan2(aux3[n][2]*inner(b1,cross(b2, b3)), inner(cross(b1, b2), cross(b2, b3)))
                                angle = 180 - abs(angle)
                                ID.append(f'-{neigh1}{aux1[m][1]}-{center1}{aux3[n][0]}-{center2}{aux3[n][1]}-{neigh2}{aux2[k][1]}-')
                                BTA.append(angle)

    Summary = f'\nSummary\n  -{neigh1}-{center1}-{center2}-{neigh2}- = ({mean(BTA):.3f} +- {std(BTA):.3f})°\n  Count = {len(BTA)} <<'

    ID = array(ID)
    BTA = array(BTA)
    AS = argsort(BTA)
    ID = ID[AS]
    BTA = BTA[AS]

    with open(f'BTA_{neigh1}-{center1}-{center2}-{neigh2}.csv', 'w') as f:
        f.write('==== Bond torsion angle in degrees ==== \n\n')
        f.write('Atoms ID\tAngle \n')
        for i in range(len(BTA)):
            f.write(f'{ID[i]}\t{BTA[i]:.3f}\n')
        f.write(f'{Summary}-')

    print(f'{Summary}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered as bond length.")

    parser.add_argument('neighboring_atom1')
    
    parser.add_argument('central_atom1')

    parser.add_argument('central_atom2')
    
    parser.add_argument('neighboring_atom2')

    parser.add_argument('--Rmin', type = float,
                        help = "Minimum distance to be considered as bond length.")

    args = parser.parse_args()

    main()
