#!/usr/local/bin/python3.9

import pandas as pd
import numpy as np

def user_input():
    '''
    PREGUNTAS AL USUARIO
    '''

    # Sistema a estudiar
    xyz = pd.read_csv('example1_xyz.xsf', header = None, delim_whitespace = True, skiprows = 2, usecols = [0,1,2,3], names=['idAt', 'rx', 'ry', 'rz'])
    nAt = 103
    idAt = ['O', 'Si', 'Pt', 'H']
    Lx, Ly, Lz = 10.5599958569, 14.9399970185, 21.9999878486
    # xyz = pd.read_csv('example1_xyz.xsf', header = None, delim_whitespace = True, skiprows = 2, usecols = [0,1,2,3], names=['idAt', 'rx', 'ry', 'rz']).to_numpy()

    # Átomos a comparar - Si con O del agua (Pt)
    at1 = 'Si'
    at2 = 'Pt'

    # Filtro el archivo con los átomos seleccionados
    xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
    xyz2 = xyz[xyz['idAt'] == at2].to_numpy()

    # Frames a medir

    # Parámetros del histograma
    dr = 0.1 # bin size (incremento)
    Rcut = 10.0 # maximum radius to be considered (max Value of the histogram)

    # Output filename
    filename = 'example1_NC'
    return Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, filename

def hist_init(dr, Rcut):
    '''
    Inicialización del histograma
    '''
    nBin = int(Rcut/dr) + 1 # number of bins
    Rcut = nBin * dr # Reajusta el máximo cuando no es un múltiplo entero del incremento.
    RDF = np.zeros(nBin) # initialize array of zeros for RDF. RDF es el histograma en sí mismo. mucho muy imporante
    return nBin, Rcut, RDF

def PBC(dist, length):
    '''
    Correct a distance dist accoring to periodic boundary conditions length.
    '''

    if dist < -0.5 * length:
        dist += length
    elif dist > 0.5 * length:
        dist -= length

    return dist

def hist_up(data, dr, H):
    '''
    Updates an existing histogram.
    '''

    binIdx = int(np.sqrt(data)/dr)
    H[binIdx] += 2 # no sé bien por qué suma 2 y no 1 ! - dice que 2 es la contribución de las partículas pero en el código decía +1 . porque son dos contribuyendo (o sea i ve a j y despues j ve a i)

def sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF):
    Rcut2 = Rcut * Rcut

    rx1 = np.array(xyz1[:,1])
    ry1 = np.array(xyz1[:,2])
    rz1 = np.array(xyz1[:,3])

    rx2 = np.array(xyz2[:,1])
    ry2 = np.array(xyz2[:,2])
    rz2 = np.array(xyz2[:,3])

    nAt1 = len(rx1)
    nAt2 = len(rx2)

    # Aplicación de PBC y actualización del histograma
    # No estoy considerando duplicados si midiera el mismo átomo contra sí mismo
    # Frenkel usa las PBC de manera diferente
    for i in range(nAt1):
        for j in range(nAt2):
            dx = rx1[i] - rx2[j]
            dx = PBC(dx, Lx)

            dy = ry1[i] - ry2[j]
            dy = PBC(dy, Ly)

            dz = rz1[i] - rz2[j]
            dz = PBC(dz, Lz)

            r2 = dx * dx + dy * dy + dz * dz

            if r2 <= Rcut2: # acá tiene 0.5 L siendo L el largo de la caja # Ver si va a ser contado y contar
                hist_up(r2, dr, RDF)

    return nAt1, nAt2

def normalize(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames, RDF, filename):

    ## Cálculo de la g(r) - determine g(r)

    # Normalizamos cada bin por su volumen ya que a mayores distancias se espera un mayor conteo, incluso para el caso ideal. ¿Porque dividimos por nAt?

    nAtTot = nAt1 + nAt2
    rho = nAtTot / (Lx * Ly * Lz) # densidad de partículas (caso ideal) - densidad uniforme
    prefact = 4 * np.pi # para no estar calculandolo todas las veces. Omito el dividido 3 porque aproximo a primer orden --- considerando que la distancia es r = dr * (i + 0.5)
    dr3 = dr * dr * dr

    with open(f'{filename}.dat', 'w') as f:
        for binIdx in range(nBin):
            # volBin = prefact * ( (binIdx+1) * (binIdx+1) * (binIdx+1) - binIdx * binIdx * binIdx ) * dr3 # frenkel
            volBin = prefact * ( (binIdx + 0.5) * (binIdx + 0.5) ) * dr3 # volume between bin i and bin i+1
            nIdeal = volBin * rho
            RDF[binIdx] /= frames * nIdeal * nAtTot
            f.write(f'{RDF[binIdx]} \n')

def main():

    Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, filename = user_input()

    nBin, Rcut, RDF = hist_init(dr, Rcut)

    frames = 1 # Contador de frames

    nAt1, nAt2 = sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF)

    normalize(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames, RDF, filename)

if __name__ == "__main__":
    main()
