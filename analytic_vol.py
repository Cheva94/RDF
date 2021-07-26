#!python3.9

'''
    Calculation: Radial Distribution Function (RDF).
    Description: using Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

from numpy import sqrt, pi, arctan

def volWedge(Rsph, a, b):
    root = sqrt(Rsph**2 - a**2 - b**2)

    T1 = Rsph**3 * (pi - 2 * arctan(a * b / (Rsph*root))) / 6
    T2 = 0.5 * (Rsph**2 * b - b**3 / 3) * (arctan(a/root) -  0.5*pi)
    T3 = 0.5 * (Rsph**2 * a - a**3 / 3) * (arctan(b/root) -  0.5*pi)
    T4 = a * b * root / 3

    return T1 + T2 + T3 + T4

def volCap(Rsph, a):
    return 0.25 * pi * (2 * Rsph**3 / 3 - Rsph**2 * a + a**3 / 3)

def volOct(Rsph, xb, yb, zb):
    if xb**2 + yb**2 + zb**2 < Rsph**2:
        return xb * yb* zb

    vol = pi * Rsph**3 / 6

    for a in [xb, yb, zb]:
        if a < Rsph:
            vol -= volCap(Rsph, a)

    for (a,b) in [(xb, yb), (xb, zb), (yb, zb)]:
        if a**2 + b**2 < Rsph**2:
            vol += volWedge(Rsph, a, b)

    return vol

def volSph(Rsph, localBox):
    vol = 0
    for xb in [localBox[0], localBox[1]]:
        for yb in [localBox[2], localBox[3]]:
            for zb in [localBox[4], localBox[5]]:
                vol += volOct(Rsph, abs(xb), abs(yb), abs(zb))

    return vol

def volShell(Rin, Rout, localBox):
    return volSph(Rout, localBox) - volSph(Rin, localBox)

def main():
    volShell(Rin, Rout, localBox)

if __name__ == "__main__":
    main()
