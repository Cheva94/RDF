from pandas import read_csv
from numpy import zeros, array

################################################################################
######################## Input processing functions
################################################################################

def userfile_mono(input_file, atom):
    '''
    Process the xsf input file given by the user in the monocomponent case.
    '''

    xsf = read_csv(input_file, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    frames_total = int(xsf.iloc[0,1])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    nAtTot = int(xsf.iloc[7,0])
    nAt = xsf.iloc[8:(nAtTot+8),0].value_counts()[atom]
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)

    return frames_total, Lx, Ly, Lz, nAtTot, nAt, xyz_all

################################################################################
######################## Histogramming functions
################################################################################

def hist_init_hdf(minimum, maximum, increment):
    '''
    Initialize 2D histogram for HDF.
    '''

    nBin = int((maximum - minimum)/increment)# + 1
    maximum = nBin * increment + minimum
    H = zeros(nBin)

    return nBin, maximum, H

def hist_up_hdf(minimum, data, increment, H):
    '''
    Updates the existing 2D histogram for HDF.
    '''

    binIdx = int((data - minimum) / increment)
    H[binIdx] += 1

def hist_init_pdf(Xmax, Ymax, increment):
    '''
    Initialize 3D histogram for PDF.
    '''

    nBinX = int(Xmax/increment) + 1
    Xmax = nBinX * increment

    nBinY = int(Ymax/increment) + 1
    Ymax = nBinY * increment

    H = zeros((nBinX, nBinY))

    return nBinX, nBinY, Xmax, Ymax, H

def hist_up_pdf(dataX, dataY, increment, H):
    '''
    Updates the existing 3D histogram for PDF.
    '''

    binIdxX = int(dataX/increment)
    binIdxY = int(dataY/increment)

    H[binIdxX, binIdxY] += 1

################################################################################
######################## Sampling functions
################################################################################

def sample_hdf(xyz, dh, HDF, nAt, Hmin):
    '''
    Determines HDF.
    '''

    rz = array(xyz[:,3])

    for i in range(nAt):
        hist_up_hdf(Hmin, rz[i], dh, HDF)

def sample_pdf(Lx, Ly, xyz, dxy, PDF, nAt):
    '''
    Determines PDF.
    '''

    rx = array(xyz[:,1])
    ry = array(xyz[:,2])

    for i in range(nAt):
        hist_up_pdf(rx[i], ry[i], dxy, PDF)

################################################################################
######################## Normalizing functions
################################################################################

def normalize_hdf(dh, nBin, frames_count, HDF, output_file, Hmin, Hmax):
    '''
    Normalize the HDF.
    '''

    HDF /= frames_count

    with open(f'{output_file}.csv', 'w') as f:
        f.write(f'{Hmin:.2f}, 0.0000\n')
        for binIdx in range(nBin):
            h = (binIdx + 0.5) * dh + Hmin
            f.write(f'{h:.2f}, {HDF[binIdx]:.4f}\n')
        f.write(f'{Hmax:.2f}, 0.0000\n')

def normalize_pdf(dxy, nBinX, nBinY, frames_count, PDF, output_file):
    '''
    Normalize the PDF.
    '''

    PDF /= frames_count

    with open(f'{output_file}.csv', 'w') as f:
        for binIdxX in range(nBinX):
            x = (binIdxX + 0.5) * dxy
            for binIdxY in range(nBinY):
                y = (binIdxY + 0.5) * dxy

                f.write(f'{x:.2f}, {y:.2f}, {PDF[binIdxX, binIdxY]:.4f}\n')

                # if PDF[binIdxX, binIdxY] == 0:
                #     continue
                # else:
                #     f.write(f'{x:.2f}, {y:.2f}, {PDF[binIdxX, binIdxY]:.4f}\n')