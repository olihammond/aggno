#Dr. Oliver S. Hammond, September 2024
#Estimate aggregation number from forward intensity
#eg. $ python3 aggno.py -f BBY0070798.dat -svd 0.844 -svld 4.922 -sd 1.057 -smw 946.652 -sld 0.77 -sc 0.008
#input args: filename, solvent density and sld, solute density, molecular weight and sld, solute conc in g/ml
#output printed and saved as '.ag0' csv file

import argparse
import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Estimate aggregation number from previously background-subtracted scattering data')

parser.add_argument('-f', '--file', dest='filename', default='None', type=str, help='Data filename')
parser.add_argument('-sd', '--soldens', dest='sol_dens', default=0.87, type=float, help='Solute density (g cm-3)')
parser.add_argument('-svd', '--solvdens', dest='solv_dens', default=0.844, type=float, help='Solvent density (g cm-3)')
parser.add_argument('-sld', '--solsld', dest='sol_sld', default=-0.08, type=float, help='Solute SLD (10^-6 Å-2 == 10^10 cm2)')
parser.add_argument('-svld', '--solvsld', dest='solv_sld', default=4.922, type=float, help='Solvent SLD (10^-6 Å-2 == 10^10 cm2)')
parser.add_argument('-sc', '--solconc', dest='sol_conc', default=0.01, type=float, help='Solute concentration (g mL-1)')
parser.add_argument('-smw', '--solmw', dest='sol_mw', default=312.5, type=float, help='Solute concentration (g mL-1)')
parser.add_argument('-n', '--npts', dest='n', default=10, type=int, help='Number of low-Q datapoints for I(0) average')

args = parser.parse_args()

def aggfunc(filename=args.filename, solv_dens=args.solv_dens, solv_sld=args.solv_sld, sol_dens=args.sol_dens, sol_mw=args.sol_mw, sol_sld=args.sol_sld, sol_conc=args.sol_conc):
    na=6.023e+23
    df = pd.read_csv(f'{filename}', sep=',')
    i0 = np.average(df['I_mod'].head(args.n))
    vol_s = 1/args.sol_dens
    delta_sld = (args.solv_sld-args.sol_sld)*10000000000
    mw_agg = (i0*na) / ( (args.sol_conc * (vol_s**2) * (delta_sld**2)) )
    nagg = mw_agg / args.sol_mw
    print(f'{args.filename}')
    print(f'I(Q)0 = {round(i0,5)}')
    print(f'Aggregate Mw = {round(mw_agg,1)}')
    print(f'Nagg = {round(nagg,1)}')
    outdata = {'IQ0': [i0], 'MwAgg': [mw_agg], 'NAgg': [nagg]}
    dfo = pd.DataFrame(outdata)
    dfo.to_csv(f"{args.filename}.ag0", decimal='.', sep=' ')

aggfunc()
