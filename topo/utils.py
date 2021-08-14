import pandas as pd
import json
import numpy as np


def read_dts(fname):
    columns = ['dt_rt', 'dt_cc', 'amp_ratio', 'cc', 'dist',
               'points', 'phase', 'sta', 'evt',
               'topo', 'dtopodx']
    dts = pd.read_csv(fname, sep=' ', header=None, names=columns,
                      converters={'points': lambda x: json.loads(x)})

    # dts = dts[(dts['dt_rt']-dts['dt_cc']).abs() < 3]
    return dts


def read_station(fname):
    stations = pd.read_csv(fname, sep=' ', header=None)
    stations.columns = ['sta', 'netwk', 'lat', 'lon', 'depth', 'elev']
    return stations


def read_event_format_1(fname='../dataset/eve-m5.8-pde'):
    events = np.genfromtxt(fname, delimiter=',',
                           dtype="U15,U25,f8,f8,f8,i4,f8",
                           names=['cat', 'date', 'lat', 'lon', 'depth', 'num',
                                  'mw'],
                           usecols=range(0, 7))
    return events


def read_event_format_2(fname):
    events = pd.read_csv(fname, sep=' ', header=None)
    events.columns = ['date', 'lat', 'lon', 'depth', 'Mw']
    return events


if __name__ == '__main__':
    df = read_dts(
        '/work/anselme/topo_eth/synthetics/dts4/'
        'correlationTimeshift_1DPREM_1Dtopolietall_SKS_cc0.95_proc.dat'
    )

    print(df.head())
