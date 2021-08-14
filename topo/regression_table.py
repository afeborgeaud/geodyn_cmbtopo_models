import pandas as pd
import glob
import os

pd.options.display.float_format = '{:.2f}'.format

if __name__ == '__main__':
    figdir = '/work/anselme/topo_eth/synthetics/figures'
    dfs = [pd.read_csv(fname, sep=' ', header=None,
                       names=['model', 'phase', 'a3D', 'b3D', 'a1D', 'b1D'])
           for fname in glob.iglob(os.path.join(figdir, 'regression_*txt'))
           ]
    df = pd.concat(dfs, axis=0)\
        .pivot(index='phase', columns='model',
               values=['a3D', 'a1D'])\
        .reindex(['ScS', 'SKS', 'SKKS', 'PcP', 'PKP', 'PKKP', 'PcS'])
    # df.style.format('{:.2f}')

    print(df.head(10))
    outname = os.path.join(figdir, 'regression_table.txt')
    df.to_csv(outname, sep=' ', float_format='%.2f')
