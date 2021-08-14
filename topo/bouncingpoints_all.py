import pyshtools
import numpy as np
import os
os.environ['PROJ_LIB'] = r'/Users/navy/anaconda3/envs/ev02/share/proj'
import pygmt
from utils import read_dts
import glob
import logging


def parse_phase(fname):
    phases = ['ScS', 'SKS', 'SKKS', 'PcP', 'PKP', 'PKKP', 'PcS']
    for phase in phases:
        if phase in fname:
            return phase
    return None


def str_data(data):
    """
    Args:
        data (list[ndarray]):
    Returns:
        str: str for use in pygmt.plot(data=str)
    """
    return '\n>\n'.join(
        [
            '\n'.join([
                f'{str} {i % (max(1, len(datum) - 1))}' for i, str in
                enumerate(
                    np.apply_along_axis(lambda x: f'{x[0]:.2f} {x[1]:.2f}',
                                    axis=1,
                                    arr=datum)
                )
                ])
            for datum in data
        ]
    )



def display(fig, dts_dict, figsize=None, frame=["af", "WSne"],
            phases=['ScS', 'PcP', 'SKS', 'SKKS', 'PKP', 'PKKP', 'PcS']):
    proj = 'W15c'
    nrows, ncols = 4, 2
    with fig.subplot(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        frame=frame,
    ):
        for i, phase in enumerate(phases):
            irow = i // ncols
            icol = i % ncols
            print(irow, icol, phase)
            color = 'blue' if 'P' in phase else 'red'
            cmap = 'blue,green' if 'P' in phase else 'red,gold'

            fig.basemap(
                region='g', projection=proj,
                panel=[irow, icol], frame=frame
            )
            fig.coast(shorelines='1/0.5p,black', projection=proj,
                      panel=[irow, icol])
            # generate text data (much faster)
            data = []
            for index in dts_dict[phase].index:
                lonlats = np.vstack(
                    [list(x) for x in
                     dts_dict[phase].loc[index, 'points'].values()]
                )
                data.append(lonlats)
            data = str_data(data)
            datapath=f'tmp.txt'
            with open(datapath, 'w') as f:
                f.write(data)
            # plot
            fig.plot(
                data=datapath,
                cmap=cmap,
                # color=color,
                style="c0.2",
                pen="black",
                panel=[irow, icol],
                projection=proj
            )
            # fig.plot(
            #     data=datapath,
            #     pen=f'0.5p,{color}',
            #     panel=[irow, icol],
            #     projection=proj
            # )
            fig.text(
                x=0,
                y=90,
                offset='-7.5c/0c',
                justify='TR',
                text=phase,
                font='18p,Helvetica-Bold,black',
                no_clip=True,
                panel=[irow, icol],
                projection=proj
            )

    return fig


if __name__ == '__main__':
    root = '/work/anselme/topo_eth/synthetics/dts4/'
    figure_root = '/work/anselme/topo_eth/synthetics/figures'
    fnames = glob.iglob(
        os.path.join(root,
                     'correlationTimeshift*1DtopoTK10*.dat'))
    dts_dict = {parse_phase(fname): read_dts(fname)
                for fname in fnames}

    logging.debug(len(dts_dict))

    fig = pygmt.Figure()
    fig = display(fig, dts_dict, figsize=('32c', '30c'))

    output = os.path.join(os.path.expanduser(figure_root),
                          'bouncing_points_all.pdf')
    fig.savefig(output)
