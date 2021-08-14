import pyshtools
import numpy as np
import os
os.environ['PROJ_LIB'] = r'/Users/navy/anaconda3/envs/ev02/share/proj'
import spharm_io
import pygmt
import xarray as xr
from utils import read_station, read_event_format_2
import logging

def display_model(
        grid, fig, panel, proj='W15c', max_cpt=None,
        show_bar=True
):
    fig.basemap(
        region='g', projection=proj,
        panel=panel
    )
    if max_cpt is None:
        topo_max = float(np.abs(grid).max())
    else:
        topo_max = max_cpt
    with pygmt.config(
            COLOR_BACKGROUND='blue',
            COLOR_FOREGROUND='red'
    ):
        pygmt.makecpt(
            cmap='polar',
            series=[-topo_max, topo_max, 2 * topo_max / 10],
        )
    fig.grdimage(
        grid=grid,
        cmap=True,
        projection=proj,
        panel=panel
    )
    fig.coast(shorelines='1/0.5p,black', projection=proj,
              panel=panel)
    if show_bar:
        fig.colorbar(
            cmap=True,
            position="JBC+o0c/1.c+w7c/0.5ch",
            box=False,
            frame=["x+ltopo(km)"],
        )


def display_dataset(
        fig, stations, events, panel, proj='W15c',
        frame=["af", "WSne"]):
    fig.basemap(
        region='g', projection=proj,
        panel=panel, frame=frame
    )
    fig.coast(shorelines='1/0.5p,black', projection=proj,
              panel=panel)
    fig.plot(
        x=stations.lon,
        y=stations.lat,
        sizes=[0.3]*len(stations),
        color='cyan',
        style="i",
        pen="black",
        panel=panel,
        projection=proj
    )
    fig.plot(
        x=events.lon,
        y=events.lat,
        sizes=[0.6]*len(events),
        color='red',
        style="a",
        pen="black",
        panel=panel,
        projection=proj
    )


def parse_title(model_name):
    a = model_name.find('CMBtopo')
    b = model_name.find('_SHcoef')
    title = model_name[a+8:b]
    return title


def display_for_paper(stations, events, models, model_root, figure_root):
    grid_dict = dict()
    for model in models:
        model_path = os.path.join(os.path.expanduser(model_root), model)
        lmaxout, coeffs_r = spharm_io.read_sh_r_block_specfem(model_path)

        coeffs_r = np.pad(coeffs_r,
                          ((0, 0), (0, 80-lmaxout), (0, 80-lmaxout)),
                          'constant')

        title = None
        if 'CMBtopo-T' in model:
            title = parse_title(model)
        else:
            if 'lietal' in model:
                title = 'LGW91'
            if 'tanaka' in model:
                title = 'TK10'

        topo = pyshtools.expand.MakeGridDH(coeffs_r, sampling=2)
        print(f'{title} max abs topo: {np.abs(topo).max()}')
        grid = pyshtools.SHGrid.from_array(topo)
        topo = xr.DataArray(topo,
                            coords=[grid.lats(), grid.lons()],
                            dims=['lon', 'lat'])
        grid_dict[title] = topo

    models_order = [
        'T1-pPv', 'TC1-pPv', 'LGW91', 'TK10', 'TK10'
    ]
    fig = pygmt.Figure()
    proj = 'W15c'
    nrows = 3
    ncols = 2
    with fig.subplot(
        nrows=nrows,
        ncols=ncols,
        figsize=('32c', '29.5c'),
        frame=["af", "WSne"],
        autolabel='(a)'
    ):
        for i in range(nrows):
            for j in range(ncols):
                if i == j == 0:
                    display_dataset(fig, stations, events,
                                    [i, j], proj)
                else:
                    max_cpt = (None if i*ncols + j == nrows * ncols - 1
                               else 10)
                    grid = grid_dict[models_order[i*ncols + j - 1]]
                    display_model(grid, fig, [i, j], proj, max_cpt)

    output = os.path.join(os.path.expanduser(figure_root),
                          'models_and_dataset.pdf')
    fig.savefig(output)


if __name__ == '__main__':
    stations = read_station('resources/dataset/STATIONS_GLOBAL')
    events = read_event_format_2('resources/dataset/events.inf')

    figure_root = '../figures'
    model_root = 'resources/models/ylm'

    models = [
        'CMBtopo-T1-pPv_SHcoef_rot_l20.ylm',
        'CMBtopo-TC1-pPv_SHcoef_rot_l20.ylm',
        'cmb_lietal_1991.ylm',
        'tanaka10_lon.ylm',
    ]

    display_for_paper(stations, events, models, model_root, figure_root)
    logging.info(f"Saved figure to {figure_root}")
