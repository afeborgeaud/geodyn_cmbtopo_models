import pyshtools
from cartopy import crs as ccrs
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def parse_title(model_name):
    a = model_name.find('CMBtopo')
    b = model_name.find('_SHcoef')
    title = model_name[a + 8:b]
    return title


def display_cartopy(coeffs_r, ax, title=None, reverse_cmap=False,
                    cmap_amp=-1.):
    topo = pyshtools.expand.MakeGridDH(coeffs_r, sampling=2)
    grid_filt = pyshtools.SHGrid.from_array(topo)
    max_topo = np.abs(topo).max()
    # print('lie{} {:.2f}'.format(title, max_topo))

    # colormap
    if cmap_amp == -1:
        max_topo = np.round(max_topo, 1)
    else:
        max_topo = cmap_amp
    cmap = cm.get_cmap('bwr')
    if reverse_cmap:
        cmap = cmap.reversed()
    norm = mpl.colors.Normalize(-max_topo, max_topo)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)

    x = grid_filt.lons()
    y = grid_filt.lats()
    xx, yy = np.meshgrid(x, y)

    if title == 'TC4':  # fix for a weird bug when displaying model TC4
        topo -= 0.025
    elif title == 'TC1-pPv_eta1e-3':
        topo += 0.025
    elif title == 'T1-pPv':
        topo -= 0.08
    elif title == 'T1-pPv':
        topo -= 0.08

    ax.contourf(xx, yy, topo,
                cmap=cmap,
                norm=norm,
                transform=ccrs.PlateCarree(),
                extent=(-180, 180, -90, 90))

    cbar = plt.colorbar(sm, ax=ax, orientation='horizontal',
                        fraction=.03, pad=.1, aspect=40,
                        ticks=[-max_topo, 0, max_topo])
    cbar.ax.tick_params(labelsize=7)

    ax.coastlines()

    if title != None:
        ax.set_title(title, fontdict={'fontsize': 8})

    return cbar


def display_subplots_cartopy(
        coeffs_dict, nrows, ncols, central_longitude=-180,
        display_func=display_cartopy,
        keys_order=None, title=None, cbar_label=None, cmap_amp=-1.):
    proj = ccrs.Mollweide(central_longitude=central_longitude)
    fig, axes = plt.subplots(
        nrows, ncols, subplot_kw={'projection': proj})
    if keys_order == None:
        keys = sorted(coeffs_dict.keys())
    else:
        keys = keys_order
    for i, (model, ax) in enumerate(zip(keys, axes.ravel())):
        coeffs = coeffs_dict[model]
        sub_title = model
        if cmap_amp != -1.:
            cbar = display_func(
                coeffs, ax, title=sub_title, cmap_amp=cmap_amp)
            # if i < (nrows-1)*ncols:
            #    cbar.remove()
        else:
            cbar = display_func(coeffs, ax, title=sub_title)
        if (cbar_label != None) and (i >= (nrows - 1) * ncols):
            cbar.ax.set_xlabel(cbar_label, fontsize=8)
    if title != None:
        fig.suptitle(title, fontsize=10)
    return fig, axes
