from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib.colors as colors

from sunpy.net import attrs as a
from sunpy.net import Fido
from sunpy.map import Map
from sunpy.io.special import srs

import astropy.units as u


def load_map_304(date:datetime):
    query = (
    a.Time(date - timedelta(hours=1), date, date), 
    a.Instrument('AIA'),
    a.Wavelength(304 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1])[0]
    return Map(mappath)

def load_map_131(date:datetime):
    query = (
    a.Time(date - timedelta(hours=1), date, date), 
    a.Instrument('AIA'),
    a.Wavelength(131 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1])[0]
    return Map(mappath)

def load_srs_table(date: datetime):
    query = (
        a.Time(date - timedelta(days=1), date, date), 
        a.Instrument.soon,
    )
    result = Fido.search(*query)
    date_string = result[0, -1]['url'].split('/')[-1].split('.')[0][:-3]
    year, month, day  = (
        int(date_string[:4]),
        int(date_string[4:6]),
        int(date_string[6:])
    )
    retreived_date = datetime(year, month, day, 0,0,0) # day start
    srs_raw = Fido.fetch(result[0, -1])[0]
    srs_table = srs.read_srs(srs_raw)
    srs_table = srs_table[np.logical_or(srs_table['ID'] == 'I', srs_table['ID'] == 'IA')]
    # srs_table = srs_table[srs_table['ID'] == 'I']
    return srs_table, retreived_date

def plot_map_with_ARs(smap):
    fig = plt.figure()
    ax = fig.add_subplot(projection=smap)

    # Log transformation is required to view the image
    smap.plot_settings['norm'] = colors.LogNorm()

    smap.plot(axes=ax)
    smap.draw_limb(axes=ax)

    # Add a text box and arrow pointing to each active region
    lat_text = -40
    transparent_white = (1, 1, 1, 0.5)
    for num, lng, lat in zip(numbers, lngs.value, lats.value):
        ax.annotate(num, (lng, lat),
                    xytext=(320, lat_text),
                    xycoords=ax.get_transform('heliographic_stonyhurst'),
                    backgroundcolor=transparent_white,
                    color='red',
                    fontweight='bold',
                    arrowprops=dict(facecolor=transparent_white, width=1, headwidth=10),
                    horizontalalignment='right', verticalalignment='top')
        lat_text += 10
    plt.colorbar()

    plt.show()


now = datetime.now()
srs_table, retreived_date = load_srs_table(now)

lats = srs_table['Latitude']
lngs = srs_table['Longitude']
numbers = srs_table['Number']

smap = load_map_304(retreived_date)
plot_map_with_ARs(smap)
smap = load_map_131(retreived_date)
plot_map_with_ARs(smap)