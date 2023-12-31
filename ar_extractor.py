import os
from datetime import timedelta, datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import TimeDelta, Time

# import sunpy.data.sample
# import sunpy.map
from sunpy.net import Fido
from sunpy.map import Map
from sunpy.coordinates import frames
from sunpy.net import attrs as a
from sunpy.io.special import srs
from sunpy.net import hek
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.time import parse_time



# implementacion singleton para el HEK client
hek_client = None

def get_hek_client_or_create():
    global hek_client
    if hek_client:
        return hek_client
    else:
        return hek.HEKClient()

def hek_search(start_time, end_time):
    global hek_client
    hek_client = get_hek_client_or_create()
    return hek_client.search(a.Time(start_time, end_time),
                             a.hek.AR, a.hek.FRM.Name == 'SPoCA') 

def load_map_304(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(304 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/304/')[0]
    return Map(mappath)

def load_map_94(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(94 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/94/')[0]
    return Map(mappath)

def load_map_335(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(335 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/335/')[0]
    return Map(mappath)

def load_map_211(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(211 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/211/')[0]
    return Map(mappath)

def load_map_193(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(193 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/193/')[0]
    return Map(mappath)

def load_map_171(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(171 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/171/')[0]
    return Map(mappath)

def load_map_131(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(131 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/131/')[0]
    return Map(mappath)

def load_map_1700(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(1700 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/1700/')[0]
    return Map(mappath)

def load_map_1600(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(1600 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/1600/')[0]
    return Map(mappath)

def load_map_4500(date):
    query = (
    a.Time( date, date + timedelta(minutes=45), date), 
    a.Instrument('AIA'),
    a.Wavelength(4500 * u.Angstrom)
    )
    result = Fido.search(*query)
    mappath = Fido.fetch(result[0, -1], path=f'./data/4500/')[0]
    return Map(mappath)

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


# cargando los datos de fulguraciones detectadas por el SWPC-GOES
df = pd.read_csv('./data/flares-goes-x-ray-unified.dat', sep='\t')
df['t-inicio'] = pd.to_datetime(df['t-inicio'])
df['t-max'] = pd.to_datetime(df['t-max'])
df['t-fin'] = pd.to_datetime(df['t-fin'])
df = df.drop(['unknown'], axis=1)

# fecha y hora de los máximos de las fulguraciones
max_dates_swpc = df['t-max']
# solamente queremos fechas para las que hay datos del SDO (2010)
max_dates_swpc = max_dates_swpc[max_dates_swpc.dt.tz_localize(None) > np.datetime64('2013-06-20 00:00:00')]

for max_date in max_dates_swpc:
    start_time = Time(max_date) 
    end_time = Time(max_date) + TimeDelta(45*u.min)
    responses  = hek_search(start_time, end_time)
    ar_set = set()

    # descarga los FITS para esta fecha
    try:
        smap_131 = load_map_131(max_date)
        smap_304 = load_map_304(max_date)
        smap_94 = load_map_94(max_date)
        smap_335 = load_map_335(max_date)
        smap_211 = load_map_211(max_date)
        smap_193 = load_map_193(max_date)
        smap_171 = load_map_171(max_date)
        smap_1600 = load_map_1600(max_date)
        smap_1700 = load_map_1700(max_date)
        smap_4500 = load_map_4500(max_date)

        srs_table, retreived_date = load_srs_table(start_time)

        lats = srs_table['Latitude']
        lngs = srs_table['Longitude']
        numbers = srs_table['Number']

        # plot_map_with_ARs(smap_131)
    except IndexError as e:
        print(f'No hay datos para esta fecha: {max_date}')
        continue
    
    # for i, response in enumerate(responses):
    #     response_index = i
    #     ar = responses[response_index]
    #     prev_size = len(ar_set)
    #     ar_id = int(ar['frm_specificid'].split('_')[-1])
    #     ar_set.add(ar_id)
    #     if prev_size == len(ar_set):
    #         continue
        
    #     p1 = ar["hpc_boundcc"][9:-2]
    #     p2 = p1.split(',')
    #     p3 = [v.split(" ") for v in p2]
    #     p4 = [(float(v[0]),float(v[1])) for v in p3]
        
    #     x_max = max([x for x,y in p4])
    #     x_min = min([x for x,y in p4])
    #     y_max = max([y for x,y in p4])
    #     y_min = min([y for x,y in p4])
    #     ar_center_x, ar_center_y = ar['event_coord1'], ar['event_coord2']
    #     x_range = x_max - x_min
    #     y_range = y_max - y_min

    #     #pp1 = list((ar_center_x + np.abs(0.5*x_range), ar_center_y + np.abs(0.5*y_range))) # upper right corner
    #     bottom_left = SkyCoord( (ar_center_x - np.abs(0.5*x_range)) * u.arcsec, (ar_center_y - np.abs(0.5*y_range)) * u.arcsec, frame=smap_131.coordinate_frame)
    #     top_right = SkyCoord( (ar_center_x + np.abs(0.5*x_range)) * u.arcsec, (ar_center_y + np.abs(0.5*y_range)) * u.arcsec, frame=smap_131.coordinate_frame)
        
        
    #     submap_131 = smap_131.submap(bottom_left, top_right=top_right)
    #     submap_304 = smap_131.submap(bottom_left, top_right=top_right)
    #     submap_94 = smap_94.submap(bottom_left, top_right=top_right)
    #     submap_171 = smap_171.submap(bottom_left, top_right=top_right)
    #     # submap_211 = smap_211.submap(bottom_left, top_right=top_right)
    #     # submap_193 = smap_193.submap(bottom_left, top_right=top_right)
    #     # submap_335 = smap_335.submap(bottom_left, top_right=top_right)
    #     # submap_1600 = smap_1600.submap(bottom_left, top_right=top_right)
    #     # submap_1700 = smap_1700.submap(bottom_left, top_right=top_right)
    #     # submap_4500 = smap_4500.submap(bottom_left, top_right=top_right)

    #     e_start = str(ar['event_starttime'])
    #     e_start = e_start.replace(' ', '_').replace(':', '_').replace('.', '_')
    #     wavelength_131 = submap_131.meta['wave_str']
    #     wavelength_304 = submap_304.meta['wave_str']
    #     wavelength_171 = submap_171.meta['wave_str']
    #     wavelength_94 = submap_94.meta['wave_str']
    #     # wavelength_193 = submap_193.meta['wave_str']
    #     # wavelength_211 = submap_211.meta['wave_str']
    #     # wavelength_335 = submap_335.meta['wave_str']
    #     # wavelength_1600 = submap_1600.meta['wave_str']
    #     # wavelength_1700 = submap_1700.meta['wave_str']
    #     # wavelength_4500 = submap_4500.meta['wave_str']

    #     submap_131.save(f'./{e_start}_AR{ar_id}_{wavelength_131}.fits', overwrite=True)
    #     submap_304.save(f'./{e_start}_AR{ar_id}_{wavelength_304}.fits', overwrite=True)
    #     submap_171.save(f'./{e_start}_AR{ar_id}_{wavelength_171}.fits', overwrite=True)
    #     submap_94.save(f'./{e_start}_AR{ar_id}_{wavelength_94}.fits', overwrite=True)



        # break
    # break
