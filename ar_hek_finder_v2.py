import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import TimeDelta

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates import frames
from sunpy.net import attrs as a
from sunpy.net import hek
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.time import parse_time

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
hek_client = hek.HEKClient()
start_time = aia_map.date - TimeDelta(2*u.hour)
end_time = aia_map.date + TimeDelta(2*u.hour)
responses = hek_client.search(a.Time(start_time, end_time),
                              a.hek.AR, a.hek.FRM.Name == 'SPoCA')


area = 0.0

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map)
aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)

ar_set = set()
for i, response in enumerate(responses):
    #if response['area_atdiskcenter'] > area and np.abs(response['hgc_y']) < 80.0:
    #    area = response['area_atdiskcenter']

    response_index = i
    ar = responses[response_index]
    prev_size = len(ar_set)
    ar_id = int(ar['frm_specificid'].split('_')[-1])
    ar_set.add(ar_id)
    if prev_size == len(ar_set):
        continue

    p1 = ar["hpc_boundcc"][9:-2]
    p2 = p1.split(',')
    p3 = [v.split(" ") for v in p2]
    p4 = [(float(v[0]),float(v[1])) for v in p3]
    x_max = max([x for x,y in p4])
    x_min = min([x for x,y in p4])
    y_max = max([y for x,y in p4])
    y_min = min([y for x,y in p4])
    ar_center_x, ar_center_y = ar['event_coord1'], ar['event_coord2']
    x_range = x_max - x_min
    y_range = y_max - y_min

    #pp1 = list((ar_center_x + np.abs(0.5*x_range), ar_center_y + np.abs(0.5*y_range))) # upper right corner
    bottom_left = SkyCoord( (ar_center_x - np.abs(0.5*x_range)) * u.arcsec, (ar_center_y - np.abs(0.5*y_range)) * u.arcsec, frame=aia_map.coordinate_frame)
    top_right = SkyCoord( (ar_center_x + np.abs(0.5*x_range)) * u.arcsec, (ar_center_y + np.abs(0.5*y_range)) * u.arcsec, frame=aia_map.coordinate_frame)
    submap = aia_map.submap(bottom_left, top_right=top_right)
    e_start = str(ar['event_starttime'])
    e_start = e_start.replace(' ', '_').replace(':', '_').replace('.', '_')
    wavelength = submap.meta['wave_str']
    
    
    figg = plt.figure()
    axx = figg.add_subplot(projection=submap)
    image = submap.plot(axes=axx)
    submap.draw_limb(axes=axx)
    submap.draw_grid(axes=axx)

    # Make some room and put the title at the top of the figure
    axx.set_position([0.1, 0.1, 0.8, 0.7])
    axx.set_title(f'{axx.get_title()} AR {ar_id}', pad=45)
    # Please be aware that if you try to save this twice,
    # it will throw an exception rather than overwriting the file.
    submap.save(f'./{e_start}_AR{ar_id}_{wavelength}.fits', overwrite=True)
    #pp2 = list((ar_center_x - np.abs((0.5*x_range)), ar_center_y + np.abs((0.5*y_range)))) # upper left corner
    #pp3 = list((ar_center_x - np.abs((0.5*x_range)), ar_center_y - np.abs((0.5*y_range)))) # lower right corner
    #pp4 = list((ar_center_x + np.abs((0.5*x_range)), ar_center_y - np.abs((0.5*y_range)))) # lower left corner

    # ar_date = parse_time(ar['event_starttime'])

#ch_boundary = SkyCoord(
    # [(float(v[0]), float(v[1])) * u.arcsec for v in p3],
#    [v * u.arcsec for v in [pp1, pp2, pp3]],
#    obstime=ar_date, observer="earth",
#    frame=frames.Helioprojective)
# rotated_ch_boundary = solar_rotate_coordinate(ch_boundary, time=aia_map.date)
    aia_map.draw_quadrangle(bottom_left, axes=ax, width=x_range * u.arcsec, height=y_range * u.arcsec, linestyle="--", linewidth=1.5, edgecolor="red") 
# ax.plot_coord(rotated_ch_boundary, color='c')
# ax.set_title('{:s}\n{:s}'.format(aia_map.name, ar['frm_specificid']))
ax.set_title('{:s}'.format(aia_map.name))
# aia_map.draw_grid(axes=ax)

# plt.colorbar()

plt.show()