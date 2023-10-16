from sunpy.net import attrs as a
from sunpy.net import Fido
import astropy.units as u
from sunpy.map import Map
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

dtime = datetime.now()

xflare_date = datetime(2011, 2, 15, 1, 55, 49)

tstart = '2017/07/02 00:00:00'
tend = '2017/07/12 00:00:00'
event_type = 'FL'
query = (
    a.Time(xflare_date - timedelta(minutes=45), xflare_date, xflare_date), 
    # a.hek.EventType(event_type),
    a.Instrument('AIA'),
    a.Wavelength(193 * u.Angstrom)
)
result = Fido.search(*query)
mappath = Fido.fetch(result[0, 0])[0]
m = Map(mappath)


figure = plt.figure(frameon=False)
ax = plt.axes([0, 0, 1, 1])
# Disable the axis
ax.set_axis_off()

# Plot the map.
# Since we are not interested in the exact map coordinates,
# we can simply use :meth:`~matplotlib.Axes.imshow`.
norm = m.plot_settings['norm']
norm.vmin, norm.vmax = np.percentile(m.data, [1, 99.9])
ax.imshow(m.data,
          norm=norm,
          cmap=m.plot_settings['cmap'],
          origin="lower")

plt.show()