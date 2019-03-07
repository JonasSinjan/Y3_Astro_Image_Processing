from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

data_points = np.array([])

with fits.open("A1_mosaic.fits") as hdulist:
    for key, val in hdulist[0].header.items():
        print(f"{key},{val}")

    data_points = hdulist[0].data

# remove the edges, first and last 150 columns

data_points = data_points[150:-150, 150:-150]

plt.hist(data_points)
plt.show()
