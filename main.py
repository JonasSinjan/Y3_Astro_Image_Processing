from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
data_points = np.array([])

with fits.open("A1_mosaic.fits") as hdulist:
    for key, val in hdulist[0].header.items():
        print(f"{key},{val}")

    data_points = hdulist[0].data

plt.hist(data_points)
plt.show()