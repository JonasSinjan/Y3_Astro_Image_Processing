from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

bleeding_edge = [
    {'tleft': (1428, 4608),
     'bright': (1445, 3509)},  # Top rectangle of main divider
    {'tleft': (1427, 2938),    # Bottom rectangle
     'bright': (1445, 0)},
    {'tleft': (1026, 370),
     'bright': (1703, 315)},
    {'tleft': (1394, 279),
     'bright': (1475, 217)},
    {'tleft': (1288, 164),
     'bright': (1521, 124)},
    {'tleft': (1021, 359),
     'bright': (1702, 316)},
    {'tleft': (1634, 61),
     'bright': (1717, 2)},
    {'tleft': (1200, 3446),
     'bright': (1659, 2967)},

]

if __name__ == "__main__":

    data_points = None

    with fits.open("A1_mosaic.fits") as hdulist:
        for key, val in hdulist[0].header.items():
            print(f"{key},{val}")

        data_points = hdulist[0].data

    for rect in bleeding_edge:
        tleft = np.array(rect["tleft"], dtype=int)
        bright = np.array(rect['bright'], dtype=int)
        # rows, columns
        data_points[bright[1]:tleft[1], tleft[0]:bright[0]] = 0 # background value

    # remove the edges, first and last 150 columns
    # remove the edges, first and last 150 columns and first and last rows

    data_points = data_points[150:-150, 150:-150]
    # data_points = np.transpose(data_points)
    # remove bleeding edges

plt.hist([x for x in data_points.flatten() if 3000 < x < 7000], bins=100)
plt.show()

fig, ax = plt.subplots()
sky = ax.imshow(data_points, origin="lower", cmap='jet')
fig.colorbar(sky)
plt.show()
