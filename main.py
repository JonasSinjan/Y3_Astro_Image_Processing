from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

bleeding_edge = [
    {'tleft': (1428, 4608),
     'bright': (1445, 3509)}, # Top rectangle of main divider
    {'tleft': (1427, 2938),
     'bright': (1445, 0)}#Bottom rectangle
]

if __name__ == "__main__":

    data_points = None

    with fits.open("A1_mosaic.fits") as hdulist:
        for key, val in hdulist[0].header.items():
            print(f"{key},{val}")

        data_points = np.transpose(hdulist[0].data)

    for rect in bleeding_edge:
        tleft = np.array(rect["tleft"],dtype=int)
        bright = np.array(rect['bright'],dtype=int)
        print(data_points[tleft[0]:bright[0],:])
        print(data_points[:,bright[1]:tleft[1]])
        print(f"B4:{data_points[tleft[0]:bright[0], bright[1]:tleft[1]]}")
        data_points[tleft[0]:bright[0], bright[1]:tleft[1]] = 3400  # background value
        print(f"After:{data_points[tleft[0]:bright[0], bright[1]:tleft[1]]}")
    # remove the edges, first and last 150 columns

    data_points = data_points[150:-150, 150:-150]
    data_points = np.transpose(data_points)
    # remove bleeding edges

    fig,ax = plt.subplots()
    sky = ax.imshow(data_points,origin="lower")
    fig.colorbar(sky)
    plt.show()
