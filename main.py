from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import threading

cut_off = 50000
background_cutoff = 4000


def flood_fill(x, y, val, data, closedset, step_size=1, threshold=0.01, always_up=False):
    if x >= data.shape[0] or y >= data.shape[1] or x < 0 or y < 0:
        return
    if (x, y) in closedset:
        return
    if np.abs(int(data[x, y]) - int(val)) / val <= threshold or (always_up and data[x, y] >= val):
        closedset.append((x, y))
    else:
        return
    flood_fill(int(x + step_size), int(y), val, data, closedset, step_size=step_size, threshold=threshold)
    flood_fill(int(x - step_size), int(y), val, data, closedset, step_size=step_size, threshold=threshold)
    flood_fill(int(x), int(y + step_size), val, data, closedset, step_size=step_size, threshold=threshold)
    flood_fill(int(x), int(y - step_size), val, data, closedset, step_size=step_size, threshold=threshold)


# masking
bleeding_edge = [
    {'tleft': (1428, 4608),
     'bright': (1445, 3509), 'name': 'Above Cen Star'},
    {'tleft': (1427, 2938),
     'bright': (1445, 0), 'name': 'Below Cen Star'},
    {'tleft': (1026, 370),
     'bright': (1703, 315), 'name': 'Xmas 1'},
    {'tleft': (1394, 279),
     'bright': (1475, 217), 'name': 'Xmas 2'},
    {'tleft': (1288, 164),
     'bright': (1521, 124), 'name': 'Xmas 3'},
    {'tleft': (1021, 359),
     'bright': (1702, 316), 'name': 'Xmas 4'},
    {'tleft': (1634, 61),
     'bright': (1717, 2), 'name': 'Xmas 5'},
    {'tleft': (1100, 442),
     'bright': (1642, 424), 'name': 'Xmas 6'},
    {'tleft': (1402, 469),
     'bright': (1482, 442), 'name': 'Xmas 7'},
    {'tleft': (1200, 3446),
     'bright': (1659, 2967), 'name': 'Central Star'},
    {'tleft': (725, 3426),
     'bright': (818, 3209), 'name': 'Other Star 1'},
    {'tleft': (865, 2358),
     'bright': (946, 2223), 'name': 'Other Star 2'},
    {'tleft': (935, 2837),
     'bright': (996, 2708), 'name': 'Other Star 3'},
    {'tleft': (19, 710),
     'bright': (112, 610), 'name': 'Other Star 4'},
    {'tleft': (2106, 3800),
     'bright': (2160, 3714), 'name': 'Other Star 5'},
]


# opening image with astropy


def main():
    def plotlogim(data):
        fig, ax = plt.subplots()
        sky = ax.imshow(np.log(data - 3420), origin="lower", cmap='viridis', aspect="equal")
        fig.colorbar(sky)
        plt.show()

    def plotlin(data):
        fig, ax = plt.subplots()
        sky = ax.imshow(data, origin="lower", cmap='jet', aspect="equal")
        fig.colorbar(sky)
        plt.show()
        # 3D plot of data points and counts

    def plotcount(data):
        fig = plt.figure()
        ax = plt.axes(projection="3d")
        X, Y = np.meshgrid(range(len(data[0])), range(len(data)))
        ax.plot_surface(X, Y, data)
        ax.set_zlim(0, 5000)
        ax.set_ylabel('Y')
        ax.set_xlabel('X')
        ax.set_zlabel('Pixel Count')
        plt.show()

        # function to convert counts to relative magnitudes and plot

    def mag(data):
        inst_mag_arr = np.zeros((len(data), len(data[0])))
        var1, var2 = len(data[0]), len(data)
        for y in range(var2):
            for x in range(var1):
                inst_mag_arr[y][x] = mag_known - 2.5 * np.log10(data[y][x])

        fig = plt.figure()
        ax = plt.axes(projection="3d")
        X, Y = np.meshgrid(range(len(data[0])), range(len(data)))
        ax.plot_surface(X, Y, inst_mag_arr)
        ax.set_ylabel('Y')
        ax.set_xlabel('X')
        ax.set_zlabel('Magnitude')
        plt.show()

    data_points = None

    with fits.open("A1_mosaic.fits") as hdulist:
        # for key, val in hdulist[0].header.items():
        #     print(f"{key},{val}")
        mag_known = hdulist[0].header['MAGZPT']
        mag_known_err = hdulist[0].header['MAGZRR']
        print(f"The known magnitude calibration is {mag_known}")
        print(f"The known magnitude calibration error is {mag_known_err}")
        data_points = hdulist[0].data

    for rect in bleeding_edge:
        tleft = np.array(rect["tleft"], dtype=int)
        bright = np.array(rect['bright'], dtype=int)
        # rows, columns
        data_points[bright[1]:tleft[1], tleft[0]:bright[0]] = 3419  # background value

    # remove the edges, first and last 150 columns
    # remove the edges, first and last 150 columns and first and last rows

    data_points = data_points[150:-150, 150:-150]
    width, height = data_points.shape
    # data_points = np.transpose(data_points)
    # remove bleeding edges
    plotlin(data_points)

    def trim_and_cut(data, cut):
        data = data[150:-150, 150:-150]

        def cutting(x):
            if x >= cut:
                return 3419
            return x

        data = np.vectorize(cutting)(data)
        return data

    data_points = data_points[150:-150, 150:-150]
    plotlin(data_points)
    plotlogim(data_points)

    # cut-off filter

    def cluster(fill_points):
        plt.figure()
        for centroid in fill_points:
            init_x, init_y = centroid[0], centroid[1]
            cluster_points = []
            flood_fill(init_x, init_y, data_points[init_x, init_y], data_points, cluster_points, step_size=10, threshold=0.2, always_up=False)
            cluster_points_trp = np.transpose(cluster_points)
            plt.scatter(cluster_points_trp[1], cluster_points_trp[0], s=1, label=f"Cluster")
        plt.legend()
        plt.show()

    #cluster(cluster_centroid)

    width, height = data_points.shape

    # # histogram of background radiation
    def histogram(data, max, min):
        n, bins, patches = plt.hist([x for x in data.flatten() if min < x < max], bins=max - min - 2)

        # fit guassian to find mean
        def gaus(x, a, x0, sigma):
            return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

        # finding the midpoints of the bins
        midpoints = [0] * (len(bins) - 1)
        for i in range(len(bins) - 1):
            midpoints[i] = (bins[i + 1] + bins[i]) / 2

        # initial guesses
        mean = 3420
        sigma = 50
        a = 2 * 10 ** 6

        x = midpoints
        y = n

        popt, pcov = curve_fit(gaus, x, y, p0=[a, mean, sigma])  # fitting

        perr = np.sqrt(np.diag(pcov))  # standard error of estimate

        print(f"p_error = {perr}")
        print(f"The amplitude is {popt[0]}, The mean is {popt[1]}, sigma = {popt[2]}")
        print(f"The error from sigma is estimated at: {popt[2] / np.sqrt(len(data))}")

        # plt.plot(x, gaus(x, *popt), 'ro:', label='Gaussian Fit')
        plt.figure(1)
        plt.plot(x, (y - gaus(x, *popt)) / y, label='Signal')
        # plt.xlim(3390, 3440)
        # plt.ylim(-0.1, 0.2)

        plt.xlabel('Pixel Value')
        plt.ylabel('Relative Frequency Offset From Gaussian')
        plt.title('Residual Nature - Highlights Local Background Regions')
        plt.legend()

        plt.show()

    #
    # showing raw image with masking elements

    def detect(data):
        loc = np.where(data == data.max())
        val = data[loc[0], loc[1]] #loc[0] is rows, loc[1] is columns
        if val.size >= 1:
            for i in range(len(loc[1])-1):  #checking row-wise, column wise to see if pixels next to each other
                if np.abs(loc[0][i]-loc[0][i+1])>=1 and np.abs(loc[1][i]-loc[1][i+1])>=1:
                    obj_arr = []
                    flood_fill(loc[1][i], loc[0][i], val, data, obj_arr, threshold=0.01)
                    tmp = 0
                    for i in obj_arr:
                        tmp += data_points[i]  # finding aperture flux of source

                    # now background around source
                    # find avg background count per pixel
                    # find total sum of pixels and subtract total amount of background contribution
                    # convert remaining flux into magnitude
                    # update boolean array to say this pixel has been dealt with

        pass

# mag(data_points)
if __name__ == "__main__":
    threading.stack_size(67108864)
    sys.setrecursionlimit(10 ** 5)
    thread = threading.Thread(target=main)
    thread.start()
