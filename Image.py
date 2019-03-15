from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

from main import flood_fill


class Image:

    def __init__(self, filename):
        self.filename = filename
        with fits.open(filename) as fits_file:
            self.known_magnitude = fits_file[0].header["MAGZPT"]
            self.known_magnitude_err = fits_file[0].header["MAGZR"]
            self.data = fits_file[0].data  # Rawdatafile in y,x
            self.height, self.width = self.data.shape
        self.boundary = 0

    def trim(self, boundary):
        self.boundary = boundary
        self.data = self.data[boundary:-boundary, boundary:-boundary]
        self.height, self.width = self.data.shape


    def cut(self, cutoff):
        self.data = np.vectorize(lambda x: self.background if x >= cutoff else x)(self.data)

    def cluster(self,fill_points):
        # Make sure that this is run on thread with additional stack memory availabile else this will likely fail!
        plt.xlim(0, self.width)
        plt.ylim(0, self.height)
        for centroid in fill_points:
            # Subtract 150 to normalise to new coordinate system
            init_x, init_y = centroid[1] - self.boundary, centroid[0] - self.boundary
            cluster_points = []
            flood_fill(init_x, init_y, self.data[init_x, init_y], self.data, cluster_points, step_size=1, threshold=0.25, always_up=True)
            cluster_points_trp = np.transpose(cluster_points)
            plt.scatter(cluster_points_trp[1], cluster_points_trp[0], s=1, label=f"Cluster")
        plt.show()

    def histogram(self, max, min):
        n, bins, patches = plt.hist([x for x in self.data.flatten() if min < x < max], bins=max - min - 2)

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
        print(f"The error from sigma is estimated at: {popt[2] / np.sqrt(len(self.data))}")
        self.background = popt[1]
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

    # Flush file to txt
    def flush_data_to_txt(self, txt_filename=None):
        filename = txt_filename if txt_filename else f"{self.filename}_dump.txt"
        with open(filename, mode="w") as txt_dump:
            txt_dump.write(f"({self.width},{self.height})\n")
            for row in self.data:
                for row_col in row:
                    txt_dump.write(str(row_col))
                txt_dump.write("\n")

    # Save changes in data table to the original file
    def update(self):
        with fits.open(self.filename, mode="update") as fits_file:
            fits_file[0].data = self.data
            fits_file.flush()

    def plotlogim(self):
        fig, ax = plt.subplots()
        sky = ax.imshow(np.log(self.data - self.background), origin="lower", cmap='viridis', aspect="equal")
        fig.colorbar(sky)
        plt.show()

    def plotlin(self):
        fig, ax = plt.subplots()
        sky = ax.imshow(self.data, origin="lower", cmap='jet', aspect="equal")
        fig.colorbar(sky)
        plt.show()

    def mag(self):
        inst_mag_arr = np.zeros((len(self.data), len(self.data[0])))
        var1, var2 = len(self.data[0]), len(self.data)
        for y in range(var2):
            for x in range(var1):
                inst_mag_arr[y][x] = self.known_magnitude - 2.5 * np.log10(self.data[y][x])

        fig = plt.figure()
        ax = plt.axes(projection="3d")
        X, Y = np.meshgrid(range(len(self.data[0])), range(len(self.data)))
        ax.plot_surface(X, Y, inst_mag_arr)
        ax.set_ylabel('Y')
        ax.set_xlabel('X')
        ax.set_zlabel('Magnitude')
        plt.show()
        # 3D plot of data points and counts
