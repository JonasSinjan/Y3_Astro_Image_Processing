from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
import sys, threading


def flood_fill(x, y, val, data, closedset, step_size=1, threshold=0.01, always_up=False, mask=None):
    if sys.getrecursionlimit() < 1000:
        raise RuntimeError("Insufficient recursion depth to allow the flood fill to commence")

    if x >= data.shape[0] or y >= data.shape[1] or x < 0 or y < 0:
        return
    if (x, y) in closedset:
        return
    if mask:
        if not mask[x, y]:
            return
    this_pt = int(data[x, y])
    val_i = int(val)
    if abs(this_pt - val_i) / val_i <= threshold or (always_up and data[x, y] >= val):
        closedset.append((x, y))
    else:
        return
    flood_fill(int(x + step_size), int(y), val, data, closedset, step_size=step_size, threshold=threshold, always_up=always_up, mask=mask)
    flood_fill(int(x - step_size), int(y), val, data, closedset, step_size=step_size, threshold=threshold, always_up=always_up, mask=mask)
    flood_fill(int(x), int(y + step_size), val, data, closedset, step_size=step_size, threshold=threshold, always_up=always_up, mask=mask)
    flood_fill(int(x), int(y - step_size), val, data, closedset, step_size=step_size, threshold=threshold, always_up=always_up, mask=mask)


class Image:

    def __init__(self, filename):
        self.filename = filename
        with fits.open(filename) as fits_file:
            self.known_magnitude = fits_file[0].header["MAGZPT"]
            self.known_magnitude_err = fits_file[0].header["MAGZR"]
            self.data = fits_file[0].data  # Rawdatafile in y,x
            self.mask = np.ones(self.data.shape, dtype=bool)
            self.height, self.width = self.data.shape
        self.boundary = 0

    def create_mask_map(self, cut_off, rect_masks=None, cluster_centroids=None):
        self.mask = (self.data <= cut_off)  # Automatically create mask map based off cut off
        if rect_masks:
            for rect in rect_masks:
                t_left = rect["tleft"] - self.boundary
                b_right = rect["bright"] - self.boundary
                self.mask[b_right[1]:t_left[1], t_left[0]:b_right[0]] = False
        if cluster_centroids:
            for centroid in cluster_centroids:
                init_x, init_y = centroid[1] - self.boundary, centroid[0] - self.boundary  # Need to flip x,y
                cluster_points = []
                flood_fill(init_y, init_y, self.data[init_x, init_y], self.data, cluster_points, 1, 0.05, True)
                for point in cluster_points:
                    self.mask[point[1], point[0]] = False

    def trim(self, boundary):
        self.boundary = boundary
        self.data = self.data[boundary:-boundary, boundary:-boundary]
        self.height, self.width = self.data.shape

    def create_catalogue(self):
        # Find brightest non masked object

        while True:
            sources = np.where(self.mask, self.data)
            if len(sources) <= 0:
                break
            peak_y, peak_x = np.unravel_index(sources.argmax(), sources.shape)
            peak_val = self.data[peak_y, peak_x]
            peak_points = []
            flood_fill(peak_y, peak_x, peak_val, self.data, peak_points, mask=self.mask)

            for point in peak_points:
                self.mask[point[1], point[0]] = False



    def cluster(self, fill_points):
        # Make sure that this is run on thread with additional stack memory available else this will likely fail!
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

        self.background = popt[1]
        self.background_sigma = popt[2]

    def filter_by_sigma(self, sigma_count=3):
        """
        Mask anything that is under sigma away from the background mean
        :param sigma:
        :return:
        """
        assert self.background
        assert self.background_sigma

        for y in range(len(self.data)):
            for x in range(len(self.data[y])):
                if self.data[y, x] - self.background < sigma_count * self.background_sigma:
                    self.mask[y, x] = False

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
        sky = ax.imshow(self.data * self.mask, origin="lower", cmap='jet', aspect="equal")
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


if __name__ == '__main__':
    def main():
        # Run all executable code here to ensure that
        # As Matplotlib is NOT thread safe running any plt commands outside of main may cause unexpected behaviour!
        pass


    sys.setrecursionlimit(10 ** 5)
    threading.stack_size(67108864)  # Largest possible stack size of 64MB on Windows
    main_thread = threading.Thread(target=main)
    main_thread.start()
