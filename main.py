from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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
    {'tleft': (1634, 61),
     'bright': (1717, 2), 'name': 'Xmas 5'},
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

if __name__ == "__main__":

    data_points = None

    with fits.open("A1_mosaic.fits") as hdulist:
        #for key, val in hdulist[0].header.items():
            #print(f"{key},{val}")

        data_points = hdulist[0].data

    for rect in bleeding_edge:
        tleft = np.array(rect["tleft"], dtype=int)
        bright = np.array(rect['bright'], dtype=int)
        # rows, columns
        data_points[bright[1]:tleft[1], tleft[0]:bright[0]] = 0  # background value

    # remove the edges, first and last 150 columns
    # remove the edges, first and last 150 columns and first and last rows

    data_points = data_points[150:-150, 150:-150]
    # data_points = np.transpose(data_points)
    # remove bleeding edges
max = 3500
min = 3350
n, bins, patches = plt.hist([x for x in data_points.flatten() if 3350 < x < 3500], bins=max-min-1)
plt.show()
#for x in data_points.flatten():
 #   if x > 50000:
  #      print(x)


# fit guassian to find mean
def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


midpoints = [0] * (len(bins) - 1)
for i in range(len(bins) - 1):
   midpoints[i] = (bins[i + 1] + bins[i]) / 2

mean = 3420
sigma = 50
a = 2 * 10 ** 6

x = midpoints
y = n

popt, pcov = curve_fit(gaus, x, y, p0=[a, mean, sigma])
perr = np.sqrt(np.diag(pcov))
print(f"p_error = {perr}")
print(f"The amplitude is {popt[0]}, The mean is {popt[1]}, sigma = {popt[2]}")
print(f"The error from sigma is estimated at: {popt[2]/np.sqrt(len(data_points))}")
#plt.plot(x, gaus(x, *popt), 'ro:', label='Gaussian Fit')
plt.figure(1)
plt.plot(x, (y-gaus(x, *popt))/y, label='Signal')
# plt.xlim(3390, 3440)
# plt.ylim(-0.1, 0.2)
plt.xlabel('Pixel Value')
plt.ylabel('Relative Frequency Offset From Gaussian')
plt.title('Residual Nature - Highlights Local Background Regions')
plt.legend()

plt.show()

#for i in y-gaus(x,*popt)

fig, ax = plt.subplots()
sky = ax.imshow(data_points, origin="lower", cmap='jet')
fig.colorbar(sky)
plt.show()
