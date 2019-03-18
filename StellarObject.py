import matplotlib.pyplot as plt
import numpy as np


class StellarObject:

    def __init__(self, points, peak_val):
        self.points = points  # Points in (y,x) format!!!
        self.peak_val = peak_val

    def plot_me(self, data, mask):
        fig, ax = plt.subplots()
        y, x = self.points[0]
        ax.imshow(data[y - 30:y + 30, x - 30:x + 30] * mask[y - 30:y + 30, x - 30:x + 30], cmap="gray", origin="lower", aspect="equal")
        ax.scatter(np.transpose(self.points)[1] - x + 30, np.transpose(self.points)[0] - y + 30, s=1)
        plt.show()

    def __str__(self):
        return str(self.points)
