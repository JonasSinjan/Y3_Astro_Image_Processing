import matplotlib.pyplot as plt


class StellarObject:

    def __init__(self, points, peak_val):
        self.points = points
        self.peak_val = peak_val

    def plot_me(self, data):
        fig, ax = plt.subplots()
        y, x = self.points[0]
        ax.imshow(data[y - 30:y + 30, x - 30:x + 30], cmap="gray",origin="lower",aspect="equal")
        plt.show()

    def __str__(self):
        return str(self.points)
