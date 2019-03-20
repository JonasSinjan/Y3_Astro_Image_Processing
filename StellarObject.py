import matplotlib.pyplot as plt
import matplotlib.patches as plt_patch
import numpy as np


class StellarObject:
    class BoundingRect:
        def __init__(self, left, bottom, right, top):
            self.left = left
            self.top = top
            self.right = right
            self.bottom = bottom
            self.height = top - bottom
            self.width = right - left

        def get_patch(self, offset_x=0, offset_y=0):
            return plt_patch.Rectangle((self.left + offset_x, self.bottom + offset_y), self.width, self.height, linewidth=1, edgecolor='r', facecolor="none")

        def get_area(self):
            return self.height*self.width

    def __init__(self, points, peak_val):
        self.points = points  # Points in (y,x) format!!!
        self.peak_val = peak_val
        self.set_bouding_rect()

    def plot_me(self, data, mask):
        fig, ax = plt.subplots()
        y, x = self.points[0]
        ax.imshow(data[y - 30:y + 30, x - 30:x + 30] * mask[y - 30:y + 30, x - 30:x + 30], cmap="gray", origin="lower", aspect="equal")
        ax.scatter(np.transpose(self.points)[1] - x + 30, np.transpose(self.points)[0] - y + 30, s=1)
        ax.add_patch(self.bounding_rect.get_patch(-x + 30, -y + 30))
        plt.show()

    def set_bouding_rect(self):
        y_points, x_points = np.transpose(self.points)
        left = min(x_points)
        bottom = min(y_points)
        right = max(x_points)
        top = max(y_points)
        self.bounding_rect = self.BoundingRect(left, bottom, right, top)



    def __str__(self):
        return str(self.points)
