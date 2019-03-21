import matplotlib.patches as plt_patch
import matplotlib.pyplot as plt
import numpy as np


class StellarObject:
    class BoundingRect:
        def __init__(self, left, bottom, right, top):
            self.left = int(left)
            self.top = int(top)
            self.right = int(right)
            self.bottom = int(bottom)
            self.height = self.top - self.bottom
            self.width = self.right - self.left

        def get_patch(self, offset_x=0, offset_y=0):
            return plt_patch.Rectangle((self.left + offset_x, self.bottom + offset_y), self.width, self.height,
                                       linewidth=1, edgecolor='r',
                                       facecolor="none")

        def get_origin(self):
            return (self.left + self.right) / 2, (self.top + self.bottom) / 2

        def get_enc_points(self):
            return [(y, x) for y in range(self.bottom, self.top + 1) for x in range(self.left, self.right + 1)]

        @staticmethod
        def scale_rect_origin(obj, sf):
            width = obj.width * sf
            height = obj.height * sf
            origin_x, origin_y = obj.get_origin()
            return StellarObject.BoundingRect(origin_x - width / 2, origin_y - height / 2, origin_x + width / 2,
                                              origin_y + height / 2)

        def get_area(self):
            return self.height * self.width

    def __init__(self, points, peak_val):
        self.points = points  # Points in (y,x) format!!!
        self.peak_val = peak_val
        self.set_bouding_rect()

    def plot_me(self, data, mask):
        fig, ax = plt.subplots()
        y, x = self.points[0]
        ax.imshow(data[y - 30:y + 30, x - 30:x + 30] * mask[y - 30:y + 30, x - 30:x + 30], cmap="gray", origin="lower",
                  aspect="equal")
        ax.scatter(np.transpose(self.points)[1] - x + 30, np.transpose(self.points)[0] - y + 30, s=1)
        ax.add_patch(self.bounding_rect.get_patch(-x + 30, -y + 30))
        ax.add_patch(self.bg_bound.get_patch(-x + 30, -y + 30))
        plt.show()

    def set_bouding_rect(self):
        y_points, x_points = np.transpose(self.points)
        left = min(x_points) - 1
        bottom = min(y_points) - 1
        right = max(x_points) + 1
        top = max(y_points) + 1
        self.bounding_rect = self.BoundingRect(left, bottom, right, top)

    def get_background_rect(self, data, relative_mag, sf=1.5):
        self.bg_bound = self.BoundingRect.scale_rect_origin(self.bounding_rect, sf)
        bg_pts = filter(lambda x: x not in self.points, self.bg_bound.get_enc_points())
        total_count = [data[i] for i in self.bg_bound.get_enc_points()]
        bg_vals = [data[val] for val in bg_pts if data[val] <= 3500]
        self.local_background = np.mean(bg_vals)
        total_background_count = self.local_background * len(
            self.bg_bound.get_enc_points())  # total background counts for all pixels
        self.source_count = total_count - total_background_count  # counts just from object
        self.mag = relative_mag - 2.5*np.log10(self.source_count)
        #plt.hist(bg_vals)
        #plt.show()

    def set_magnitudes(self, instr_zero_pt):
        self.relative_magnitude = instr_zero_pt - 2.5 * np.log10(self.peak_val)

    def __str__(self):
        return str(self.points)
