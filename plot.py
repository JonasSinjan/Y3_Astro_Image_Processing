from Image import *
import matplotlib.patches as patches
import ast
if __name__ == '__main__':
    def main(filename):
        # open csv and turn it into an array called catalogue
        # retrieve magnitudes and turn it into arr/list
        img = Image("A1_mosaic.fits")
        img.trim(150)
        fig,ax = plt.subplots()
        image_map = ax.imshow(img.data * img.mask,origin="lower",aspect="equal")

        df = pd.read_csv(filename)
        mag = df["Relative Magnitude"]

        star_points = df["Points"]
        bound_rects = []

        for star in star_points:
            star = np.array(ast.literal_eval(star))
            y_points, x_points = np.transpose(star)
            left = min(x_points) - 1
            bottom = min(y_points) - 1
            right = max(x_points) + 1
            top = max(y_points) + 1
            ax.add_patch(patches.Rectangle((left,bottom),right-left,top-bottom))
        fig.colorbar(image_map)
        plt.show()

        m = np.arange(9, 16, 0.5)
        N = [(len(list(filter(lambda x: x < m_i, mag)))) for m_i in m]
        plt.plot(m, N)
        plt.show()

        plt.figure()
        slope, intercept, rvalue, pvalue, stderr = linregress(m, np.log10(N))
        fit_str = f"Linear Regression Fit\nSlope:{round(slope, 3)}±{round(stderr, 3)}\nR^2:{round(rvalue ** 2, 3)}"
        plt.plot(m, [i * slope + intercept for i in m], 'b-', label=fit_str)
        plt.plot(m, np.log10(N), 'r.', linestyle='--', label='Raw Data')
        plt.legend()
        plt.show()


    sys.setrecursionlimit(10 ** 5)
    threading.stack_size(67108864)  # Largest possible stack size of 64MB on Windows
    main_thread = threading.Thread(target=main, args=('survey_5sig.cat',))
    main_thread.start()
