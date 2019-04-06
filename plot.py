from Image import *
import matplotlib.patches as patches
import ast

if __name__ == '__main__':
    def main(filename):
        # open csv and turn it into an array called catalogue
        # retrieve magnitudes and turn it into arr/list
        img = Image("A1_mosaic.fits")
        img.create_mask_map(50000, rect_masks=bleeding_edge)
        img.trim(150)
        img.filter_by_sigma(5)
        fig, ax = plt.subplots(figsize=(6, 11), dpi=800)
        image_map = ax.imshow(np.log10(img.data * img.mask), origin="lower", cmap='gray', aspect="equal")

        df = pd.read_csv(filename)
        mag = df["Relative Magnitude"]

        star_points = df["Points"]

        for star in star_points:
            star = np.array(ast.literal_eval(star))
            y_points, x_points = np.transpose(star)
            #y, x = star[0]
            left = min(x_points) - 1
            bottom = min(y_points) - 1
            right = max(x_points) + 1
            top = max(y_points) + 1
            ax.add_patch(patches.Rectangle((left, bottom), right - left, top - bottom, linewidth=0.8, edgecolor='r',
                                           facecolor="none"))
            #ax.scatter(np.transpose(x_points) - x + 30, np.transpose(y_points) - y + 30, s=1)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        # ax.plot(1450, 3100, color='red', label='Detected Galaxies')
        # ax.legend()
        # fig.colorbar(image_map)
        ax.set_ylim(3920, 3975)
        ax.set_xlim(380, 440)
        plt.show()

        # m = np.arange(10, 25, 0.25)  # need to be wary of this range - could change for different sigma
        # N = [(len(list(filter(lambda x: x < m_i, mag)))) for m_i in m]
        # plt.plot(m, N)
        # plt.xlabel('Magnitude Limit')
        # plt.ylabel('Number of Galaxies')
        # plt.title(f'Number of galaxies against magnitude limit\n{filename}')
        # plt.show()
        #
        # m = np.arange(10, 17, 0.5)  # need to be wary of this range - could change for different sigma
        # N = [(len(list(filter(lambda x: x < m_i, mag)))) for m_i in m]
        # plt.figure()
        # slope, intercept, rvalue, pvalue, stderr = linregress(m, np.log10(N))
        # fit_str = f"Linear Regression Fit\nSlope:{round(slope, 3)}±{round(stderr, 3)}\nR^2:{round(rvalue ** 2, 3)}"
        # plt.plot(m, [i * slope + intercept for i in m], 'b-', label=fit_str)
        # plt.plot(m, np.log10(N), 'r.', linestyle='--', label='Raw Data')
        # plt.xlabel('Magnitude Limit')
        # plt.ylabel('Log_10 (Number of Galaxies)')
        #
        # plt.title(f'{filename}')
        # plt.legend()
        # plt.show()


    sys.setrecursionlimit(10 ** 5)
    threading.stack_size(67108864)  # Largest possible stack size of 64MB on Windows
    main_thread = threading.Thread(target=main, args=('survey_5sig_0.9_bugfix_circ_com_test.cat',))
    main_thread.start()
