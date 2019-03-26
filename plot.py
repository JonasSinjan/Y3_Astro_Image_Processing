from Image import *
import matplotlib.patches as patches
import ast

if __name__ == '__main__':
    def main(filename, filename2, filename3, filename4):
        # open csv and turn it into an array called catalogue
        # retrieve magnitudes and turn it into arr/list
        img = Image("A1_mosaic.fits")
        img.trim(150)
        fig, ax = plt.subplots(figsize=(6, 11), dpi=800)
        image_map = ax.imshow(np.arcsinh(img.data*img.mask), origin="lower", cmap="gray", aspect="equal")

        df = pd.read_csv(filename)
        mag = df["Relative Magnitude"]
        df2 = pd.read_csv(filename2)
        mag2 = df2["Relative Magnitude"]
        df3 = pd.read_csv(filename3)
        mag3 = df3["Relative Magnitude"]
        df4 = pd.read_csv(filename4)
        mag4 = df4["Relative Magnitude"]

        df_list = [df, df2, df3, df4]
        color_list = ['r', 'b', 'y', 'g']
        sig_list = [2.5, 3, 4, 5]
        for count, df in enumerate(df_list):
            star_points = df["Points"]
            color = color_list[count]
            sig = sig_list[count]
            ax.plot(1250, 3050, color=color, label=f'{sig}')
            for star in star_points:
                star = np.array(ast.literal_eval(star))
                y_points, x_points = np.transpose(star)
                left = min(x_points) - 1
                bottom = min(y_points) - 1
                right = max(x_points) + 1
                top = max(y_points) + 1
                ax.add_patch(patches.Rectangle((left, bottom), right - left, top - bottom, linewidth=1, edgecolor=color,
                                               facecolor="none"))
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_xlim(1500, 1700)
        ax.set_ylim(2800, 3300)
        plt.legend()
        #fig.colorbar(image_map)
        plt.show()

        #
        # plt.show()

        # need to be wary of this range - could change for different sigma
        # N = [(len(list(filter(lambda x: x < m_i, mag)))) for m_i in m]
        # N2 = [(len(list(filter(lambda x: x < m_i, mag2)))) for m_i in m]
        # N3 = [(len(list(filter(lambda x: x < m_i, mag3)))) for m_i in m]
        # N4 = [(len(list(filter(lambda x: x < m_i, mag4)))) for m_i in m]
        # plt.plot(m, N)
        # plt.plot(m, N2)
        # plt.plot(m, N3)
        # plt.plot(m, N4)
        # plt.show()
        # mag_list = [mag, mag2, mag3, mag4]
        # file_list = [filename, filename2, filename3, filename4]
        # m = np.arange(10, 16.5, 0.5)
        #
        # for count, mag in enumerate(mag_list):
        #     plt.figure()
        #     N = [(len(list(filter(lambda x: x < m_i, mag)))) for m_i in m]
        #     slope, intercept, rvalue, pvalue, stderr = linregress(m, np.log10(N))
        #     fit_str = f"Linear Regression Fit\nSlope:{round(slope, 3)}±{round(stderr, 3)}\nR^2:{round(rvalue ** 2, 3)}"
        #     plt.plot(m, [i * slope + intercept for i in m], 'b-', label=fit_str)
        #     filename = file_list[count]
        #     plt.plot(m, np.log10(N), 'r.', linestyle='--', label=f'{filename}')
        #     plt.xlabel('m')
        #     plt.ylabel('Log10(N)')
        #     plt.legend()
        # plt.show()


    sys.setrecursionlimit(10 ** 5)
    threading.stack_size(67108864)  # Largest possible stack size of 64MB on Windows
    main_thread = threading.Thread(target=main, args=(
        'survey_2.5sig_0.85.cat', 'survey_3sig_0.85.cat', 'survey_4sig_0.85.cat', 'survey_5sig_0.85.cat'))
    main_thread.start()
