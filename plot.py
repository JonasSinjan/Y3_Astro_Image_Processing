from Image import *
import matplotlib.patches as patches
import ast

if __name__ == '__main__':
    def main(filename):
        # open csv and turn it into an array called catalogue
        # eve magnitudes and turn it into arr/list
        img = Image("A1_mosaic.fits")
        img.header_dump()
        img.create_mask_map(50000, rect_masks=bleeding_edge)
        img.trim(150)
        img.filter_by_sigma(5)
        fig, ax = plt.subplots(figsize=(6, 11), dpi=800)
        image_map = ax.imshow(np.log10(img.data * img.mask), origin="lower", cmap='gray', aspect="equal")

        df = pd.read_csv(filename)
        mag = df["Relative Magnitude"]
        mag_err = df["Magnitude Error"]
        ratio = mag_err / mag * 100
        print(max(ratio), np.mean(ratio))
        star_points = df["Points"]

        # for star in star_points:
        #     star = np.array(ast.literal_eval(star))
        #     y_points, x_points = np.transpose(star)
        #     # y, x = star[0]
        #     left = min(x_points) - 1
        #     bottom = min(y_points) - 1
        #     right = max(x_points) + 1
        #     top = max(y_points) + 1
        #     ax.add_patch(patches.Rectangle((left, bottom), right - left, top - bottom, linewidth=0.4, edgecolor='g',
        #                                    facecolor="none"))
        #     # ax.scatter(np.transpose(x_points) - x + 30, np.transpose(y_points) - y + 30, s=1)
        # ax.set_xlabel("X")
        # ax.set_ylabel("Y")
        # # ax.plot(1450, 3100, color='red', label='Detected Galaxies')
        # # ax.legend()
        # # fig.colorbar(image_map)
        # ax.set_ylim(3915, 3980)
        # ax.set_xlim(380, 440)
        # # ax.set_ylim(3600, 4200)
        # # ax.set_xlim(0,600)
        # plt.show()
        # plt.figure()
        m = np.arange(9, 25, 0.25)  # need to be wary of this range - could change for different sigma

        mag_max = [0] * len(mag)
        mag_min = [0] * len(mag)
        for count, err in enumerate(mag_err):
            mag_max[count] = (mag[count] + err)
            mag_min[count] = (mag[count] - err)

        N_max = [(len(list(filter(lambda x: x < m_i, mag_min)))) for m_i in m]
        N_min = [(len(list(filter(lambda x: x < m_i, mag_max)))) for m_i in m]
        N = [(len(list(filter(lambda x: x < m_i, mag)))) for m_i in m]
        N_range = [(max_N - N_min[count]) / 2 for count, max_N in enumerate(N_max)]
        N_err = [np.sqrt(i) for i in N]
        tot_err_1 = [i + N_range[count] for count, i in enumerate(N_err)]
        plt.errorbar(m, N, yerr=tot_err_1, marker='.', markersize='2', color='blue', ecolor='red',
                     markerfacecolor='red')
        plt.xlabel('Magnitude Limit')
        plt.ylabel('Number of Galaxies')
        plt.title(f'Number of galaxies against magnitude limit\n{filename}')
        plt.show()
        #
        m_2 = np.arange(9.5, 17.0, 0.25)  # need to be wary of this range - could change for different sigma
        mag_max_2 = [0] * len(mag)
        mag_min_2 = [0] * len(mag)
        for count, err in enumerate(mag_err):
            mag_max_2[count] = (mag[count] + err)
            mag_min_2[count] = (mag[count] - err)

        N_max_2 = [(len(list(filter(lambda x: x < m_i, mag_min_2)))) for m_i in m_2]
        N_min_2 = [(len(list(filter(lambda x: x < m_i, mag_max_2)))) for m_i in m_2]
        N_2 = [(len(list(filter(lambda x: x < m_i, mag)))) for m_i in m_2]
        logN_2_err = [np.sqrt((1 / (np.log(10) * np.sqrt(i)))**2 + (1/np.log(10)*0.10)**2) for i in N_2]
        # print(len(N_2), len(logN_2_err))
        N_range_2 = [(np.log10(max_N) - np.log10(N_min[count])) / 2 for count, max_N in enumerate(N_max)]
        tot_err = [i + N_range_2[count] for count, i in enumerate(logN_2_err)]

        # here the problem is taking log10 of zero leads to infnite/nan for throws an error and code exits here
        # err_row_arr = np.stack(np.log10(N_max_2), np.log10(N_min_2))
        # log_err = np.log10(N_range)
        plt.figure()
        slope, intercept, rvalue, pvalue, stderr = linregress(m_2, np.log10(N_2))
        fit_str = f"Linear Regression Fit\nSlope:{round(slope, 3)}±{round(stderr, 3)}\nR^2:{round(rvalue ** 2, 3)}"
        plt.plot(m_2, [i * slope + intercept for i in m_2], 'b-', label=fit_str)
        plt.plot(m_2, [i * 0.33 + intercept - 0.3 for i in m_2], 'y-', label='max')
        plt.plot(m_2, [i * 0.285 + intercept + 0.2 for i in m_2], 'g-', label='min')
        # plt.plot(m_2, np.log10(N_2), 'r.', label='Data')
        plt.errorbar(m_2, np.log10(N_2), fmt='r.', yerr=tot_err, linestyle='--', label='Raw Data')
        plt.xlabel('Magnitude Limit')
        plt.ylabel('Log_10 (Number of Galaxies)')

        plt.title(f'{filename}')
        plt.legend()
        plt.show()


    def main_2(*filename):
        # just to plot gradient for different thresholds on same plot
        mag_list, mag_err_list = [], []
        for count, file in enumerate(filename):
            df = pd.read_csv(file)
            mag_list.append(df["Relative Magnitude"])
            # mag_err_list[count] = df["Magnitude Error"]

        name_list = [0.75, 0.80, 0.85, 0.90, 0.95]
        color_list = ['r', 'b', 'g', 'orange', 'black']
        m_2 = np.arange(9.5, 17.0, 0.25)  # need to be wary of this range - could change for different sigma

        for count, i in enumerate(mag_list):
            N = [(len(list(filter(lambda x: x < m_i, i)))) for m_i in m_2]
            slope, intercept, rvalue, pvalue, stderr = linregress(m_2, np.log10(N))
            fit_str = f"{name_list[count]}\nSlope:{round(slope, 3)}±{round(stderr, 3)}\nR^2:{round(rvalue ** 2, 3)}"
            plt.plot(m_2, [i * slope + intercept for i in m_2], color=color_list[count], label=fit_str)

        plt.xlabel('Magnitude Limit')
        plt.ylabel('Log_10 (Number of Galaxies)')
        plt.ylim(0.8, 3.0)
        plt.title('Variation of Threshold Parameter')
        plt.legend()
        plt.show()


    sys.setrecursionlimit(10 ** 5)
    threading.stack_size(67108864)  # Largest possible stack size of 64MB on Windows
    # main_thread = threading.Thread(target=main, args=(
    # 'survey_5sig_0.75.cat', 'survey_5sig_0.8.cat', 'survey_5sig_0.85.cat', 'survey_5sig_0.9_with_err.cat',
    # 'survey_5sig_0.95_with_err_newmask.cat',))
    main_thread = threading.Thread(target=main, args=('survey_5sig_0.95_with_err_newmask.cat',))
    main_thread.start()
