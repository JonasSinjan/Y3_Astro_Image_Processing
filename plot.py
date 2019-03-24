from Image import *

if __name__ == 'main':
    def main():
        # open csv and turn it into an array called catalogue
        # retrieve magnitudes and turn it into arr/list

        m = np.arange(9, 16, 0.5)
        N = [(len(list(filter(lambda x: x.mag < m_i, catalogue)))) for m_i in m]
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
    main_thread = threading.Thread(target=main)
    main_thread.start()
