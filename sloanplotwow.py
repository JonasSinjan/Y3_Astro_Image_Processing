import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

a = np.arange(10.75, 15.75, 0.5)
b = [ 1, 1, 5, 3, 14, 32, 75, 114, 279, 548, 1024]
cum = [sum(b[:i]) for i in range(1,len(b))]

plt.plot(a, np.log10(cum))
plt.show()

slope,intercept,r_val,p_val,err = linregress(a, np.log10(cum))
print(slope)