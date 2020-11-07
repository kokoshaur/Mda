import numpy as np
import matplotlib.pyplot as plt

y = np.loadtxt("ansver.txt", delimiter='\t', dtype=np.float)
x = list(range(0, 8000))

plt.plot(x, y)
plt.show()
