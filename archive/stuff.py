import numpy as np
import matplotlib.pyplot as plt

# function f(x)
f = lambda x: 1.9 * x / (1.9 * x + 0.1 * (1 - x))

# # plot function
x = np.linspace(0, 1, 10)
# plt.plot(x, f(x))
# # 0.05 ticks and grid
# plt.xticks(np.arange(0, 1.1, 0.1))
# plt.yticks(np.arange(0, 1.1, 0.05))

# plt.grid()
# plt.show()

# plot the same function as histogram with height of bars representing the function value
# a green histogram on the bottom, a red on top of the green bars, pointing down, with (1-f(x)) height

# plot histogram
fig, ax = plt.subplots()
ax.bar(x, f(x), width=0.1, align='center', color='g')
ax.bar(x, 1-f(x), width=0.1, align='center', color='r', bottom=f(x))
# ticks and grid
plt.xticks(np.arange(0, 1.1, 0.1))
plt.yticks(np.arange(0, 1.1, 0.05))
plt.grid()
plt.show()

# make 6 subplots plots with Tmean, Tmin, Tmax


