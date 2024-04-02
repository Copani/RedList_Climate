import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import numpy as np

np.random.seed(19680801)
fig = plt.figure()

axs = fig.subfigures(2, 2)

pts = np.random.rand(30)*.2
# Now let's make two outlier points which are far away from everything.
pts[[3, 14]] += .8

# If we were to simply plot pts, we'd lose most of the interesting
# details due to the outliers. So let's 'break' or 'cut-out' the y-axis
# into two portions - use the top (ax1) for the outliers, and the bottom
# (ax2) for the details of the majority of our data
ax = axs[0,0].subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)  # adjust space between axes

# plot the same data on both axes
ax[0].plot(pts)
ax[1].plot(pts)

# zoom-in / limit the view to different portions of the data
ax[0].set_ylim(.78, 1.)  # outliers only
ax[1].set_ylim(0, .22)  # most of the data

# hide the spines between ax and ax2
ax[0].spines.bottom.set_visible(False)
ax[1].spines.top.set_visible(False)
ax[0].xaxis.tick_top()
ax[0].tick_params(labeltop=False)  # don't put tick labels at the top
ax[1].xaxis.tick_bottom()

# Now, let's turn towards the cut-out slanted lines.
# We create line objects in axes coordinates, in which (0,0), (0,1),
# (1,0), and (1,1) are the four corners of the axes.
# The slanted lines themselves are markers at those locations, such that the
# lines keep their angle and position, independent of the axes size or scale
# Finally, we need to disable clipping.

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax[0].plot([0, 1], [0, 0], transform=ax[0].transAxes, **kwargs)
ax[1].plot([0, 1], [1, 1], transform=ax[1].transAxes, **kwargs)

ax01 = axs[0,1].subplots(1, 1)
ax01.plot(pts)

plt.show()