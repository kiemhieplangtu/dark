import matplotlib.pyplot as plt
import numpy             as np

circle2=plt.Circle((5,5),.3,color='b',fill=False)

fig = plt.gcf()
ax = plt.gca()

# change default range so that new circles will work
ax.set_xlim((0,10))
ax.set_ylim((0,10))
# some data

fig.gca().add_artist(circle2)

plt.show()