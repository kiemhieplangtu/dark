import matplotlib.pyplot as plt

# plt.axvline(x=0.22058956, linewidth=4, color='r')
# plt.axvline(x=0.33088437)
# plt.axvline(x=2.20589566)

xcoords = [0.22058956, 0.33088437, 2.20589566]
for xc in xcoords:
    plt.axvline(x=xc)

plt.show()