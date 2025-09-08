import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('overlap.csv', delimiter=',', skip_header=0,
                     skip_footer=0, names=['n', 'h'])
data2 = np.genfromtxt('flight.csv', delimiter=',', skip_header=0,
                     skip_footer=0, names=['n', 'h'])

fig = plt.figure()

ax1 = fig.add_subplot(111)
ax1.semilogy(data['n'], data['h'], color='r', marker='.', linestyle='None', label='overlap')
ax1.semilogy(data2['n'], data2['h'], color='b', marker='.', linestyle = 'None', label='flight')
ax1.set_xlabel('N principal qn')
ax1.set_ylabel('N. or atoms [a.u.]')

plt.legend(bbox_to_anchor=(0.8, 0.5), loc=2, borderaxespad=0.)

plt.grid(True)

plt.show()
