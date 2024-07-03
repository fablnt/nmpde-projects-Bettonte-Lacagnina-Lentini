import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../PlotData.csv')

# Plot the data
plt.plot(data['Time'], data['avg(u)'], label='Data from Plot Over Data Filter')

# Add labels and title
plt.xlabel('Time step')
plt.ylabel('avg(u)')
plt.title('Replicated Plot from ParaView')

# Add legend
plt.legend()

# Show the plot
plt.show()

plt.savefig('PlotOverTime.png')
