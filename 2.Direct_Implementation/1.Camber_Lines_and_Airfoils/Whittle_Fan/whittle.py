# Whittle fan file for random computations

import pandas as pd

# Add half the local blade pitch to the coordinates of the blade
filename_lower = "Whittle_lower.txt"
filename_upper = "Whittle_upper_alignedendpoints.txt"

lower = pd.read_csv(filename_lower, sep="\s+", names=["x", "y"])
upper = pd.read_csv(filename_upper, sep="\s+", names=["x", "y"])

pitch = 0.05105

lower["y"] = lower["y"] + pitch/2.0
upper["y"] = upper["y"] + pitch/2.0

lower_file = open("whittle_lower_shifted.txt", 'w')
upper_file = open("whittle_upper_shifted.txt", 'w')

lower_file.write(lower)
upper_file.write(upper)

lower_file.close()
upper_file.close()
