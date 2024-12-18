import sys
import numpy as np
from kneed import KneeLocator
import pandas as pd

# Get the file path from command line arguments
data_file_path = sys.argv[1]

df = pd.read_csv(data_file_path, header=None)  # Read without a header
pca = df[0].tolist()  # Convert the first column to a list

# Using the KneeLocator to find the elbow point
kl = KneeLocator(range(1, len(pca) + 1), pca, curve='convex', direction='decreasing')
elbow_point = kl.elbow

print(elbow_point)