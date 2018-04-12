# vectorize_diagrams.py

import numpy as np
import sys, json
from math import sqrt, atan, pi

def ComputePhi(X_grid, Y_grid, pd, p, sigma):
  # Number of points in the diagram
  num_points = len(pd);

  # Life span of generators
  x = np.array([ b          for (b,d) in pd])   # birth
  y = np.array([ abs(d - b) for (b,d) in pd])   # lifespan

  # Define the weight of y
  omega_y = np.power(y, p)

  Phi = np.zeros(X_grid.shape);

  for k in range(0, num_points):
    Phi = Phi + omega_y[k] * np.exp( -(np.pow(X_grid - x[k],2.0) + np.pow(Y_grid - y[k],2.0)) / (2.0 * sigma ** 2.0))

  # Weight function to set to zero at the boundary and
  # to smoothly increase to 1 far from the boundary
  W = (2.0 / pi) * atan(Y_grid);                    # Set to 0 if y is 0
  # W = (4 / pi^2) * atan(X_grid) .* atan(Y_grid);  # Set to 0 if x or y are 0

  # Multiply by W
  Phi = np.multiply(W,Phi)  # elementwise
  return Phi



helpstring = """
vectorize_diagrams input_diagram_filename output_vector_filename
"""
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(helpstring)
        exit(1)
    input_diagram_filename = sys.argv[1]
    output_vector_filename = sys.argv[2]
    # Number of bins
    #x_bins = 64; y_bins = 64;
    x_bins = 128; y_bins = 128;
    # x_bins = 256; y_bins = 256;

    p = 1.0;  # Exponent
    dim = 2;  # Dimension

    # File name for diagram
    x_min = 0
    y_min = 0
    x_max = 20
    y_max = 10

    # Compute x and y values
    x = np.linspace(x_min, x_max, x_bins + 1)
    y = np.linspace(y_min, y_max, y_bins + 1)

    # Compute the center points of the grid intervals
    x_grid = (x[0:-1] + x[1:]) / 2.0
    y_grid = (y[0:-1] + y[1:]) / 2.0

    # Compute 2D grids X and Y
    (X, Y) = np.meshgrid(x_grid, y_grid);

    # Grid sizes in x and y
    x_grid_h = x[1]-x[0];
    y_grid_h = y[1]-y[0];

    # Sigma for Gaussian
    sigma2_x = 1.0 * x_grid_h
    sigma2_y = 1.0 * y_grid_h
    sigma = sqrt(sigma2_x) * sqrt(sigma2_y);
    
    # Read diagram from file
    pd = json.load(input_diagram_filename);

    # Compute Gaussian
    Phi = ComputePhi(X, Y, pd[dim], p, sigma);

    # Transform Phi into a vector
    pers_vec = np.flatten(Phi)   # row-major

    with open(output_vector_filename, 'w') as outfile:
      json.dump(pers_vec, outfile)

