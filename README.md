# KMeans++ Clustering

This project implements the k-means++ clustering algorithm using a combination of Python and C for improved performance. 
The code includes a Python-C API that allows the k-means algorithm to be called from Python, leveraging the speed of C.

## Project Structure

KMeans-2022/

├── kmeans.c # C code for the k-means algorithm

├── kmeans_plus_algo.py # Python script for preprocessing and calling the C function

├── setup.py # Script to build the Python-C extension

└── README.md # This file

## Prerequisites

- Python 3.x
- NumPy
- Pandas
- C compiler (e.g., gcc)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/MIKIHERSHCOVITZ/KMeans-2022.git
cd KMeans-2022
```

2. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate    # On Windows, use `venv\Scripts\activate`
```

3. Install the required Python packages:
```bash
pip install numpy pandas matplotlib
```

4. Build the Python-C extension:
```sh
python setup.py build_ext --inplace
```

## Usage

1. Prepare your input data files. Ensure you have two CSV files with your data. Each row should represent a data point, and each column should represent a feature. Example input files are provided in the example directory.

2. Run the kmeans_plus_algo.py script:

```sh
python kmeans_plus_algo.py k [max_iter] [eps] file_1 file_2
```

k: Number of clusters.
max_iter (optional): Maximum number of iterations (default is 300).
eps (optional): Epsilon value for convergence (default is 0.001).
file_1: Path to the first input file.
file_2: Path to the second input file.

## Example

To run the k-means++ algorithm with example data:
1. Create your input data files (see below for an example).
2. Execute the following command:
```sh
python kmeans_plus_algo.py 3 300 0.001 example/data1.csv example/data2.csv
```

## Results

# General Output
The script outputs the following information:

1. Initial Centroids Indices: The first line of the output shows the indices of the initial centroids selected by the k-means++ algorithm.
2. Final Centroids Coordinates: The subsequent lines display the final coordinates of the centroids after the k-means clustering algorithm has converged.

# Expected Result from Example Data
```sh
5,2,4
8.2000,9.3000,5.2000,5.9000
3.5000,3.1667,2.3667,1.5333
6.0000,7.1000,4.0000,4.8000
```

Initial Centroids Indices: 5,2,4 indicates that the initial centroids are chosen from the data points with indices 5, 2, and 4.

Final Centroids Coordinates: The final centroids are:
```sh
8.2000,9.3000,5.2000,5.9000
3.5000,3.1667,2.3667,1.5333
6.0000,7.1000,4.0000,4.8000
```