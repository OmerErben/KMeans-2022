import numpy as np
import pandas as pd
import sys
import mykmeanssp
import matplotlib.pyplot as plt


def print_file(filename):
    with open(filename) as f:
        for line in f:
            print(line.strip())
        print()


def create_datafile_from_merged_inputs(filename, data):
    with open(filename, "w") as file:
        temp = data.reset_index()
        for index, row in temp.iterrows():
            lst_row = list(row[1:])
            cluster = [(str(dim)) for dim in lst_row]
            file.write(','.join(cluster) + "\n")


def create_clusters_file(filename, data):
    with open(filename, "w") as file:
        for vec in data:
            cluster = [(str(dim)) for dim in vec]
            file.write(",".join(cluster) + "\n")


def submit_args():
    if len(sys.argv) != 5 and len(sys.argv) != 6:
        print("Invalid Input! 1")
        return 1
    try:
        k = int(sys.argv[1])
        if len(sys.argv) == 6:
            max_iter = int(sys.argv[2])
            eps = float(sys.argv[3])
            file_1_name = sys.argv[4]
            file_2_name = sys.argv[5]
        else:
            max_iter = 300
            eps = float(sys.argv[2])
            file_1_name = sys.argv[3]
            file_2_name = sys.argv[4]

        if type(k) != int or type(max_iter) != int or max_iter <= 0 or k <= 0 or type(eps) not in [int, float] or eps < 0.0:
            print("Invalid Input!")
            return 1

        with open(file_1_name) as file_1, open(file_2_name) as file_2:
            pass  # Just to check if files can be opened
    except (ValueError, OSError):
        print("Invalid Input! 2")
        return 1
    return k, max_iter, eps, file_1_name, file_2_name


def extract_vec(df, num):
    choice = df[df.index == num]
    vec = choice.to_numpy().flatten()
    return vec


def euc_nor(vec):
    return np.sum(vec**2)


def sort_df(file_1, file_2):
    # Read the dataFrames
    try:
        df1 = pd.read_csv(file_1, header=None)
        df2 = pd.read_csv(file_2, header=None)
        # Set the first column to 'key'
        df1.rename(columns={0: 'key'}, inplace=True)
        df2.rename(columns={0: 'key'}, inplace=True)
        # Join the dfs on 'key'
        df = df1.merge(df2, on='key', how='inner')
        df.set_index('key', inplace=True)
        # sort the indicis
        df.sort_index(inplace=True)
    except Exception as e:
        print("Invalid Input! 3", e)
        return 1
    return df


def kmeans_pp(k, max_iter, eps, df):
    i = 1
    np.random.seed(0)
    create_datafile_from_merged_inputs("merged_input.txt", df)
    # Choose the index randomly and set it to a numpy array
    rows = df.index.to_numpy()

    if k > len(rows):
        print("invalid Input! 4")
        return 1

    rnd_num = np.random.choice(rows, 1)[0]
    rnd_vec = extract_vec(df, rnd_num)
    means = [rnd_vec]
    means_indices = [rnd_num]
    n = len(df.index)
    while i != k:
        D_lst = [0 for m in range(n)]
        P_lst = [0 for m in range(n)]
        # Create the d list containing all the Dl
        for l in range(n):
            min_dist = np.inf
            vec = extract_vec(df, rows[l])
            for j in range(i):
                cur_norm = euc_nor(vec - means[j])
                if cur_norm < min_dist:
                    min_dist = cur_norm
            D_lst[l] = min_dist
        # Calculate the probabilities
        D_sum = sum(D_lst)
        if D_sum == 0:
            print("Error: Sum of distances is zero. Check input data.")
            return 1
        for l in range(n):
            P_lst[l] = D_lst[l] / D_sum
        i += 1
        rnd_num = np.random.choice(rows, 1, p=P_lst)[0]
        rnd_vec = extract_vec(df, rnd_num)
        means.append(rnd_vec)
        means_indices.append(rnd_num)
    create_clusters_file("cluster_file.txt", means)
    return_lst = [str(i) for i in means_indices]
    return_str = ','.join(return_lst)
    print(return_str)
    return means  # Return the centroids


def plot_clusters(data, centroids, title='Clusters'):
    plt.scatter(data[:, 0], data[:, 1], c='blue', marker='o', label='Data Points')
    plt.scatter(centroids[:, 0], centroids[:, 1], c='red', marker='x', label='Centroids')
    plt.title(title)
    plt.xlabel('Feature 1')
    plt.ylabel('Feature 2')
    plt.legend()
    plt.show()


def main():
    args = submit_args()
    cluster_filename = "cluster_file.txt"
    data_filename = "merged_input.txt"
    if args == 1:
        return 1
    k, max_iter, eps, file_1, file_2 = args
    df = sort_df(file_1, file_2)
    if isinstance(df, int):
        return 1
    centroids = kmeans_pp(k, max_iter, eps, df)
    if centroids == 1:
        return 1
    try:
        kmeans_success = mykmeanssp.k_means(k, max_iter, eps, data_filename, cluster_filename)
    except Exception as e:
        print("An Error Has Occurred", e)
        return 1
    if kmeans_success == 0:
        print_file(cluster_filename)
        # Plot the data points and centroids
        data = df.to_numpy()[:, 1:]  # Assuming data points are in all columns except the first
        centroids = np.array(centroids)
        plot_clusters(data, centroids)
        return 0
    else:
        return 1


if __name__ == '__main__':
    main()
