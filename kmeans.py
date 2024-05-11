import sys

def submit_args():
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print("Invalid Input!")
        return 0
    try:
        k = int(sys.argv[1])
        if len(sys.argv) == 5:

            max_iter = int(sys.argv[2])
            input_file_name = sys.argv[3]
            output_file_name = sys.argv[4]
        else:
            max_iter = 200
            input_file_name = sys.argv[2]
            output_file_name = sys.argv[3]

        if type(k) != int or type(max_iter) != int or max_iter <= 0 or k <= 0:
            return 0

        f_input = open(input_file_name)
        f_output = open(output_file_name, 'w')
        f_input.close()
        f_output.close()

    except (ValueError,OSError):
        print("Invalid Input!")
        return 0

    return k, max_iter, input_file_name, output_file_name


def main():
    args = submit_args()
    if args == 0:
        return 0
    k, max_iter, input_file_name, output_file_name = args
    return k, input_file_name, output_file_name, max_iter

args = submit_args();
if args == 0:
    exit()
k, max_iter, input_filename,output_filename = args


# Determine the points from the file, the clusters and the initial means
points = []
clusters = [[] for i in range(k)]
try:
    f = open(input_filename)
except OSError:
    print ("Invalid input")
    exit()
with f:
    line = "temp"
    while line != "":
        line = f.readline()
        if line == "":
            continue
        line = line.split(",")
        points.append([float(i) for i in line])
f.close()
if len(points) == 0:
    print("Invalid input")
    exit()
if k >= len(points):
    print ("Invalid input")
    exit()
means = [points[i] for i in range(k)]


# Define helper functions
def euc_nor(vec):
    return (sum([i**2 for i in vec]))**0.5


def vec_sub_pow(vec1, vec2):
    return [(vec1[i] - vec2[i]) for i in range(max(len(vec1), len(vec2)))]


def cent_updt(vec_lst):
    if len(vec_lst) == 0:
        print ("An Error Has Occurred")
        exit()
    res = [0 for i in range(len(vec_lst[0]))]
    for dim in range(len(vec_lst[0])):
        res[dim] = sum([i[dim] for i in vec_lst])
        res[dim] = res[dim] / float(len(vec_lst))
    return res


# Define the k_means  function
def k_means(points, means, clusters):
    try:
        for i in range(len(points)):
            argmin = 10000
            bool_var = True
            idx = -1
            for j in range(len(clusters)):
                if (points[i] in clusters[j]):
                    clusters[j].remove(points[i])
            for j in range(len(means)):
                vec = vec_sub_pow(points[i], means[j])
                if euc_nor(vec) < argmin or bool_var:
                    bool_var = False
                    argmin = euc_nor(vec)
                    idx = j
            clusters[idx].append(points[i])
        for i in range(len(means)):
            means[i] = cent_updt(clusters[i])
        return [means, clusters]
    except IndexError:
        print ("An Error Has Occurred")
        exit()


# Run the function until we reach our goal
if __name__ == '_main_':
    main()

iter = 0
eps = 0.001
while iter < max_iter:
    cont = False
    iter += 1
    res = k_means(points, means, clusters)
    means = res[0]
    clusters = res[1]
    for i in means:
        if euc_nor(i) >= eps:
            cont = True
            break
    if cont:
        continue
    else:
        break

# Write to the file
with open(output_filename, 'w') as f:
    for i in means:
        output = ""
        for j in range(len(i)):
            temp = "%.4f" % i[j]
            output += (str(temp) + ",")
        output = output[:len(output) - 1]
        f.write(output)
        f.write('\n')
f.close()

