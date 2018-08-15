import matplotlib.pyplot as plt
import numpy as np

def comp():
    comp_scores = np.array([[5.054586795,5.542373418,5.648130398,5.308235868,5.808882526,5.493326454,5.089804852,5.435990281,4.642283189,4.528585479],
                            [5.001403227,5.386256501,5.377802597,5.269542766,5.269542766,5.111714747,4.66123133,5.328848907,4.575306807,4.457904625],
                            [3.777257486,5.04788977,4.514859489,4.865686213,4.865686213,4.728818194,4.63411965,5.071269272,4.66123133,4.492756753],
                            [4.898174467,5.56529646,4.971907235,5.163270196,4.49892713,4.119300131,4.277906557,4.565163446,3.740357066,3.779086377]])

    kernels = ["linear", "rbf", "polynomial", "sigmoid"]
    num_features = ["2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]

    fig, ax = plt.subplots(figsize=(20,10))
    im = ax.imshow(comp_scores)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(num_features)))
    ax.set_yticks(np.arange(len(kernels)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(num_features, fontsize=16)
    ax.set_yticklabels(kernels, fontsize=16)

    plt.xlabel("number of top features used", fontsize=20)
    plt.ylabel("type of kernel function", fontsize=20)

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("composite score", rotation=-90, va="bottom", fontsize=20)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), ha="center", rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(kernels)):
        for j in range(len(num_features)):
            text = ax.text(j, i, comp_scores[i, j], ha="center", va="center", color="w")

    #ax.set_title("Harvest of local farmers (in tons/year)")
    #fig.tight_layout()
    # plt.show()
    plt.savefig("clf_comp_performance.svg")

def r_comp():
    recur_scores = np.array([[5.505899922,6.111004657,6.111004657,5.925858277,6.253974618,6.10233506,5.948009507,5.219872871,5.263517528,5.263517528],
                             [5.59,6.10233506,6.10233506,5.919346808,5.934861258,6.468977739,6.114038654,5.934861258,5.550389216,5.48482536],
                             [3.319269696,5.355860251,5.355860251,5.546900732,5.546900732,4.958884739,5.546900732,5.919346808,5.162733477,5.35854358],
                             [4.926322044,6.111004657,5.433497327,5.59,4.853968206,4.462002844,4.896658064,4.410816229,3.79442724,4.137424703]])

    kernels = ["linear", "rbf", "polynomial", "sigmoid"]
    num_features = ["2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]

    fig, ax = plt.subplots(figsize=(20,10))
    im = ax.imshow(recur_scores)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(num_features)))
    ax.set_yticks(np.arange(len(kernels)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(num_features, fontsize=16)
    ax.set_yticklabels(kernels, fontsize=16)

    plt.xlabel("number of top features used", fontsize=20)
    plt.ylabel("type of kernel function", fontsize=20)

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("composite score", rotation=-90, va="bottom", fontsize=20)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), ha="center", rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(kernels)):
        for j in range(len(num_features)):
            text = ax.text(j, i, recur_scores[i, j], ha="center", va="center", color="w")

    #ax.set_title("Harvest of local farmers (in tons/year)")
    #fig.tight_layout()
    # plt.show()
    plt.savefig("clf_comp_performance_recurrent.svg")

comp()
r_comp()
