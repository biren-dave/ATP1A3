import seaborn
import pandas as pd
import numpy as np

# names of columns
drug_conditions = ["marimastat 3 uM", "ibudilast 20 uM", "ibudilast 200 uM", 
                   "cabozantinib 0.1 uM", "cabozantinib 1 uM", "sorafenib 2 uM", "sorafenib 20 uM", 
                   "axitinib 2 uM", "axitinib 20 uM", "tofacitinib 2 uM", "tofacitinib 20 uM", 
                   "thalidomide 0.5 uM", "thalidomide 5 uM", "icatibant 0.1 uM", "icatibant 1 uM"]

# names of rows
cell_lines = ["G523", "G885", "G729", "G564", "G861"]              

# corresponding numerical value for each cell
viability_scores = np.array([[103.2445191,66.64440593,4.739848128,97.90725205,94.42274748,116.9796604,3.906406317,14.56179271,23.33680375,114.614175,117.3194148,92.78114457,104.8391006,99.06876532,99.90973465],
                             [99.54387203,95.89649106,4.679616554,99.55142491,99.892707,98.84118352,3.74480514,5.318624446,6.386772279,99.10143741,66.16415959,99.87811211,99.80021925,103.0671789,103.4136019],
                             [111.7364657,101.337013,13.41437035,115.9142502,109.4053095,116.1193251,9.227651647,84.10387886,55.32228376,115.9941438,115.0353046,93.96947089,112.1451494,94.90497327,109.030808],
                             [100.5542982,100.0080483,54.18581791,95.78376446,100.0588735,98.01859304,2.305208161,57.75735247,33.72282194,98.08674223,98.5647624,116.1383385,100.4812294,112.8439546,100.0594984],
                             [103.7184012,101.7273883,89.18527618,100.145608,61.24720474,103.5849188,2.019981229,46.41365913,4.292586744,101.9055264,89.33413613,102.7615483,100.2794762,111.1432579,105.007406]])
 
# compiles all of the above information into a pandas dataframe
data = pd.DataFrame(data=viability_scores, index=cell_lines, columns=drug_conditions)

print(data)

# make heirarchically clustered heatmap using seaborn
seaborn.set(color_codes=True)
clustered_map = seaborn.clustermap(data)

# save the generated heatmap
clustered_map.savefig("clustered.jpg")
