#!/usr/bin/env python3
import sys
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch

in_pickle = sys.argv[1]
data_dict = pickle.load(open(in_pickle, 'rb'))
print(data_dict)
fig = plt.figure()
ax = fig.add_axes([0,0,1,1], frameon=False, aspect=1.)
ax.text(0.5, 0.9, "PDBTM non-redundant", ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))
ax.text(0.5, 0.8, "No of proteins: {}\nNo of TM segments: {}".format(data_dict["total_prots"], data_dict["total_segments"]), ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))
ax.text(0.5, 0.7, "Only with TM segments of length >= 17", ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))
ax.text(0.5, 0.6, "No of proteins: {}\nNo of TM segments: {}".format(data_dict["total_prots_17"], data_dict["total_segments_17"]), ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))

ax.text(0.3, 0.5, "Full TM segments", ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))
ax.text(0.7, 0.5, "Trimmed TM segments", ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))

ax.text(0.3, 0.35, "With charged amino acids\nProteins: {}\nSegments: {}".format(data_dict["prots_with_charge"], data_dict["seg_with_charge"]), ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))
ax.text(0.7, 0.35, "With charged amino acids\nProteins: {}\nSegments: {}".format(data_dict["prots_with_charge_trimmed"], data_dict["seg_with_charge_trimmed"]), ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))

ax.text(0.3, 0.20, "".format(data_dict["prots_with_charge"], data_dict["seg_with_charge"]), ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))
ax.text(0.7, 0.20, "With charged amino acids\nProteins: {}\nSegments: {}".format(data_dict["prots_with_charge_trimmed"], data_dict["seg_with_charge_trimmed"]), ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))
# ax.text(0.5, 0.7, "Proteins with charged amino acids in TM segments:\n {}".format(data_dict["total_prots"], data_dict["total_segments"]), ha="center", bbox=dict(boxstyle="round", fc="w",ec="k"))

# p = mpatch.Circle((0.5, 0.6))
# ax.annotate("fancy", (0.5, 0.1), (0.5, 0.5), ha="center", va="center", arrowprops=dict(arrowstyle="fancy",patchB=p, shrinkA=5, shrinkB=5, fc="k", ec="k", connectionstyle="arc3,rad=-0.05",), bbox=dict(boxstyle="square", fc="w"))
plt.draw()
plt.show()
