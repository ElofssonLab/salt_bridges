#!/usr/bin/env python3
import sys
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch

in_pickle = sys.argv[1]
title = sys.argv[2]
save_file = sys.argv[3]
data_dict = pickle.load(open(in_pickle, 'rb'))
# print(data_dict)
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0,0,1,1], frameon=False, aspect=1.)
# ax.arrow(0.48, 0.8, -0.18, -0.1)
# ax.text(0.5, 0.95, "PDBTM non-redundant", ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))





b15 = ax.annotate( "Same\ncharged pairs in step:\n1: $\\bf{{{}}}$\n2: $\\bf{{{}}}$\n3: $\\bf{{{}}}$\n4: $\\bf{{{}}}$\n5: $\\bf{{{}}}$".format(data_dict["num_same_pairs_in_1"], data_dict["num_same_pairs_in_2"], data_dict["num_same_pairs_in_3"], data_dict["num_same_pairs_in_4"], data_dict["num_same_pairs_in_5"]), (0.7, 0.05),ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
b14 = ax.annotate( "Oppositely\ncharged pairs in step:\n1: $\\bf{{{}}}$\n2: $\\bf{{{}}}$\n3: $\\bf{{{}}}$\n4: $\\bf{{{}}}$\n5: $\\bf{{{}}}$".format(data_dict["num_opp_pairs_in_1"], data_dict["num_opp_pairs_in_2"], data_dict["num_opp_pairs_in_3"], data_dict["num_opp_pairs_in_4"], data_dict["num_opp_pairs_in_5"]), (0.3, 0.05),ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
b13 = ax.annotate( "Charged pairs", (0.5, 0.20),ha="center", bbox=dict(boxstyle="round", fc="wheat",ec="k"))
b12 = ax.annotate( "Number of segments with #charge:\n1 charge:  $\\bf{{{}}}$\n2 charges: $\\bf{{{}}}$\n3 charges: $\\bf{{{}}}$\n4 charges: $\\bf{{{}}}$\n5 charges: $\\bf{{{}}}$".format(data_dict["seg_with_charge_in_1_trimmed"], data_dict["seg_with_charge_in_2_trimmed"], data_dict["seg_with_charge_in_3_trimmed"], data_dict["seg_with_charge_in_4_trimmed"], data_dict["seg_with_charge_in_5_trimmed"]), (0.7, 0.26),ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
b11 = ax.annotate( "Number of segments with #charge:\n1 charge:  $\\bf{{{}}}$\n2 charges: $\\bf{{{}}}$\n3 charges: $\\bf{{{}}}$\n4 charges: $\\bf{{{}}}$\n5 charges: $\\bf{{{}}}$".format(data_dict["seg_with_charge_in_1"], data_dict["seg_with_charge_in_2"], data_dict["seg_with_charge_in_3"], data_dict["seg_with_charge_in_4"], data_dict["seg_with_charge_in_5"]), (0.3, 0.26),ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
b10 = ax.annotate( "With multiple charged amino acids\nProteins: $\\bf{{{}}}$\nSegments: $\\bf{{{}}}$".format(data_dict["prots_with_multi_charge_trimmed"], data_dict["seg_with_multi_charge_trimmed"]), xytext=(0.7, 0.45),xy=(0.5,1), xycoords=b12, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"), arrowprops=dict(arrowstyle='fancy'))
b9 = ax.annotate( "With multiple charged amino acids\nProteins: $\\bf{{{}}}$\nSegments: $\\bf{{{}}}$".format(data_dict["prots_with_multi_charge"], data_dict["seg_with_multi_charge"]), xytext=(0.3, 0.45),xy=(0.5,1), xycoords=b11, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"), arrowprops=dict(arrowstyle='fancy'))
b8 = ax.annotate( "With charged amino acids\nProteins: $\\bf{{{}}}$\nSegments: $\\bf{{{}}}$".format(data_dict["prots_with_charge_trimmed"], data_dict["seg_with_charge_trimmed"]), xytext=(0.7, 0.55),xy=(0.5,1), xycoords=b10, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"), arrowprops=dict(arrowstyle='fancy'))
b7 = ax.annotate( "With charged amino acids\nProteins: $\\bf{{{}}}$\nSegments: $\\bf{{{}}}$".format(data_dict["prots_with_charge"], data_dict["seg_with_charge"]), xytext=(0.3, 0.55),xy=(0.5,1), xycoords=b9, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"), arrowprops=dict(arrowstyle='fancy'))
b6 = ax.annotate( "Trimmed TM segments", xytext=(0.7, 0.65),xy=(0.5,1), xycoords=b8, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="wheat",ec="k"), arrowprops=dict(arrowstyle='fancy'))
b5 = ax.annotate( "Full TM segments", xytext=(0.3, 0.65), xy=(0.5,1), xycoords=b7, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="wheat",ec="k"), arrowprops=dict(arrowstyle='fancy'))
b4b = ax.annotate( "No of proteins: $\\bf{{{}}}$\nNo of TM segments: $\\bf{{{}}}$".format(data_dict["total_prots_17"], data_dict["total_segments_17"]), xytext=(0.5, 0.71),xy=(0.5,1), xycoords=b6, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"), arrowprops=dict(arrowstyle='fancy',patchB=b6))
b4a = ax.annotate( "No of proteins: $\\bf{{{}}}$\nNo of TM segments: $\\bf{{{}}}$".format(data_dict["total_prots_17"], data_dict["total_segments_17"]), xytext=(0.5, 0.71),xy=(0.5,1), xycoords=b5, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"), arrowprops=dict(arrowstyle='fancy',patchB=b5))
b4 = ax.annotate( "No of proteins: $\\bf{{{}}}$\nNo of TM segments: $\\bf{{{}}}$".format(data_dict["total_prots_17"], data_dict["total_segments_17"]), (0.5, 0.71) ,ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
b3 = ax.annotate( "Only with TM segments of length >= 17", xytext=(0.5, 0.80),xy=(0.5,1), xycoords=b4, textcoords='data' ,ha="center", bbox=dict(boxstyle="round", fc="wheat",ec="k"), arrowprops=dict(arrowstyle='fancy'))
# ax.text(0.3, 0.20, "Segments with 1 charged residue: ()".format(data_dict["prots_with_charge"], data_dict["seg_with_charge"]), ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
# ax.text(0.7, 0.20, "With charged amino acids\nProteins: {}\nSegments: {}".format(data_dict["prots_with_charge_trimmed"], data_dict["seg_with_charge_trimmed"]), ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
# ax.text(0.5, 0.7, "Proteins with charged amino acids in TM segments:\n {}".format(data_dict["total_prots"], data_dict["total_segments"]), ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"))
b2 = ax.annotate( "No of proteins: $\\bf{{{}}}$\nNo of TM segments: $\\bf{{{}}}$".format(data_dict["total_prots"], data_dict["total_segments"]), xytext=(0.5, 0.86),xy=(0.5,1.1),xycoords=b3, textcoords='data', ha="center", bbox=dict(boxstyle="round", fc="lightyellow",ec="k"), arrowprops=dict(arrowstyle='fancy'))
b1 = ax.annotate(title, xytext=(0.5, 0.95), xy=(0.5,1), xycoords=b2, textcoords='data', ha="center", bbox=dict(boxstyle="round", fc="wheat",ec="k"), arrowprops=dict(arrowstyle='fancy'))

# p = mpatch.Circle((0.5, 0.6))
# ax.annotate("fancy", (0.5, 0.1), (0.5, 0.5), ha="center", va="center", arrowprops=dict(arrowstyle="fancy",patchB=p, shrinkA=5, shrinkB=5, fc="k", ec="k", connectionstyle="arc3,rad=-0.05",), bbox=dict(boxstyle="square", fc="lightyellow"))
# plt.draw()
# plt.show()
plt.savefig(save_file)
