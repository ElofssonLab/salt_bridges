#!/usr/bin/env python3
import math
import sys

c0 = 0.27045
c1 = 9.29274167549645
c2 = -0.64513139783394
c3 = 0.00822196628688
piover180 = math.atan2(1, 1) / 45

aas = 'ABCDEFGHIKLMNPQRSTUVWXYZ'
aaProfiles = {'A': [0.1267255, 0.0215152],
              'B': [1.5556351, 0.0133523],
              'C': [-0.0765051, 0.0994228],
              'D': [1.7939795, 0.0172643],
              'E': [1.4193720, 0.0089351],
              'F': [-0.2766953, 0.0010297],
              'G': [0.4813492, 0.0047210],
              'H': [1.1998590, 0.0080127],
              'I': [-0.4597384, 0.0181495],
              'K': [1.8485768, 0.0218446],
              'L': [-0.4282992, 0.0023804],
              'M': [-0.0774786, 0.0984413],
              'N': [1.3266132, 0.0092375],
              'P': [1.0860888, 0.0100568],
              'Q': [1.3336109, 0.0111996],
              'R': [1.6492534, 0.0512044],
              'S': [0.7023921, 0.0077661],
              'T': [0.5266550, 0.0311973],
              'U': [-0.0774786, 0.0984413],
              'V': [-0.2447218, 0.0979201],
              'W': [0.2909390, 0.0189282, -0.5479140, 0.0930222, 6.4736619],
              'X': [0.6348884, 0.0180273],
              'Y': [0.6275249, 0.0103896, -0.5744404, 0.0947821, 6.9164963],
              'Z': [1.3761092, 0.0099898]}


def pos_spec_contrib(aa, i, L):
    pos = 9 * (2 * (i - 1) / (L - 1) - 1)
    if aa in 'WY':
        return aaProfiles[aa][0] *\
            math.exp(-1 * aaProfiles[aa][1] * pos**2) +\
            aaProfiles[aa][2] *\
            (math.exp(-1 * aaProfiles[aa][3] * (pos - aaProfiles[aa][4])**2) +
             math.exp(-1 * aaProfiles[aa][3] * (pos + aaProfiles[aa][4])**2))
    else:
        return aaProfiles[aa][0] * math.exp(-1 * aaProfiles[aa][1] * pos**2)


def calc_segment_DG(segment):
    DG_sum = 0.0
    DG_sin_sum = 0.0
    DG_cos_sum = 0.0
    j = 1
    helix_length = len(segment)
    for aa in segment:
        DG = pos_spec_contrib(aa, j, helix_length)
        DG_sum += DG
        DG_sin_sum += (DG * math.sin(100 * (j - 1) * piover180))
        DG_cos_sum += (DG * math.cos(100 * (j - 1) * piover180))
        j += 1
    segment_DG = DG_sum + c0 * \
        (DG_sin_sum**2 + DG_cos_sum**2)**0.5 + c1 + c2 * \
        helix_length + c3 * helix_length**2
    return segment_DG


def running_segment_dg(query_seq, Lmin, Lmax):
    """ Inner running_segment_dg function """
    returnS = ""
    queryLength = len(query_seq)
    for L in range(Lmin, Lmax + 1):
        pre_flank = (L + 1) // 2
        post_flank = L // 2
        returnS += "#L: " + str(L) + "\n"
        returnS += "#Number of sliding windows: " +\
                   str(len(query_seq) - (L - 1)) + "\n"
        for k in range(pre_flank, queryLength - post_flank + 1):
            curr_segment = query_seq[k - pre_flank: k + post_flank]
            helix_DG = calc_segment_DG(curr_segment)
            returnS += str(k)
            returnS += " {0:.3f}".format(helix_DG) + " " + curr_segment + "\n"
    return returnS


def dgCalc(fastafilename, lmin=19, lmax=19):
    """ Calculates dG for a fasta file.
    Inputs - fasta file name, lmin=19, lmax=19"""
    seqs = {}
    returnString = ''
    with open(fastafilename, 'r') as inPipe:
        seq = ''
        pid = ''
        for line in inPipe:
            if line[0] == '>':
                if len(seq) > 0:
                    seqs[pid] = seq
                    seq = ''
                pid = line[1:].strip()
            else:
                seq += line.strip()
        seqs[pid] = seq
    for pid, seq in seqs.items():
        returnString = "#DG calculated from single sequence\n"
        if '|' in pid:
            pid = pid.split('|')[1]
        returnString += "#SeqID: " + pid + "\n"
        returnString += "#>{}".format(pid.strip()) + "\n"
        returnString += "#Seq: " + seq + "\n"
        returnString += "#SeqLength: " + str(len(seq)) + "\n"
        returnString += "#Number of window sizes: " +\
                        str(lmax - lmin + 1) + "\n"
        returnString += running_segment_dg(seq, lmin, lmax)
        returnString += "//"
    return returnString


if __name__ == '__main__':
    assert calc_segment_DG('VTRMVIIMVIAFLICWLPYA') == -2.2780845419344122
    # print(calc_segment_DG('VTRMVIIMVIAFLICWLPYA'))
    # print(calc_segment_DG('TKAMAFIGVSFGITFAIAMVLGPII'))
    # print(calc_segment_DG('TKAMAFIGVSFGITFAIAMVLGPI'))
    # print(calc_segment_DG('TKAMAFIGVSFGITFAIAMVLGP'))
    # print(calc_segment_DG('TKAMAFIGVSFGITFAIAMVLG'))
    inFile = sys.argv[1]
    lmin = 19
    lmax = 19
    if len(sys.argv) > 2:
        lmin = int(sys.argv[2])
        lmax = int(sys.argv[3])
        print(dgCalc(inFile, lmin, lmax))
    else:
        print(dgCalc(inFile))
