#!/usr/bin/env python3
import sys

# salt_bridge_cutoff = 4  # distance in Ångström to consider a salt bridge

def parseLine(line):
    num = int(line[6:11])
    atom = line[13:16].strip()
    alt_loc = line[16].strip()
    residue = line[17:20]
    chain = line[21]
    chainnum = int(line[22:26]) - 1  # Change to 0 based 
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    # print(num, chainnum, x, y, z)
    return (x, y, z, num, atom, residue, chain, chainnum, alt_loc)


def calcEuclidianDist(x1, y1, z1, x2, y2, z2):
    return ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5


def calcSaltBridges(filename, chain, salt_bridge_cutoff = 4, con_req=2):
    raw_bridges = {}
    bridges = []
    otherBridges = []
    currChain = []
    otherChains = []
    POS = ['ARG', 'LYS', 'HIS']
    NEG = ['ASP', 'GLU']
    with open(filename, 'r') as pdbFile:
        # chain = 'C'
        multi_model = False
        for line in pdbFile:
            #### Only use the first model in case of multiple
            if line.startswith("MODEL"):
                if multi_model:
                    break
                else:
                    multi_model = True
            if line[:4] == 'ATOM'\
               and line[17:20] in ['ASP', 'GLU', 'ARG', 'LYS', 'HIS']\
               and line[13:15] in ['OD', 'OE', 'NH', 'NZ', 'NE', 'ND']:
                if line[21] == chain:
                    currChain.append(parseLine(line))
                    # print(line)
                    # print(parseLine(line))
                else:
                    otherChains.append(parseLine(line))

    for i, res in enumerate(currChain):
        for secRes in currChain[i + 1:]:
            # Check ALL pairings for salt bridges
            # if secRes[7] - res[7] > 7:
            #     break
            # If it is the same residue OR the first and second residues are NOT opposite charged
            # continue to the next
            if res[7] == secRes[7]\
               or not ((res[5] in POS and secRes[5] in NEG)
               or (res[5] in NEG and secRes[5] in POS)):
                continue
            dist = 10
            # Check what connections
            if res[5] in POS:
                if res[4] in ["NH1", "NH2", "NZ", "ND1", "NE2"] and secRes[4] in ["OD1", "OD2", "OE1", "OE2"]:
                    dist = calcEuclidianDist(res[0], res[1], res[2],
                                             secRes[0], secRes[1], secRes[2])
            elif secRes[5] in POS:
                if secRes[4] in ["NH1", "NH2", "NZ", "ND1", "NE2"] and res[4] in ["OD1", "OD2", "OE1", "OE2"]:
                    dist = calcEuclidianDist(secRes[0], secRes[1], secRes[2],
                                             res[0], res[1], res[2])
            else:
                print("Something went wrong with distance calculations")
                print(res, secRes, sep='\n')
                sys.exit()
            # key = res[6] + str(res[7]) + '-' + secRes[6] + str(secRes[7])
            if dist < salt_bridge_cutoff:  # and key not in bridges:
                # bridges[key] = [res, secRes]
                bridge = [res[4:], secRes[4:], dist]
                # print(bridge)
                # sys.exit()
                # bridges.append(bridge)
                b_id = str(res[7]) + "->" + str(secRes[7])
                if b_id in raw_bridges:
                    raw_bridges[b_id].append(bridge)
                else:
                    raw_bridges[b_id] = [bridge]
                # bridges.append(bridge)
                # print(bridge)
                # print("Hit (dist ", dist, ") between: ")
                # print(res)
                # print(secRes)
        # for secRes in otherChains:
        #     if res[7] == secRes[7]\
        #        and not ((res[5] in POS and secRes[5] in NEG)
        #        or (res[5] in NEG and secRes[5] in POS)):
        #         continue
        #     dist = calcEuclidianDist(res[0], res[1], res[2],
        #                              secRes[0], secRes[1], secRes[2])
        #     # key = res[6] + str(res[7]) + '-' + secRes[6] + str(secRes[7])
        #     if dist < salt_bridge_cutoff:  # and key not in bridges:
        #         # bridges[key] = [res, secRes]
        #         bridge = [res[4:], secRes[4:], dist]
        #         otherBridges.append(bridge)

    # bridges.extend(otherBridges)
    # print(raw_bridges)
    for b_id, brids in raw_bridges.items():
        if len(brids) >= con_req:
            bridges.extend(brids)
    return bridges
    # for each in bridges:
    # print(each)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Error: Takes filename and chain as argument")
        sys.exit()
    filename = sys.argv[1]
    chain = sys.argv[2]
    for b in calcSaltBridges(filename, chain, 4, 1):
        print(b)
# print(len(currChain))
# print(len(otherChains))
# print(calcEuclidianDist(0, 0, 0, 1, 1, 1))
