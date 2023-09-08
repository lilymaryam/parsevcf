import os
import argparse
import numpy as np
import ete3
import tqdm as progressbar

parser = argparse.ArgumentParser()
parser.add_argument('tree', type=str,help='path to MAT')
parser.add_argument('-s', '--samples', required=False, type=str,help='comma separated list of samples')
parser.add_argument('-stdout', action='store_true', help='print matrix to stdout instead of a file')

args = parser.parse_args()
tree = args.tree
t = ete3.Tree(tree, format=1)
if args.samples:
    samps = args.samples.split(',')
else:
    samps = sorted([leaf.name for leaf in t])

def path_to_root(ete_tree, node):
    # Browse the tree from a specific leaf to the root
    node = ete_tree.search_nodes(name=node)[0]
    path = [node]
    #print(node)
    while node:
        node = node.up
        path.append(node)
    #print(path)
    return path


def dist_matrix(tree, samples):
    samp_ancs = {}
    samp_dist = {}
    #t = ete3.Tree(tree, format=1)
    #for each input sample, find path to root and branch lengths
    for s in progressbar.tqdm(samples, desc="Finding roots and branch lengths"):
        s_ancs = path_to_root(tree, s)
        samp_ancs[s] = s_ancs
    #create matrix for samples
    #matrix = np.empty([len(samples),len(samples)])
    matrix = np.full((len(samples),len(samples)), -1)

    #print(matrix)
    for i in progressbar.trange(len(samples), desc="Creating matrix"): # trange is a tqdm optimized version of range
        s = samples[i]

        for j in range(len(samples)):
            os = samples[j]
            #Future goal: add catch to prevent reiteration of already checked pairs 
            if os == s:
                matrix[i][j] = '0'
                #print(matrix)
            if os != s:
                if matrix[i][j] == -1:
                    #find lca, add up branch lengths
                    s_path = 0 
                    os_path = 0 
                    for a in samp_ancs[s]:
                        
                        s_path += a.dist
                        if a in samp_ancs[os]:
                            
                            lca = a
                            s_path -= a.dist
                            #print(t.get_node(a))
                            break
                    for a in samp_ancs[os]:
                        
                        os_path += a.dist
                        if a == lca:
                            #print('os_path', os_path)
                            os_path -= a.dist
                            break
                    #print(s,os)
                    #print(s_path, os_path)
                    matrix[i][j] = int(s_path + os_path)
                    matrix[j][i] = int(s_path + os_path)
    return samples, matrix                    

samps, mat = dist_matrix(t, samps)

'''
for i in range(len(mat)):
    for j in range(len(mat[i])):
        if mat[i][j] != mat[j][i]:
            print(i,j)
'''


#print(f'sample\t'+'\t'.join(samps))
for i in range(len(samps)):
    #strng = np.array2string(mat[i], separator='\t')[1:-1]
    line = [ str(int(count)) for count in mat[i]]
    if args.stdout:
        print(f'{samps[i]}\t' + '\t'.join(line))
    else:
        with open(f"{os.path.basename(tree)}distance_matrix.tsv", "w") as outfile:
            outfile.write(f'{samps[i]}\t' + '\t'.join(line))
    