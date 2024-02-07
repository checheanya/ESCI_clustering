'''
Prunning algorithm for phylogenetic tree implementation
'''

from Node import Leaf
from Node import Branch
from Node import InternalNode
from math import log


def main():
    transition_probability = {"A":{"G":0.04,"C":0.04,"T":0.04,"A":0.88}, "C":{"G":0.04,"C":0.88,"T":0.04,"A":0.04}, "G":{"G":0.88,"C":0.04,"T":0.04,"A":0.04}, "T":{"G":0.04,"C":0.04,"T":0.88,"A":0.04}}
    s1 = Leaf("C", transition_probability)
    s2 = Leaf("C", transition_probability)
    n = InternalNode(transition_probability)
    n.set_node_at_i(Branch.LEFT,s1)
    n.set_node_at_i(Branch.RIGHT,s2)
    l = n.compute_likelihood_by_nucleotide()
    p = 0
    for k in l:
        p += l[k]*0.25
    print(log(p))


if __name__ == '__main__':
    main()
