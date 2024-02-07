'''
Tree-related classes for the prunning algorithm. 
'''

from enum import Enum
from abc import ABC,abstractmethod


class Branch(Enum):
    LEFT = 0
    RIGHT = 1

    
class TransitionMatrixOfBranch(object):    
        def __init__(self, branch, transition_matrix):
            self.branch = branch    
            self.transition_matrix = transition_matrix
        
        def get_transition_matrix(self):
            return self.transition_matrix
        
        def get_branch(self):
            return self.branch    
    

class Node(ABC):    
    def __init__(self, transition_matrix):
        self.transition_matrix = transition_matrix    
    
    @abstractmethod
    def compute_likelihood_by_nucleotide(self):
        pass  
    
    def set_node_at_i(self, branch, node):
        if isinstance(branch, TransitionMatrixOfBranch):
            self.descendent_nodes[Branch[branch.get_branch()]] = node
        
    def get_transition_matrix(self):
        return self.transition_matrix
        

class InternalNode(Node):    
    def __init__(self, transition_matrix):
        '''
        Constructor
        '''
        self.descendent_nodes = [None]*2
        self.transition_matrix = transition_matrix
     
    def set_node_at_i(self, i, node):
        if isinstance(i, Branch):
            self.descendent_nodes[i.value] = node

    def compute_likelihood_by_nucleotide(self):
        # initialize the probabilities
        probability_each_nucleotide = {"A":1,"C":1,"T":1,"G":1} 
        for n in self.descendent_nodes:
            current_prob_descent = n.compute_likelihood_by_nucleotide()
            transition_matrix_n = n.get_transition_matrix() 
            for nuc in probability_each_nucleotide: 
                p = 0 
                for nec in current_prob_descent:
                    p += transition_matrix_n[nuc][nec]*current_prob_descent[nec]
                probability_each_nucleotide[nuc] *= p
        return probability_each_nucleotide 


class Leaf(Node):
    
    def __init__(self,nucleotide, transition_matrix):
        self.set_nucleotide(nucleotide)
        self.transition_matrix = transition_matrix    
    
    def set_nucleotide(self, nucleotide):
        self.probability_each_nucleotide = {"A":0,"C":0,"T":0,"G":0}
        self.probability_each_nucleotide[nucleotide] = 1       
    
    def compute_likelihood_by_nucleotide(self):
        return self.probability_each_nucleotide    
    
