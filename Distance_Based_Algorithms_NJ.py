'''
/ Draft code / 

Note that i've modified the sequence (elongating it) and reduced
the probability of mutations.

I could not produce a newick tree in time, I am naively overwriting
the names of the nodes whenever i redefine the self.m matrix.

At a conceptual level I would have arbitrarily selected 
some inner node, defined it as root, and proceeded to generate 
the phylogenetic relationships by utilizing the Distances to other
neighbours as calculated by the neighour_joining method. 

'''

from numpy import matrix as mat
from ToolsToWorkWithSequences import ToolsToWorkWithSequences as tools
from DistanceMatrix import DistanceMatrix as dist_mat
from Evolution import Evolution as evol 
from Sequence import Sequence as seq

class DistanceBasedAlgorithms(object):
    def __init__(self, d_matrix):
        # include distance matrix in self for access in methods 
        self.m = []
        for i in range( d_matrix.n_rows() ):
            self.m.append([])
            for j in range( d_matrix.n_rows() ):
                self.m[ -1 ].append( d_matrix.get_value_i_j( i, j ) )
        self.tree = ""


    def neighbour_joining(self):
        # base case where two neighbours can be included into the tree
        # without
        if len(self.m) <= 2:
            #print( self.tree )
            return()

        # Compute row_sums
        total_distances = []
        for i in range( len( self.m ) ):
            row_sum = 0
            for j in range( len( self.m ) ):
                row_sum += self.m[ i ][ j ] 
            total_distances.append( row_sum )
        
        # Find minimum element in D*ij
        smallest_neighbours = [ float( "inf" ) , [ None, None ] ]
        d_aux = []
        for i in range( len( self.m ) ):
            d_aux.append( [] )
            for j in range( len( self.m ) ):
                if j == i:
                    d_aux[ -1 ].append( 0 )
                    continue
                d_aux[ -1 ].append( ( len( self.m ) - 2 ) * self.m[ i ][ j ] - total_distances[ i ] - total_distances[ j ] )
                if d_aux[ -1 ][ -1 ] < smallest_neighbours[ 0 ]:
                    smallest_neighbours[ 1 ][ 0 ] = i
                    smallest_neighbours[ 1 ][ 1 ] = j
                    smallest_neighbours[ 0 ] = d_aux[ -1 ][ -1 ]

        # Compute distance from i,j
        delta_ij = ( total_distances[ smallest_neighbours[ 1 ][ 0 ] ] - total_distances[ smallest_neighbours[ 1 ][ 1 ] ] ) / ( len( self.m ) - 2  )
        
        # Define limblength/distance between nodes
        si, sj = smallest_neighbours[ 1 ][ 0 ], smallest_neighbours[ 1 ][ 1 ]
        limb_i = 0.5 * ( self.m[ si ][ sj ] + delta_ij )
        limb_j = 0.5 * ( self.m[ si ][ sj ] - delta_ij )
        
        print(f" Limb lengths of our currently considred neighbours are: {limb_i}, {limb_j} " )
        self.tree += "(" + str(si) + ":" + str(limb_i) + "," + str(sj) + ":" + str(limb_j) + ")"
        new_node_idx = len(self.m)
        # form matrix D' by removing i-th and j-th row (in our case constructing a new matrix and populating it)
        d_prime = [ [ 0 for _ in range( len( self.m )- 1 ) ] for _ in range( len( self.m )  - 1) ]
        
        # UNDESIRED BEHAVIOUR; WILL FIX
        for a in range( len( self.m ) ):
            for k in range(  len( self.m ) ):
                if k == si or k == sj or a == k:
                    continue
                d_prime[ a - 1 ][ k - 1 ] = ( self.m[si][k] + self.m[k][sj] - self.m[si][sj] ) / 2
       
        self.m = d_prime
        m = mat( d_prime ) 
        print( m )
        # recursively call function until only to neighbours are present
        self.neighbour_joining()
                

# this function returns a list of four objects of type Sequence
def simulate_four_sequences():
    ''' 
    Recreate the diagram in slide 9 of "Practical Session 3_GOOD_Francesc.pptx"
 
              S_a <--------  A  --------> S_c                 # 0'th generation
               |                           |                  
               v                           | 
       S_a <--S_a--> S_b                   |                  # 100'th generation 
        |             |                    |  
        |             |                    v  
        |             |            S_c <--S_c--> S_d          # 150'th generation
        |             |             |             |                 
        v             v             v             v                 
    Species_a    Species_b     Species_c     Species_d       # 300'th generation     
    '''    

    # create sequence object of ancestral sequence and associated transition matrix
    S_a = seq( "S_a", "AAAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATAAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAATACTATATTTCGCTCCGCGTTCGAGACGAATAACACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAAATACTATATTTCGCTCCGCGTTCGAGACGAATAACAATACTATATTTCGCTCCGCGTTCGAGACGAATAAC" )
    transition_probability = {
                              "A":{"G":0.000004,"C":0.000004,"T":0.000004,"A":0.999988},
                              "C":{"G":0.000004,"C":0.999988,"T":0.000004,"A":0.000004}, 
                              "G":{"G":0.999988,"C":0.000004,"T":0.000004,"A":0.000004}, 
                              "T":{"G":0.000004,"C":0.000004,"T":0.999988,"A":0.000004}
                              }

    evolution = evol( S_a, transition_probability )
    evolution.split_species_in_two( "S_a", "S_c" )
    evolution.evolve( 100 )
    evolution.split_species_in_two( "S_a", "S_b" )
    evolution.evolve( 50 )
    evolution.split_species_in_two( "S_c", "S_d" )
    evolution.evolve( 150 )
    return( evolution )

def main():
    # generate four sequences according to Session 3's tree (see function)
    four_sequences = simulate_four_sequences()
   
    # retreive the names of the four sequences and populate the distance matrix
    # with 0's corresponding to the distance between all possible species/sequence pairs
    names = sorted( four_sequences.get_list_of_species_name() )
    distance_matrix = dist_mat( [ name for name in names ] )
   
    # populate each i, j pair (i.e each possible sequence pair amongst our four generated
    # sequences ) with the corresponding distance between them utilizing the
    # observed_pairwise_nucleotide_distance method to count the number of 
    # differing nucleotides at the same position between pairs
    i = 0
    while i < len( names ):
        seq_i = four_sequences.get_sequence_species( names[ i ] ).__str__()
        seq_i = ( seq_i[ seq_i.find( "[" ) + 1 : -1 ] ).strip()
        seq_i = seq_i.replace( "'", "" )
        seq_i = seq_i.replace( ", ", "" )

        j = i + 1
        while j < len( names ):
            seq_j = four_sequences.get_sequence_species( names[ j ] ).__str__()
            seq_j = ( seq_j[ seq_j.find( "[" ) + 1 : -1 ] ).strip()
            seq_j = seq_j.replace( "'", "" )
            seq_j = seq_j.replace( ", ", "" )
            
            # calculate distance between seq_i and seq_j 
            distance_ij = tools.observed_pairwise_nucleotide_distance( seq_i, seq_j )
            
            # introduce the distance between seq_i and seq_j into the distance matrix
            distance_matrix.add_value_i_j( i, j, distance_ij )
            #print( distance_matrix.get_value_i_j( i, j ))
            
            j += 1
        i += 1

    # initialize a DistanceBasedAlgorithms object with the distance_matrix and
    # call NJ to generate a Newick Formatted tree
    original_matrix = DistanceBasedAlgorithms( distance_matrix )
    original_matrix.neighbour_joining() 

if __name__ == "__main__":
    main()
