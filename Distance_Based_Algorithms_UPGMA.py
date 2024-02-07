'''
/ Draft code /
I've noticed quite bizarre behaviour where the newick tree
may not necessarily reflect the true relationship between species
but i've attributed this issue primarily to the small length
of sequence and unoptimized rate of mutation.

Thus I have modified the frequency of mutation (made them more rare)
AND increased the length of sequence. Thus, a less frequent mutation 
rate and longer sequences more closely resembles a situtation origin
of evolution our sequences.

These modifications are quite illustrative of UPGMA's pitfalls, since 
it is only a good modeller of phylogenetic relationships when the 
molecular clock hypothesis corresponds to the sequences in question 
(typically longer sequences) )

'''

from ToolsToWorkWithSequences import ToolsToWorkWithSequences as tools
from DistanceMatrix import DistanceMatrix as dist_mat
from Evolution import Evolution as evol 
from Sequence import Sequence as seq

class DistanceBasedAlgorithms(object):
    
    def __init__(self, d_matrix):
        # include distance matrix in self for access in methods 
        self.d_matrix = d_matrix
        # initialize newick tree as a dictionary
        self.newick_tree = dict()
        self.largest_cluster = None

  
    def count_subclusters ( self, clusters ):
        if isinstance( clusters, str ):
            return( 1 ) 
        return( sum( self.count_subclusters( clu ) for clu in clusters ) )

  
    def UPGMA(self):
        for i in range( self.d_matrix.n_rows() ):
            for j in range( self.d_matrix.n_rows() ):
                print( self.d_matrix.get_value_i_j( i, j ), end = " , " )
            print()

        # if less than 2 rows remain, base case has been reached and no more clusters may be collapsed
        if self.d_matrix.n_rows() < 2:
            tree = str(self.Newick( self.largest_cluster ) )
            tree = tree.replace( "',", "':" )
            tree = tree.replace( "],", "]:" )

            print( tree )
            return()
        
        # find smallest distance
        n = self.d_matrix.n_rows() 
        s_d, s_i, s_j = float( "inf" ), None, None
        for i in range( n ):
            for j in range( n ):
                if i == j:
                    continue

                d_ij = self.d_matrix.get_value_i_j( i, j )
                if d_ij < s_d:
                    s_d, s_i, s_j = d_ij, i, j
        
        # introduce cluster into newick tree 
        new_d = s_d / 2
        
        species_i = self.d_matrix.get_name_of_species()[ s_i ]
        species_j = self.d_matrix.get_name_of_species()[ s_j ]
        new_name = ( ( species_i , species_j) )
        
        self.newick_tree[ new_name ] = new_d 
        if isinstance( species_i, str ):
            self.newick_tree[ species_i ] = new_d 
        if isinstance( species_j, str ):
            self.newick_tree[ species_j ] = new_d 

        self.largest_cluster = new_name

        # create a copy of the matrix, with appropriate clusters
        aux_matrix = self.d_matrix.copy()
        aux_matrix.change_name_species_i( s_i, new_name ) 
        aux_matrix.remove_species_i( s_j )
        
        # calculate the cardinality of our considered clusters 
        i_cardinality = self.count_subclusters( species_i )
        j_cardinality = self.count_subclusters( species_j )

        # calculate and introduce the averaged pairwise distance between the new cluster and the 
        # other clusters/species
        for j in range( n - 1 ):
            if j == s_j:
                continue
            new_value = self.d_matrix.get_value_i_j( j, s_j ) * j_cardinality + self.d_matrix.get_value_i_j( s_i, j ) * i_cardinality
            new_value /= i_cardinality + j_cardinality 
            aux_matrix.add_value_i_j( j, s_i, new_value )
        
        # overwrite d_matrix with auixliary matrix
        self.d_matrix = aux_matrix
        
        # recursively call itself until it reaches base case of one cluster left
        self.UPGMA()

  
    def Newick( self, ancestor ):
        # recursively find subclusters witihn clusters and append members and distance 
        if isinstance( ancestor, str ):
            return( ancestor ) #, self.newick_tree[ ancestor ] ] )

        current_tree = []
        for child in ancestor:
            current_tree += [ self.Newick( child ), self.newick_tree[ ancestor ] ] 

        return( current_tree )


# this function returns a list of four sequences of type Sequence
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
    Species_a    Species_b     Species_c     Species_d       # 300'th generation     '''    

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
    # call UPGMA to generate a Newick Formatted tree
    original_matrix = DistanceBasedAlgorithms( distance_matrix )
    original_matrix.UPGMA() 


if __name__ == "__main__":
    main()
  
