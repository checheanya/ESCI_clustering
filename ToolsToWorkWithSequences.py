from Evolution import Evolution as evol 
from Sequence import Sequence as seq

class ToolsToWorkWithSequences:
    def __init__( self ):
        pass
    
    def nucleotide_statistics( sequence ) :
        # verify that sequence is an object of type Sequence, raise excpetion if not
        if not isinstance( sequence , seq ):
            raise TypeError(f"Sequence must be of type Sequence, instead got type {type(sequence)}")
       
        # initialize counts for each nucleotide and start 
        counts = { n : 0 for n in ("A", "T", "C", "G") }

        # iterate across sequence and increment the count associated to each nucleotide
        i = 0
        while i < sequence.sequence_length():
            counts[ sequence.nucleotide_at_position( i ) ] += 1
            i += 1
       
        # scale all counts by the sequence length to convert count to relative frequency, then return
        # dictionary of frequencies with key/value pairs corresponding to nucleotide/frequency-of-nucleotide  
        counts = { n : counts[ n ] / sequence.sequence_length() for n in counts }
        return( counts )

    def observed_pairwise_nucleotide_distance( seq1, seq2 ):
        # inizialize counter of differences
        count_diff = 0

        # for each nucleotide pair at each position, if they differ, increment counter by one
        for n1, n2 in zip( seq1, seq2 ):
            if n1 != n2:
                count_diff += 1

        # return counter
        return( count_diff )

def main():
    ''' Create a method called nucleotide_statistics that uses as parameter
    an object of type sequence and returns a dictionary with the
    percentage of A, C, T, G found in the sequence. Apply the method to
    each of the evolved sequences '''
    
    # here we recreate the diagram in slide 9 of 'Practical Session 3_GOOD_Francesc.pptx'
#
#             S_a <--------  A  --------> S_c                 # 0'th generation
#              |                           |                  
#              v                           | 
#      S_a <--S_a--> S_b                   |                  # 100'th generation 
#       |             |                    |  
#       |             |                    v  
#       |             |            S_c <--S_c--> S_d          # 150'th generation
#       |             |             |             |                 
#       v             v             v             v                 
#   Species_a    Species_b     Species_c     Species_d       # 300'th generation     
    
    # create sequence object of ancestral sequence and associated transition matrix
    S_a = seq( "S_a", "AAATACTATATTTCGCTCCGCGTTCGAGACGAATAAC" )
    transition_probability = {
                              "A":{"G":0.03,"C":0.03,"T":0.03,"A":0.91},
                              "C":{"G":0.03,"C":0.91,"T":0.03,"A":0.03}, 
                              "G":{"G":0.91,"C":0.03,"T":0.03,"A":0.03}, 
                              "T":{"G":0.03,"C":0.03,"T":0.91,"A":0.03}
                              }

    # initialize the evolution class with ancestral sequence S_a
    evolution = evol( S_a, transition_probability )
    
    # split the 0th generation and evolve until first speciation event for 100 generations
    evolution.split_species_in_two( "S_a", "S_c" )
    evolution.evolve( 100 )
    
    # split 100'th generation of S_a into S_a and S_b to accumulate mutations in new species
    evolution.split_species_in_two( "S_a", "S_b" )
    # evolve until second speciation event
    evolution.evolve( 50 )

    # split 150'th generation of S_c into S_c and S_d to accumulate mutations in new species
    evolution.split_species_in_two( "S_c", "S_d" )
    # evolve until the 300'th generation 
    evolution.evolve( 150 )

    # total time equates to 300 generations ( 100 + 50 + 150 )
    # now we can analyze the four species within the evolution object through the nucleotide_statistics method
    for species_name in sorted( evolution.get_list_of_species_name() ):
        # species_seq is the sequence of the specified species (S_a, S_b, S_c, S_d in our case),
        species_seq = evolution.get_sequence_species( species_name ) 
        
        # stats is a dictionary corresponding to the percentage of A, C, T, G found in the sequence
        stats = ToolsToWorkWithSequences.nucleotide_statistics( species_seq )

        print( f"Species { species_name } is characterized by nucleotide frequencies of:" )
        print( ", ".join( [ f"{nuc} : {stats[ nuc ]:.3f}" for nuc in stats ] ), end = "\n\n" )


    ''' Create a method called observed_pairwise_nucleotide_distance that
    takes as parameters two sequences. Returns the number of
    nucleotides that for the same position are different in the two
    sequences '''
    
    # define two sequences to test with
    test_sequence_1 = "AAATACTATATTTCGCTCCGCGTTCGAGACGAATAAC"
    test_sequence_2 = "AATTACTATCGTTCGCTCCACGTTTGAGTAGAATATG"
   
    # testing the observer_pairwise_nucleotide_distance method
    print( f"Sequence 1 {test_sequence_1}" )
    print( f"Sequence 2 {test_sequence_2}" )
    print( "Number of nucleotides that differ in the same position is: ", end = "" ) 
    print( ToolsToWorkWithSequences.observed_pairwise_nucleotide_distance( test_sequence_1, test_sequence_2 ), end = "\n\n" )

    # testing nucleotide_statistics method with first test sequence
    print( f"Sequence 1 {test_sequence_1}" )
    print( ToolsToWorkWithSequences.nucleotide_statistics( seq( "test1", test_sequence_1 ) ) )

if __name__ == "__main__":
    main()
