import os
import sys

exec(open('sequence_functions.py').read())

# Dictionary of restriction site sequences
restriction_sites = {'MfeI': 'CAATTG', 'AccI':  'GTATAC', 'NotI': 'GCGGCCGC', 'SnaBI': 'TACGTA', 'SpeI': 'ACTAGT' }
# These are constant sequences that do not include restriction sites or chain sequences.
constant_regions = ['tctcgagta'] # 0
constant_regions.append('catgggctccaggctgctctgttgggtgctgctttgtctcctgggagcaggcccagtaaaggctgga') # 1
constant_regions.append('aagatc') # 2
constant_regions.append('acgtgacaccacccaaagtctcactgtttgagcctagcaaggcagaaattgccaacaagcagaaggccaccctggtgtgcctggcaagagggttctttccagatcacgtggagctgtcctggtgggtcaacggcaaagaagtgcattctggggtctgcaccgacccccaggcttacaaggagagtaattactcatattgtctgtcaagccggctgagagtgtccgccacattctggcacaaccctaggaatcatttccgctgccaggtccagtttcacggcctgagtgaggaagataaatggccagaggggtcacctaagccagtgacacagaacatcagcgcagaagcctggggacgagcagactgtggcattactagcgcctcctatcatcagggcgtgctgagcgccactatcctgtacgagattctgctgggaaaggccaccctgtatgctgtgctggtctccggcctggtgctgatggccatggtcaagaaaaagaactctgggagtggagccacaaatttctctctgctgaaacaggctggagatgtggaggaaaaccccggccctatgaagagcctgcgcgtgctgctggtcatcctgtggctg') # 3
constant_regions.append('tcgtgggtctggagccaa') # 4
constant_regions.append('gacattcagaacccggaaccggct') # 5
constant_regions.append('cagctgaaggacccccgatctcaggatagtactctgtgcctgttcaccgactttgatagtcagatcaatgtgcctaaaaccatggaatccggaacttttattaccgacaagtgcgtgctggatatgaaagccatggacagtaagtcaaacggcgccatcgcttggagcaatcagacatccttcacttgccaggatatcttcaaggagaccaacgcaacatacccatcctctgacgtgccctgtgatgccaccctgacagagaagtctttcgaaacagacatgaacctgaattttcagaatctgagcgtgatgggcctgagaatcctgctgctgaaggtcgctgggtttaatctgctgatgacactgcggctgtggtcctcatgaattcggaccgtgtccaatgtagc') # 6
constant_regions.append('gtcgacaatcaacctctgga') # 7

#### Create a lentiviral insert by concatening constant regions (including restriction sites) and A and B chain sequences
# Trim the full_nt_sequence provided by Mixcr so that the ultimate sequence is in the correct reading frame.
def trim_by_remainder( seq, factor = 3, remainder=0):
    r = len(seq) % factor
    trimlen = r - remainder
    if trimlen < 0:
        trimlen += 3
    if trimlen == 0:
        return seq
    return seq[:-trimlen]

# Generate the insert sequence given an alpha and beta full_nt_sequence
def generate_insert( achain, bchain ):
    # Read in the a and b chain sequence
    # The a chain sequence needs to have a remainder of 1 when divided by 3 to keep the insert in frame.
    achain = trim_by_remainder( achain, 3, 0 )

    # The b chain sequence needs to have a remainder of 0 when divided by 3 to keep the insert in frame.
    bchain = trim_by_remainder( bchain, 3, 1 )

    # The insert is built in this particular order to make sure we can keep track of
    # The cumulative lengths of the sequence between restriction sites so we know where precisely the
    # restriction sites should be. 
    expected_restriction_sites = {}
    insert_nuc = constant_regions[0]
    expected_restriction_sites['NotI'] = len(insert_nuc)
    insert_nuc += restriction_sites['NotI'] + constant_regions[1] + bchain + constant_regions[2]
    expected_restriction_sites['SnaBI'] = len(insert_nuc)
    insert_nuc += restriction_sites['SnaBI'] + constant_regions[3] 
    expected_restriction_sites['MfeI'] = len(insert_nuc)
    insert_nuc += restriction_sites['MfeI'] + constant_regions[4] + achain + constant_regions[5] 
    expected_restriction_sites['AccI'] = len(insert_nuc)
    insert_nuc += restriction_sites['AccI'] + constant_regions[6] 
    expected_restriction_sites['SpeI'] = len(insert_nuc)
    insert_nuc += restriction_sites['SpeI'] + constant_regions[7]


    # Check if the sequence aaggctgga is duplicated, resulting in a repeat of KAGKAG
    # when it should just be KAG. If so, replace it so that it only happens once.
    # A duplication happened in the past, which may mean that mixcr gave the KAG in this instance
    # whereas it didn't in previous instances.
    insert_nuc = insert_nuc.replace('aaggctggaaaggcta','aaggctgga')

    insert_aa = translate(insert_nuc,0)
    print (expected_restriction_sites)

    # Check nucleotide sequence s for errant restriction sites
    # and then fix them without changing AA sequence. 
    # Assumes s is in reading frame 0

    #rs = restriction_sites.values()[0]
    #insert_nuc = insert_nuc[:30] + rs + insert_nuc[ (30+len(rs)):]
    while True:
        print ('\nChecking restriction sites \n')
        # Go through the restriction sites and check if there are any in incorrect places
        correct = True # number of restriction sites whose positioning might be correct
        rsite_locations = find_restriction_sites( restriction_sites, insert_nuc ) 
        for s in expected_restriction_sites:
            print ('Expected:', s, expected_restriction_sites[s])
            print ('Actual:', s, rsite_locations[s])
            if not s in rsite_locations or rsite_locations[s] == []:
                print ('Expected site ' + s + ' not present.')
                break
            for i in rsite_locations[s]:
                if i != expected_restriction_sites[s]:
                    correct = False
                    print (s + ' incorrect')
                    insert_nuc = swap_codons( insert_nuc, s, i )
                else:
                    print (s + ' correct')
        if correct:
            break

    if len( insert_aa.split('_') ) > 2 :
        print ('ERROR: EXTRA STOP CODON')
        insert_aa = insert_aa.replace('_','****_****')
    else:
        print ('Correct number of stop codons')

    if True:
        # Check if there are errant stop codons
        # and then fix them if there are
        print(insert_aa)
        stopsite = insert_aa.index('IRTVSNVALVSTINLW') - 1 # Location of expected stop codon
        stopsites = find_subsequence_instances( insert_aa, '_' )
        print ('Stopsites', stopsites, stopsite)

        insert_nuc = 'a' + insert_nuc

    print ('\nInsert nucleotide sequence:')
    print (insert_nuc)
    print ('\nInsert amino acid sequence:')
    print (insert_aa)

    return( insert_nuc, insert_aa)
