import os
import sys

exec(open("functions.py").read())
exec(open("generate_lentiviral_insert.py").read())

## To run this, open python and from python navigate to the folder that contains this file
## 0. Make your input file name. And then change the inputfile = XXXX to your input file name and outputfile = YYYYYY to your output file name. Put the input file in the GenerateInsertScripts folder.
inputfile = "Treg4_seq.txt" # Replace the sys.argv[1] with the input text file name (tab-delimited)
outputfile = "Treg_output_seq.txt" # Replace the sys.argv[2] with the output text file name (doesn't exist yet, will be created). It will be comma separated (csv)
## 1. Open python. Type "import os"
## 2. Type this: os.chdir('C:\\Users\\rtewari\\Desktop\\GenerateInsertScripts\\GenerateInsertScripts\\')
## 3. Type this: exec(open("lentiviral_insert_job_windows.py").read())

# Takes as input the 'full_nt_sequence' for both the alpha and beta junctions that were output
# by mixcr. Outputs a lentiviral insert sequence for that TCR pair. Can also use nucleotide sequences 
# provided by other pipelines (10X and Rhapsody TCRs have been used successfully)

# Takes as input a tab-delimited text file where each row is a TCR pair, and the columns are
# sample (name for the insert) | description (optional) | cdr3_a_nt ( alpha full_nt_sequence as provided by Mixcr) | cdr3_b_nt  (beta full_nt_sequence as provided by Mixcr) | cdr3_a (just the alpha junction amino acid sequence) | cdr3_b (just the beta junction amino acid sequence)

pairs = read_tab_delimited( inputfile )

print( pairs.keys() )

# Here specify the output file
fout = open(outputfile,'w')
# Here specify the information for the TCR pairs. The information is in the form of
# an arbitrary sample name, full nucleotide sequences for the alpha and beta chains, 
# as well as corresponding amino acid sequences (just as a check to make sure the
# nucleotide sequences were processed correctly.
# All of these things are in the form of arrays to allow for the generation of
# multiple inserts for multiple TCR pairs.
samples = pairs['sample']  # sample names for the inserts
achains = pairs['cdr3_a_nt']
bchains = pairs['cdr3_b_nt']
ajunctions = pairs['cdr3_a']
bjunctions= pairs['cdr3_b']

# Loop through the input sequences and generate an insert for each entry
# Insert sequence info is outputted to a file
# Other information is printed out
i = 0
success = 0
failures = 0
for (achain,bchain,aj,bj) in zip(achains,bchains,ajunctions,bjunctions):
    print ( 'Row ', str(i) )
    print (achain,bchain,aj,bj)
    insert = generate_insert( achain, bchain )
    if aj in insert[1] and bj in insert[1]:
        print ('AC:', achain )
        print ('Insert:', insert)
        print ('AJ:', aj)
        print ('BJ:', bj)
        print ('The right one')
        success = success + 1
    else:
        print ('Not the right one')


    print ('Successes: ' + str(success) + '/' + str(len(bchains)) )
    fout.write( str(i+1) + '. alpha junction: ' + aj + ', '+ 'beta junction: ' + bj + '\n')
    s = str(samples[i]) + '\n'
    if 'description' in pairs:
        s = s + pairs['description'][i] + '\n'
    fout.write( s )
    fout.write('Insert nucleotide sequence: ' + insert[0] + '\n' )
    fout.write('Insert amino acid sequence: ' + insert[1] + '\n\n' )
    if len(insert[1].split('_')) > 2:
        failures = failures+1
    i = i + 1
fout.close()

print (fout.name)

print ('Num failures:', failures)
