import os
import sys

# Assume the dictionary is a dictionary of lists
# that each have the same length
def write_tab_delimited( d, filename, sep='\t', header=True ):
    keys = d.keys()
    keys.sort()
    fout = open(filename,'w')
    fout.write( str(keys[0]) )
    for k in keys[1:]:
        fout.write(sep+str(k))
    fout.write('\n')
    N = len( d.values()[0] )
    for i in range(N):
        fout.write( str(d[keys[0]][i]) )
        for k in keys[1:]:
            fout.write(sep+str(d[k][i]))
        fout.write('\n')
    fout.close()

# Read in a .csv file that has the desired a and b sequences, then write out the clonal information
def read_tab_delimited( f, header=True, sep='\t' ): # returns a dictionary
    fin = open(f,'r')
    if header:	
        header = fin.readline()
        if header == '':
            return {}
        if header[-1] == '\n':
            header = header[:-1]
        header = header.split(sep)
        data = { h:[] for h in header }
        lh = len(header)
        for line in fin:
            if line[-1] == '\n':
                line = line[:-1]
            tokens = line.split(sep)
            if lh != len(tokens):
                print ('Mismatch!', lh, len(tokens))
                print (header)
                print (tokens)
            for i in range(0,len(tokens)):
                try: 
                    t = float(tokens[i])
                    if int(t) == t:
                        t = int(t)
                    data[header[i]].append(t)
                except ValueError:
                    data[header[i]].append(tokens[i])				
                except OverflowError:
                    data[header[i]].append(tokens[i])
    else:
        data = [ ]	
        for line in fin:
            if line == '':
                return {}
            tokens = line.rstrip().split()
            for i in range(len(tokens)):
                if len(data) <= i:
                    data.append([])
                try: 
                    t = float(tokens[i])
                    if int(t) == t:
                        t = int(t)
                    data[i].append(t)
                except ValueError:
		                    data[i].append(tokens[i])				
                except OverflowError:
		                    data[i].append(tokens[i])				
	
        fin.close()
    return data


# d is a dictionary of lists
# sort each list in order according keytosort
def sort_dict( d, keytosort = 'pid' ):
    order = range(len(d[keytosort]))
    order.sort( key = lambda i: d[keytosort][i] )
    for key in d.keys():
        d[key] = [ d[key][i] for i in order ]
    return(order)
