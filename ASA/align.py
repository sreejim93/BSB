#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Dynamic Programming assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <SREEJITA MAZUMDER 2702877>
"""



import argparse
import pickle



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = argparse.ArgumentParser(prog = 'python3 align.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Aligns the first two sequences in a specified FASTA\n'
        '  file with a chosen strategy and parameters.\n'
        '\n'
        'defaults:\n'
        '  strategy = global\n'
        '  substitution matrix = pam250\n'
        '  gap penalty = 2')
        
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('output', nargs='*', 
        help='path to an output file where the alignment is saved\n'
             '  (if a second output file is given,\n'
             '   save the score matrix in there)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='print the score matrix and alignment on screen', default=False)
    parser.add_argument('-s', '--strategy', dest='strategy',
        choices=['global','semiglobal','local'], default="global")
    parser.add_argument('-m', '--matrix', dest='substitution_matrix',
        choices=['pam250','blosum62','identity'], default='pam250')
    parser.add_argument('-g', '--gap_penalty', dest='gap_penalty', type=int,
        help='must be a positive integer', default=2)

    args = parser.parse_args()

    args.align_out = args.output[0] if args.output else False
    args.matrix_out = args.output[1] if len(args.output) >= 2 else False
                      # Fancy inline if-else statements. Use cautiously!
                      
    if args.gap_penalty <= 0:
        parser.error('gap penalty must be a positive integer')

    return args



def load_substitution_matrix(name):
    "Loads and returns the specified substitution matrix from a pickle (.pkl) file."
    # Substitution matrices have been prepared as nested dictionaries:
    # the score of substituting A for Z can be found with subst['A']['Z']
    # NOTE: Only works if working directory contains the correct folder and file!
    
    with open('substitution_matrices/%s.pkl' % name, 'rb') as f:
        subst = pickle.load(f)
    return subst
    
    

def load_sequences(filepath):
    "Reads a FASTA file and returns the first two sequences it contains."
    
    seq1 = []
    seq2 = []
    with open(filepath,'r') as f:
        for line in f:
            if line.startswith('>'):
                if not seq1:
                    current_seq = seq1
                elif not seq2:
                    current_seq = seq2
                else:
                    break # Stop if a 3rd sequence is encountered
            else:
                current_seq.append(line.strip())
    
    if not seq2:
        raise Exception('Error: Not enough sequences in specified FASTA file.')
    
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2



def align(seq1, seq2, strategy, substitution_matrix, gap_penalty):
    "Do pairwise alignment using the specified strategy and parameters."
    # This function consists of 3 parts:
    #
    #   1) Initialize a score matrix as a "list of lists" of the appropriate length.
    #      Fill in the correct values for the first row and column given the strategy.
    #        (local / semiglobal = 0  --  global = stacking gap penalties)
    #   2) Fill in the rest of the score matrix using Dynamic Programming, accounting
    #      for the selected alignment strategy, substitution matrix and gap penalty.
    #   3) Perform the correct traceback routine on your filled in score matrix.
    #
    # Both the resulting alignment (sequences with gaps and the corresponding score)
    # and the filled in score matrix are returned as outputs.
    #
    # NOTE: You are strongly encouraged to think about how you can reuse (parts of)
    #       your code between steps 2 and 3 for the different strategies!
    
    
    ### 1: Initialize
    M = len(seq1)+1
    N = len(seq2)+1
    score_matrix = []

    for i in range(M):
        row = []
        score_matrix.append(row)
        for j in range(N):
            row.append(0)

        #####################
        # START CODING HERE #
        #####################

    start_gap_penalty = 0
    if strategy == 'global':
        start_gap_penalty = gap_penalty
        print(start_gap_penalty)
    if strategy == 'semiglobal':
        start_gap_penalty = 0
    score_matrix[0][0] = 0   # Change the zeroes in the first row and column to the correct values.

    for i in range(1,N):
        score_matrix[0][i] = score_matrix[0][i-1] - start_gap_penalty

    for j in range(1,M):
        score_matrix[j][0] = score_matrix[j-1][0] - start_gap_penalty


        #####################
        #  END CODING HERE  #
        #####################

    

    ### 2: Fill in Score Matrix
 
    #####################
    # START CODING HERE #
    #####################
    def dp_function(diagonal, up, left, amino_seq1, amino_seq2,i,j):

        print(diagonal,up,left,amino_seq1, amino_seq2,i,j)
        #Diagonal
        diagonal_value = diagonal + substitution_matrix[amino_seq1][amino_seq2]

        #Up
        upper_value    = up - gap_penalty

        #Left
        left_value     = left - gap_penalty


        if diagonal_value == upper_value and diagonal_value > left_value:
            return diagonal_value, ("diagonal","up")

        elif diagonal_value > upper_value and diagonal_value == left_value:
            return diagonal_value, ("diagonal","left")

        elif diagonal_value > upper_value and diagonal_value > left_value:
            return diagonal_value, "diagonal"

        elif upper_value > diagonal_value and upper_value == left_value:
            return upper_value, ("diagonal","up")

        elif upper_value > diagonal_value and upper_value > left_value:
            return upper_value, "up"

        else:
            return left_value, "left"

    Directions_Track = {}

    for i in range(1,M):
        for j in range(1,N):
            score_matrix[i][j], direction = dp_function(score_matrix[i-1][j-1], score_matrix[i-1][j], score_matrix[i][j-1],seq1[i-1], seq2[j-1], i, j)
            Directions_Track[f'{i},{j}'] = direction


    #####################
    #  END CODING HERE  #
    #####################   

    ### 3: Traceback
    #####################
    # START CODING HERE #
    #####################   
    if strategy == 'global':
        a = M - 1
        b = N - 1
        finalmax=score_matrix[a][b]
        maximum_value=finalmax
        for i in range(a-1,0,-1):
            if score_matrix[i][b]>maximum_value:
                maximum_value=score_matrix[i][b]
                a=i
        L1 = [maximum_value]
        L2 = [str(a) + str(b)]
        x = b - 1
        y = a - 1
        r1 = 0
        c1 = 0

        i = x
        while i > 0:
            j = y
            flag = ''
            max1 = 0
            while j > 0:
                if score_matrix[j][i] > max1:
                    max1 = score_matrix[j][i]
                    r1 = j
                    c1 = i
                    if j != y:
                        flag = str(y) + '-'
                j -= 1
            if flag != '':
                L2.append(flag)
            L1.append(max1)
            y = r1 - 1
            L2.append(str(r1) + str(c1))
            i -= 1
        print(L1,L2)
        aligned_seq1 = ''
        aligned_seq2 = ''
        p, q = 0, 0
        for i in range(len(L2) - 1, -1, -1):
            x = int(L2[i][0])
            y = L2[i][1]
            aligned_seq1 += seq1[x - 1]
            p = x - 1
            if y == '-':
                aligned_seq2 += '-'
            else:
                aligned_seq2 += seq2[int(y) - 1]
                q = int(y) - 1
        print(aligned_seq1, aligned_seq2)
        for i in range(len(seq2[q + 1:])):
            aligned_seq1 += '-'
            aligned_seq2 += seq2[q + 1:][i]
        for i in range(len(seq1[p + 1:])):
            aligned_seq2 += '-'
            aligned_seq1 += seq1[p + 1:][i]
        align_score = finalmax
        print(aligned_seq1, aligned_seq2)
        alignment = (aligned_seq1, aligned_seq2, align_score)
        return (alignment, score_matrix)
    

    if strategy == 'semiglobal':

        maximum_value = 0
        column = len(seq2)
        for row in range(0,M):
            value = score_matrix[row][column]
            if value > maximum_value:
                    maximum_value = value
                    a = row
                    b = column
        finalmax=maximum_value
        for c in range(0,b):
            if(score_matrix[a][c]>maximum_value):
                maximum_value = score_matrix[a][c]
                b=c
        
        L1 = [maximum_value]
        L2 = [str(a) + str(b)]
        x = b - 1
        y = a - 1
        r1 = 0
        c1 = 0

        i = x
        while i > 0:
            j = y
            flag = ''
            max1 = 0
            while j > 0:
                if score_matrix[j][i] > max1:
                    max1 = score_matrix[j][i]
                    r1 = j
                    c1 = i
                    if j != y:
                        flag = str(y) + '-'
                j -= 1
            if flag != '':
                L2.append(flag)
            L1.append(max1)
            y = r1 - 1
            L2.append(str(r1) + str(c1))
            i -= 1
        aligned_seq1 = ''
        aligned_seq2 = ''
        p, q = 0, 0
        for i in range(len(L2) - 1, -1, -1):
            x = int(L2[i][0])
            y = L2[i][1]
            aligned_seq1 += seq1[x - 1]
            p = x - 1
            if y == '-':
                aligned_seq2 += '-'
            else:
                aligned_seq2 += seq2[int(y) - 1]
                q = int(y) - 1
        print(aligned_seq1, aligned_seq2)
        for i in range(len(seq2[q + 1:])):
            aligned_seq1 += '-'
            aligned_seq2 += seq2[q + 1:][i]
        for i in range(len(seq1[p + 1:])):
            aligned_seq2 += '-'
            aligned_seq1 += seq1[p + 1:][i]
        align_score = finalmax
        print(aligned_seq1, aligned_seq2)
        alignment = (aligned_seq1, aligned_seq2, align_score)
        return (alignment, score_matrix)

    if strategy == 'local':
        maximum_value = 0

        for column in range(len(seq2), -1, -1):
            for row in range(M - 1, -1, -1):
                if score_matrix[row][column] < 0:
                    score_matrix[row][column] = 0
                if score_matrix[row][column] > maximum_value:
                    maximum_value = score_matrix[row][column]
                    a = row
                    b = column
        L1 = [maximum_value]
        L2 = [str(a)+str(b)]
        x = b - 1
        y = a - 1
        r1 = 0
        c1 = 0

        i = x
        while i > 0:
            j = y
            flag=''
            max1 = 0
            while j > 0:
                if score_matrix[j][i] > max1:
                    max1 = score_matrix[j][i]
                    r1 = j
                    c1 = i
                    if j!=y:
                        flag=str(y)+'-'
                j -= 1
            if flag!='':
                L2.append(flag)
            L1.append(max1)
            y = r1 - 1
            L2.append(str(r1) + str(c1))
            i -= 1
        aligned_seq1 = ''
        aligned_seq2 = ''
        for i in range(len(L2)-1,-1,-1):
            x = int(L2[i][0])
            y = L2[i][1]
            aligned_seq1+=seq1[x-1]
            if y=='-':
                aligned_seq2+='-'
            else:
                aligned_seq2+=seq2[int(y)-1]
        align_score=maximum_value
        print(aligned_seq1,aligned_seq2)
        alignment = (aligned_seq1, aligned_seq2, align_score)
        return (alignment, score_matrix)

    #####################
    #  END CODING HERE  #
    #####################   






def print_score_matrix(s1,s2,mat):
    "Pretty print function for a score matrix."
    
    # Prepend filler characters to seq1 and seq2
    s1 = '-' + s1
    s2 = ' -' + s2
    
    # Print them around the score matrix, in columns of 5 characters
    print(''.join(['%5s' % aa for aa in s2])) # Convert s2 to a list of length 5 strings, then join it back into a string
    for i,row in enumerate(mat):               # Iterate through the rows of your score matrix (and keep count with 'i').
        vals = ['%5i' % val for val in row]    # Convert this row's scores to a list of strings.
        vals.insert(0,'%5s' % s1[i])           # Add this row's character from s2 to the front of the list
        print(''.join(vals))                   # Join the list elements into a single string, and print the line.



def print_alignment(a):
    "Pretty print function for an alignment (and alignment score)."
    
    # Unpack the alignment tuple
    seq1 = a[0]
    seq2 = a[1]
    score = a[2]
    
    # Check which positions are identical
    match = ''
    for i in range(len(seq1)): # Remember: Aligned sequences have the same length!
        match += '|' if seq1[i] == seq2[i] else ' ' # Fancy inline if-else statement. Use cautiously!
            
    # Concatenate lines into a list, and join them together with newline characters.
    print('\n'.join([seq1,match,seq2,'','Score = %i' % score]))



def save_alignment(a,f):
    "Saves two aligned sequences and their alignment score to a file."
    with open(f,'w') as out:
        out.write(a[0] + '\n') # Aligned sequence 1
        out.write(a[1] + '\n') # Aligned sequence 2
        out.write('Score: %i' % a[2]) # Alignment score


    
def save_score_matrix(m,f):
    "Saves a score matrix to a file in tab-separated format."
    with open(f,'w') as out:
        for row in m:
            vals = [str(val) for val in row]
            out.write('\t'.join(vals)+'\n')
    


def main(args = False):
    # Process arguments and load required data
    if not args: args = parse_args()
    
    sub_mat = load_substitution_matrix(args.substitution_matrix)
    seq1, seq2 = load_sequences(args.fasta)

    # Perform specified alignment
    strat = args.strategy
    gp = args.gap_penalty
    alignment, score_matrix = align(seq1, seq2, strat, sub_mat, gp)

    # If running in "verbose" mode, print additional output
    if args.verbose:
        print_score_matrix(seq1,seq2,score_matrix)
        print('') # Insert a blank line in between
        print_alignment(alignment)
    
    # Save results
    if args.align_out: save_alignment(alignment, args.align_out)
    if args.matrix_out: save_score_matrix(score_matrix, args.matrix_out)



if __name__ == '__main__':
    main()