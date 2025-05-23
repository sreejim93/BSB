#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Hidden Markov Models assignment in the Algorithms in Sequence Analysis course at the VU.

INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    <Sreejita Mazumder, Stud. nr.: 2702877, VUnetID: smr373>
"""

import os.path as op

from os import makedirs
from math import log10
from hmm_utility import parse_args, load_fasta, load_tsv, print_trellis, print_params, serialize



def viterbi(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the most probable state path, the corresponding P(X), and trellis."""

    allStates = A.keys()
    
    emittingStates = E.keys()
    L = len(X) + 2

    # Initialize
    V = {k:[0] * L for k in allStates} # The Viterbi trellis
    V['B'][0] = 1.

    # Middle columns
    for i,s in enumerate(X):
        for l in emittingStates:
            terms = [V[k][i] * A[k][l] for k in allStates]
            V[l][i+1] = max(terms) * E[l][s]

    # Last column
    for k in allStates:
        term = V[k][i+1] * A[k]['E'] 
        if term > V['E'][-1]:
            V['E'][-1] = term
            pi = k # Last state of the State Path

    # FOR VITERBI ONLY: Trace back the State Path
    l = pi
    i = L-2
    while i:
        i -= 1
        for k in emittingStates:
            if V[k][i] * A[k][l] * E[l][X[i]] == V[l][i+1]:
                pi = k + pi
                l = k
                break

    P = V['E'][-1] # The Viterbi probability: P(X,pi|A,E)
    return(pi,P,V) # Return the state path, Viterbi probability, and Viterbi trellis

def forward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Forward probability and corresponding trellis."""

    allStates = A.keys()
    L = len(X) + 2
    
    # Initialize
    F = {k:[0] * L for k in allStates}
    F['B'][0] = 1

    #####################
    # START CODING HERE #
    #####################
    # HINT: The Viterbi and Forward algorithm are very similar! 
    # Adapt the viterbi() function to account for the differences.

    # Middle columns
    # for ...

    # Last columns
    # for ...:
    #     F['E'][-1] += ...
    emittingStates = E.keys()
   

    def add_forward(i,l,s):
        terms = [F[k][i] * A[k][l] for k in allStates]
        F[l][i+1] = sum(terms) * E[l][s]
        # Last columns
        term = [F[k][i+1] * A[k]['E'] for k in allStates]
        F['E'][-1] = sum(term)
    
    for i,s in enumerate(X):
        for l in emittingStates:
            add_forward(i,l,s)

    
    #####################
    #  END CODING HERE  #
    #####################

    P = F['E'][-1] # The Forward probability: P(X|A,E)
    return(P,F)

def backward(X,A,E):
    """Given a single sequence, with Transition and Emission probabilities,
    return the Backward probability and corresponding trellis."""

    allStates = A.keys()
    emittingStates = E.keys()
    L = len(X) + 2 ## L = 12+2 = 14

    # Initialize
    B = {k:[0] * L for k in allStates} # The Backward trellis
    for k in allStates:
        B[k][-2] = A[k]['E']

    #####################
    # START CODING HERE #
    #####################
    # Remaining columns
    # for i in range(L-3,-1,-1):
    #     s = seq[i]
    #     ...
    

    def add_backward(i,k,s):
        terms = [A[k][l]*E[l][s]*B[l][i+1] for l in emittingStates]
        B[k][i] = sum(terms)

    for i in range(L-3, -1, -1):
        seq=X
        s = seq[i]
        for k in allStates:
            add_backward(i,k,s)
            
    #####################
    #  END CODING HERE  #
    #####################

    P = B['B'][0] # The Backward probability -- should be identical to Forward!
    return(P,B)

def baumwelch(set_X,A,E):
    """Given a set of sequences X and priors A and E,
    return the Sum Log Likelihood of X given the priors,
    along with the calculated posteriors for A and E."""

    allStates = A.keys()
    print("All states")
    print(allStates)
    emittingStates = E.keys()
    print("emitting states")
    print(emittingStates)
    
    # Initialize a new (posterior) Transition and Emission matrix
    new_A = {}
    for k in A:
        new_A[k] = {l:0 for l in A[k]}
    print("new matrix a")
    print(new_A)

    new_E = {}
    for k in E:
        new_E[k] = {s:0 for s in E[k]}
    print("new matrix E")
    print(new_E)
    # Iterate through all sequences in X
    SLL = 0 # Sum Log-Likelihood
    for X in set_X:
        P,F = forward(X,A,E)  # Save both the forward probability and the forward trellis
        _,B = backward(X,A,E) # Forward P == Backward P, so only save the backward trellis
        SLL += log10(P)
        print("SLL")
        print(SLL)
        #####################
        # START CODING HERE #
        #####################
        # Inside the for loop: Expectation
        # Count how often you observe each transition and emission.
        # Add the counts to your posterior matrices. (new_A, new_E)
        # Remember to normalize to the sequence's probability P!

        #transitions
        #for k in allStates:
            #transitions for the last state E
            #new_A[k]['E'] += (F[k][-2] * A[k]['E'] / P)
               #other transitions
        #for l in emittingStates:
            #for i in range(len(X)):
               # new_A[k][l] += F[k][i] * A[k][l] * E[l][X[i]] * B[l][i+1] / P
        #(lambda k : new_A[k]['E'] += F[k][-2] * A[k]['E'] / P for k in allStates)(k)
        #(lambda words, rules: sum([[word[:-len(rule)]] if word.endswith(rule) else [] for word in words for rule in rules], []))(str_test.split(), stem_rules)
        
        def newa(n):
           new_A[k]['E'] += (F[k][-2] * A[k]['E'] / P)
           for l in emittingStates:
               for i in range(len(X)):
                   new_A[k][l] += F[k][i] * A[k][l] * E[l][X[i]] * B[l][i+1] / P
        for k in allStates:
           newa(k)
            
        #emissions
        '''
            for k in emittingStates:
                #new_E[k][s] += F[k][i+1]3 * B[k][i+1] / P
        
        
        #print("new matrix ####### E")
        #print(new_E)
        '''
        '''
        index_seq=[]
        symbols_seq=[]
        for i,s in enumerate(X):
            index_seq.append(i)
            symbols_seq.append(s)
            
        states_k=[]
        for k in emittingStates:
            states_k.append(k)
            
        
        for s in X:
            for k in emittingStates:
                a = (lambda i,s,k: (F[k][i+1] * B[k][i+1] / P),states_k,index_seq, symbols_seq)
                print("tuple")
                print(a)
                new_E[k][s] += (a) 
    '''
        def newE(i,s):
            for k in emittingStates:
                new_E[k][s] += F[k][i+1] * B[k][i+1] / P
            
        for i,s in enumerate(X):
            newE(i,s)

    # Outside the for loop: Maximization
    # Normalize row sums to 1 (except for one row in the Transition matrix!)
    # new_A = ...
    # new_E = ...
    '''
    a=[]
    for l in emittingStates:
        a.append(l)
        print("emmitingStates:" + l)
        sumOfValues = list(map(lambda l: sum(new_E[l].values()),a))'''
        
    def sum_all():
        return [sum(new_E[l].values()) for l in emittingStates]
    sumOfValues=sum_all()
    
    
    
    #print("sumOfValues:" )
    #print(sumOfValues[0])
   # print(sumOfValues[1])
    
    
    #apply if else
    for l in emittingStates:
        if l=='L':
            for emission, prob in new_E[l].items():
                new_E[l][emission] = (prob / sumOfValues[0])
            else:
                new_E[l][emission] = (prob / sumOfValues[1])
            
    '''
    for l in emittingStates:
        for emission, prob in new_E[l].items():
            print("Emission")
            print(emission)
            print("prob")
            print(prob)
            new_E[l][emission] = (prob / sumOfValues[0])
        print(new_E)
'''
    for k in allStates:
        sumOfValues = sum(new_A[k].values())
        if sumOfValues > 0:
            for transition, prob in new_A[k].items():
                new_A[k][transition] = (prob / sumOfValues)

    #####################
    #  END CODING HERE  #
    #####################
    
    return(SLL,new_A,new_E)
    

def main(args = False):
    "Perform the specified algorithm, for a given set of sequences and parameters."
    
    # Process arguments and load specified files
    if not args: args = parse_args()

    cmd = args.command            # viterbi, forward, backward or baumwelch
    verbosity = args.verbosity
    set_X, labels = load_fasta(args.fasta)  # List of sequences, list of labels
    A = load_tsv(args.transition) # Nested Q -> Q dictionary
    E = load_tsv(args.emission)   # Nested Q -> S dictionary
    
    def save(filename, contents):
        if args.out_dir:
            makedirs(args.out_dir, exist_ok=True) # Make sure the output directory exists.
            path = op.join(args.out_dir,filename)
            with open(path,'w') as f: f.write(contents)
        # Note this function does nothing if no out_dir is specified!



    # VITERBI
    if cmd == 'viterbi':
        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the most probable state path, with the corresponding probability and matrix
            Q, P, T = viterbi(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.path' % label, Q)
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            print('>%s\n Path = %s' % (label,Q))
            if verbosity: print(' Seq  = %s\n P    = %1.2e\n' % (X,P))
            if verbosity >= 2: print_trellis(T, X)
            


    # FORWARD or BACKWARD
    elif cmd in ['forward','backward']:
        if cmd == 'forward':
            algorithm = forward
        elif cmd == 'backward':
            algorithm = backward

        for j,X in enumerate(set_X): # For every sequence:
            # Calculate the Forward/Backward probability and corresponding matrix
            P, T = algorithm(X,A,E)

            # Save and/or print relevant output
            label = labels[j]
            save('%s.matrix' % label, serialize(T,X))
            save('%s.p' % label, '%1.2e' % P)
            if verbosity >= 2:
                print('\n>%s\n P = %1.2e\n' % (label,P))
                print_trellis(T, X)
            elif verbosity: print('>%-10s\tP = %1.2e' % (label,P))



    # BAUM-WELCH TRAINING
    elif cmd == 'baumwelch':
        # Initialize
        i = 1
        i_max = args.max_iter
        threshold = args.conv_thresh

        current_SLL, A, E = baumwelch(set_X,A,E)
        if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
        if verbosity >= 2: print_params(A,E)
        
        last_SLL = current_SLL - threshold - 1 # Iterate at least once

        # Iterate until convergence or limit
        while i < i_max and current_SLL - last_SLL > threshold:
            i += 1
            last_SLL = current_SLL

            # Calculate the Sum Log-Likelihood of X given A and E,
            # and update the estimates (posteriors) for A and E.
            current_SLL, A, E = baumwelch(set_X,A,E)

            if verbosity: print('Iteration %i, prior SLL = %1.2e' % (i,current_SLL))
            if verbosity >= 2: print_params(A,E)

        converged = current_SLL - last_SLL <= threshold
        final_SLL = sum([log10(forward(X,A,E)[0]) if forward(X,A,E)[0] > 0 else 0 for X in set_X])

        # Save and/or print relevant output
        save('SLL','%1.2e\t%i\t%s' % (final_SLL, i, converged))
        save('posterior_A',serialize(A))
        save('posterior_E',serialize(E))
        if verbosity: print('========================================\n')

        if converged:
            print('Converged after %i iterations.' % i)
        else:
            print('Failed to converge after %i iterations.' % i_max)

        if verbosity:
            print('Final SLL: %1.2e' % final_SLL)
            print('Final parameters:')
            print_params(A,E)



if __name__ == '__main__':
	main()
