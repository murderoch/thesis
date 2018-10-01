from itertools import repeat, product, combinations
from collections import Counter

term_letters = 'SPDFGHIKLMNOQRTUVWXYZ'

def get_term_symbols(l, r):
    """Return a list of term symbols for the configuration l^r."""

    # Total number of (ml, ms) pairs for this subshell.
    n = (2*l+1)*2
    # All possible values of ml = -l, -l+1, ..., l-1, l.
    ml = list(range(-l,l+1))
    # All possible values of 2ms = -1, 1. That is, ms = -1/2, +1/2. We work
    # with 2ms instead of ms so that we can handle integers only.
    ms2 = [-1,1]
    # All possible (ml, 2ms) pairs for this subshell.
    ml_ms2 = list(product(ml, ms2))

    # All possible microstates for r electrons in this subshell.
    microstates = list(combinations(range(n), r))
    # The totals ML = sum(ml) and MS2 = sum(2ms) for each microstate
    ML = [sum([ml_ms2[microstate[j]][0] for j in range(r)])
                                    for microstate in microstates]
    MS2 = [sum([ml_ms2[microstate[j]][1] for j in range(r)])
                                    for microstate in microstates]
    # Count the microstates (MS, ML). Store them this way round so we can
    # pick off the ground state term (maximum S) first.
    MS2_ML = Counter(zip(MS2,ML))
    N = len(microstates)

    # Extract the term symbols by starting at the minimum (ML, MS) value and
    # removing microstates corresponding to the (L, S) term it belongs to.
    # Repeat until we're out of microstates.
    terms = []
    while N>0:
        S, L = min(MS2_ML)
        terms.append('{}{}'.format(-S+1, term_letters[-L]))
        for ML in range(L, -L+1):
            for MS in range(S, -S+1,2):
                MS2_ML[MS,ML] -= 1
                if MS2_ML[MS,ML] == 0:
                    del MS2_ML[MS,ML]
                N -= 1
    return terms


print(get_term_symbols(1, 4))