#name: Scatter Plot Python
#language: python
#tags: demo, viewers, hide-suggestions
#input: dataframe t
#input: column ColumnName
#output: graphics
import matplotlib.pyplot as plt
import numpy as np
def computeLPSArray(pat, M, lps):
    len = 0 # length of the previous longest prefix suffix

    lps[0] # lps[0] is always 0
    i = 1

    # the loop calculates lps[i] for i = 1 to M-1
    while i < M:
        if pat[i]== pat[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            # This is tricky. Consider the example.
            # AAACAAAA and i = 7. The idea is similar
            # to search step.
            if len != 0:
                len = lps[len-1]

                # Also, note that we do not increment i here
            else:
                lps[i] = 0
                i += 1
# Python program for KMP Algorithm
def KMPSearch(pat, txt):
    count = 0
    M = len(pat)
    N = len(txt)

    # create lps[] that will hold the longest prefix suffix
    # values for pattern
    lps = [0]*M
    j = 0 # index for pat[]

    # Preprocess the pattern (calculate lps[] array)
    computeLPSArray(pat, M, lps)

    i = 0 # index for txt[]
    while i < N:
        if pat[j] == txt[i]:
            i += 1
            j += 1

        if j == M:
            count+=1
            j = lps[j-1]

        # mismatch after j matches
        elif i < N and pat[j] != txt[i]:
            # Do not match lps[0..lps[j-1]] characters,
            # they will match anyway
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
    return count


def count_subsequense_in_sequnse(sequencesColumn, subsequence):
    count = 0
    for i in range(len(sequencesColumn)):
       count += KMPSearch(subsequence, sequencesColumn[i])
    return count
import numpy as np
import matplotlib.pyplot as plt

sequencesColumn = t[ColumnName].values
counts = np.zeros(64)
x = np.arange(64)
numberSequence = 0
for char1 in list(['g', 'a', 't', 'c']):
    for char2 in list(['g', 'a', 't', 'c']):
        for char3 in list(['g', 'a', 't', 'c']):
            # create string for char 1, char 2, char 3
            subsequence = char1 + char2 + char3
            # count number of subsequence in sequences
            count = count_subsequense_in_sequnse(sequencesColumn, subsequence)
            counts[numberSequence] = count
            numberSequence += 1


fig, ax = plt.subplots()
ax.hist(x, bins= np.arange(1,65), weights=counts)
plt.show()
