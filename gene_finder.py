# -*- coding: utf-8 -*-
"""
Finds the protein sequence generated from a sequence of
nucleotides.

@author: Emily Lepert

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ##


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'

    We need to have two extra unit tests so that all nucleotides are covered
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == "A":
        return("T")
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "G":
        return "C"


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        Current unit tests are sufficient
        because they are of different lengths and use all letters
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    newDna = ""
    k = 0
    for i in dna:
        newDna = get_complement(i) + newDna
        k += 1
    return(newDna)


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string. (TAG, TGA, TAA)

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        We should test all of the stop codon sequences as well as
        making sure that the program knows not to read any random 3
        sequence of letters as stop codons, but only groups of 3 as codons
    >>> rest_of_ORF("ATGAAATAA")
    'ATGAAA'
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # keeping track of the sequence
    ORF = ""
    i = 0
    # for every nucleotide in the sequence
    while i < len(dna):
        # update the sequence with the new nucleotide
        ORF = ORF + dna[i]
        # when we have a new codon, analyze it to see if it's a stop codon
        if (i % 3) == 0 and (dna[i:i+3] == "TAG" or dna[i:i+3] == "TGA"
                             or dna[i:i+3] == "TAA"):
            # bc ORF was keeping track of all of the sequence, remove the
            # stop codon from it
            return(ORF[0:i])
        i += 1
    return(ORF)


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        To make sure the program works I need to add an example of a non-nested
        ORF
    >>> find_all_ORFs_oneframe("ATGCATGAAATGTGTAGATAGTGCCCATAA")
    ['ATGCATGAAATGTGTAGA']
    >>> find_all_ORFs_oneframe("ATGCATGAAATGTGTAGATAGTGCATGCCATAA")
    ['ATGCATGAAATGTGTAGA', 'ATGCCA']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("GATGCATGAATGTAGATAGATATGTGCCC")
    ['ATG', 'ATGTGCCC']
    """
    # index
    i = 0
    # list containing all ORFs
    ORF = []
    # the current ORF
    currentORF = ""
    while i < len(dna):
        # if we're at the end of a codon and it is a start codon
        if i % 3 == 0 and dna[i:i+3] == "ATG":
            # the current ORF is the rest of the DNA or until
            # a stop codon is found
            currentORF = rest_of_ORF(dna[i:])
            # if there isn't a stop codon in the sequence, then
            # the ORF is the rest of the DNA
            if currentORF is None:
                currentORF = dna[i:]
                i = len(dna)
            else:
                i += len(currentORF)
            # add the ORF strand to the list
            ORF.append(currentORF)
        else:
            i += 1
    return(ORF)


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        We need to check that the function only finds non-nest ORFs

    >>> find_all_ORFs("ATGCATATGGAATAGTAG")
    ['ATGCATATGGAA']
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    i = 0
    # store all ORFs found
    allORFs = []
    while i < 3:
        ORF = find_all_ORFs_oneframe(dna[i:])
        for j in ORF:
            allORFs.append(j)
        i += 1
    return(allORFs)


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        This test case makes sure that there are only non-nest
        ORFs for both strands
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    allORF = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return(allORF)


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFs = find_all_ORFs_both_strands(dna)
    # keeps track of the length of the longest strand
    length = 0
    # keeps track of the current longest strand
    longest = ''
    for i in ORFs:
        if len(i) >= length:
            longest = i
            length = len(i)
    return(longest)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        """
    i = 0
    # keeps track of the longest ORF string
    longest = 0
    # each loop is a new trial
    while i < num_trials:
        # shuffle the string
        shuffle = shuffle_string(dna)
        # get the length of the longest string
        current_longest = len(longest_ORF(shuffle))
        # update the longest string length if necessary
        if current_longest > longest:
            longest = current_longest
        i += 1
    return(longest)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        We can add more tests that test for more of the codons
        including stop codons
        >>> coding_strand_to_AA("ATGTAAATGCATTGCGGT")
        'M|MHCG'
        >>> coding_strand_to_AA("ATGCCACGTACTGAGTTCATC")
        'MPRTEFI'
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    i = 0
    protein = ''
    while i < len(dna):
        # if we're at the start of a codon and there is a full codon
        if i % 3 == 0 and len(dna[i:i+3]) == 3:
            # translate the codon to an amino acid
            amino_acid = aa_table[dna[i:i+3]]
            protein += amino_acid
        i += 1
    return(protein)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # find the longest ORF over 1500 trials
    threshold = longest_ORF_noncoding(dna, 1500)
    # find all the ORFs in this strand
    ORF = find_all_ORFs_both_strands(dna)
    amino_sequences = []
    for i in ORF:
        if len(i) >= threshold:
            amino_sequences.append(coding_strand_to_AA(i))
    return(amino_sequences)


dna = load_seq("./data/X73525.fa")
print(gene_finder(dna))


if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    #doctest.run_docstring_examples(coding_strand_to_AA, globals())
