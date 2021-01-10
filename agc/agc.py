#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = ["Maelle Boucher, Matthieu Planté"]
__copyright__ = "Universite CY Tech"
__credits__ = ["Maelle Boucher", "Matthieu Planté"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Maelle Boucher"
__email__ = "bouchermae@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """
    Take compressed fasta file and return sequences yield
    with a length greater or equal to minseqlen
    """
    if amplicon_file.endswith("gz"):
        with gzip.open(amplicon_file, "rb") as filin:
            seq = b""
            for line in filin:
                if line.startswith(b">"):
                    if len(seq) >= minseqlen:
                        yield seq.decode('ascii')
                    seq = b""
                else:
                    seq += line.strip()

        yield seq.decode('ascii')
    else:
        with open(amplicon_file, "r") as filin:
            seq = ""
            for line in filin:
                if line.startswith(">"):
                    if len(seq) >= minseqlen:
                        yield seq
                    seq = ""
                else:
                    seq += line.strip()

            yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Take a compressed fasta file, min length of sequences, and their min count
    return a yield of the sequence and its count
    """
    gen_fasta = read_fasta(amplicon_file, minseqlen)

    seq_count = Counter()
    for seq in gen_fasta:
        seq_count[seq] += 1

    for count in seq_count.most_common():
        if count[1] >= mincount:
            yield count


def get_chunks(sequence, chunk_size):
    """
    Take a sequence and a chunk size
    return a chunk yield
    """
    chunk = []
    for i in range(0, len(sequence), chunk_size):
        if chunk:
            if len(sequence[i:i + chunk_size]) != len(chunk[0]):
                break
        chunk.append(sequence[i:i + chunk_size])

    if len(chunk) >= 4:
        return chunk
    else:
        return None


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    """
    Take a sequence and a kmer_size
    return a kmer yield
    """
    for i in range(0, len(sequence)-kmer_size+1):
        yield sequence[i:i+kmer_size]


def get_identity(alignment_list):
    """
    Take an alignment list of two aligned sequences
    return the identity percentage of them
    """
    count_base = 0
    for i in range(0, len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            count_base += 1
    return count_base / len(alignment_list[0]) * 100


               
def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Retourne générateur séquences non chimérique, Format: yield [seq, count]"""
    subseq_dict = {}
    parents = []
    gen = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    i = 0
    for s, c in gen :
        if i < 2 :
            parents.append(s)
            i+=1
        if i >= 2 :
            chunks = get_chunks(s, chunk_size)
            if chunks != None :
                for chunk1 in get_chunks(parents[0], chunk_size):
                    if chunk1 != None :
                        parent1_kmers = cut_kmer(chunk1, kmer_size)
                        for chunk2 in get_chunks(parents[1], chunk_size):
                             if chunk2 != None :
                                 parent2_kmers = cut_kmer(chunk2, kmer_size)
                                 for chunk in chunks :
                                     if len(subseq_dict) < 9 :
                                         if (chunk in parent1_kmers) or (chunk in parent2_kmers) :
                                             subseq_dict[s] +=1
                                                                      
    seq_parents = []
    for subseq in subseq_dict:
        if len(seq_parents<=2):
            if subseq in parents[0] and subseq in parents[1]:
                seq_parents.append(subseq)
        else:
            break

    identity_list = []
    parents_chunks = [get_chunks(parent, chunk_size) for parent in parents]
    for elm in seq_parents:
        for i in range(len(get_chunks(elm, chunk_size))):
            alignment_list = nw.global_align(get_chunks(elm[i], chunk_size), parents_chunks[i], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identityc1 = get_identity(alignment_list)
            alignment_list = nw.global_align(get_chunks(elm[i], chunk_size), parents_chunks[i], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identityc2 = get_identity(alignment_list)
            identity_list.append([identityc1,identityc2])
    for elm in dereplication_fulllength(amplicon_file,minseqlen,mincount):
        yield [elm[0], elm[1]]
        
        
def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Take a compressed fasta file, the minimum length for sequences, their min occurences,
    the chunk size and the kmer size
    return the coccurence of each non chimera sequences
    """
    return [elm for elm in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)]


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """
    Take OTU list and write it in a fasta file
    """
    with open(output_file, "w") as filout:
        for i in range(0,len(OTU_list)):
            filout.write(">OTU_{} occurrence:{}\n{}\n".format(i+1,
                OTU_list[i][1], fill(OTU_list[i][0])))

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    


if __name__ == '__main__':
    main()