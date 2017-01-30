#!/usr/bin/python


from __future__ import print_function
import argparse
import sys
import re
#import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import random

# License

'''
MIT License

Copyright (c) 2017 Johan Bengtsson-Palme

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

# Usage
'''
This little script fragments DNA sequences and outputs sequences with simulated errors.
Please type 'python seqsim.py -h' for usage and options!
'''

def main():

    ## Read arguments
    args = read_arguments()
    #print(args)

    if args.target_length == 0 and args.length_dropoff == 0:
        sys.exit('Target length (-l) and dropoff (-d) cannot both be zero!')

    ## Set up substitution matrix
    if (args.matrixfile == ''):
        subm = build_substitution_matrix(args.sub_rate,args.ins_rate,args.del_rate,args.hp_rate,args.qd_rate)
    else:
        subm = load_substitution_matrix(args.matrixfile)

    if (args.matrixexportfile != ""):
        exportsubm(subm,args.matrixexportfile)
    
    
    ## Read input file (genome data)
    sequence_data = read_input(args.infile)

    fragments = {}

    eprint("Generating fragments...\n")
    
    ## Select a fragment
    for i in range(0,args.nreads):
        if ((i+1) / args.nreads * 10 == int((i+1) / args.nreads * 10)):
            eprint(str((i+1) / args.nreads * 100) + "% ... ")
        m = len(sequence_data.keys())
        #r = random.randint(0, m-1)
        #keys = sequence_data.keys()
        #seqID = keys.get(r)
        seqID = random.choice(list(sequence_data.keys()))
        seq = sequence_data[seqID]

        ## Derive length parameters
        l = args.target_length
        if (args.length_distribution > 0):
            l = int(random.gauss(args.target_length,args.length_distribution))
            ml = l
        if args.target_length == 0:
            ml = len(seq)
        else:
            ml = l
        if (args.length_dropoff > 0):
            for b in range(0,ml):
                if (random.random() < args.length_dropoff):
                    break
            l = b
        
        ## Find start point (must be low enough to contain full length)
        ms = len(seq) - l
        
        if ms < 0:
            ms = 0
        s = random.randint(0,ms)

        if (s+l) > len(seq):
            l = len(seq) - s
            
        seqFragment = seq[s:(s+l)]
        fragments[seqID + "." + str(i+1)] = seqFragment
        
    ## Output original (to file in the future)
    fo = open(args.outfile+".fragments.fasta", 'w')
    for fragID in fragments.keys():
        print(">"+fragID+" "+str(len(fragments[fragID]))+" bp", file = fo)
        print(fragments[fragID], file = fo)

    eprint("\nIntroducing errors according to model...\n")
    ## Introduce errors
    processed, qualities = model_errors(fragments, subm, args.target_length, args.length_distribution)

    ## Output processed sequence
    if (args.outformat == "fasta"):
        fp = open(args.outfile+".processed.fasta", 'w')
        for fragID in processed.keys():
            print(">"+fragID+" "+str(len(processed[fragID]))+" bp", file = fp)
            print(processed[fragID], file = fp)
    else:
        fp = open(args.outfile+".processed.fastq", 'w')
        for fragID in processed.keys():
            print("@"+fragID+" "+str(len(processed[fragID]))+" bp", file = fp)
            print(processed[fragID], file = fp)
            print("+", file = fp)
            qstring = encode_quality(qualities[fragID], args.outformat)
            print(qstring, file = fp)

    eprint("Finished!\n")
    
def read_arguments():
    ## Needed args: target length, length distr, read dropoff, substition, insert, deletion, homeopolymer, substitution matrix
    
    parser = argparse.ArgumentParser(description='Fragment DNA sequence(s) and output sequence with simulated errors.')

    parser.add_argument('-i','--input', dest='infile', action='store', type=str, default='', help='input file in FASTA format', required = True)
    parser.add_argument('-o','--output', dest='outfile', action='store', type=str, default='', help='prefix for the output files', required = True)
    parser.add_argument('-n','--output_reads', dest='nreads', action='store', type=int, default=1, help='number of simulated reads to output')
    parser.add_argument('-f','--output_format', dest='outformat', action='store', type=str, default="fasta", help='sets the output format. Can be "fasta", or FASTQ in "sanger" (offset 33), "solexa" (offset 59), "illumina64" (offset 64) or "illumina33" (offset 33)')

        
    parser.add_argument('-l','--target_length', dest='target_length', action='store', type=int, default=0, help='the target length of the output reads')
    parser.add_argument('--dist','--length_distribution', dest='length_distribution', action='store', type=float, default=0.0, help='the distribution around the target length of the output reads')
    parser.add_argument('-d','--length_dropoff', dest='length_dropoff', action='store', type=float, default=0.0, help='the probability for read dropoff')
    parser.add_argument('-s','--substitution_rate', dest='sub_rate', action='store', type=float, default=0.0, help='sets the global substitution rate')
    parser.add_argument('-z','--insertion_rate', dest='ins_rate', action='store', type=float, default=0.0, help='sets the global insertion rate')
    parser.add_argument('-x','--deletion_rate', dest='del_rate', action='store', type=float, default=0.0, help='sets the global deletion rate')
    parser.add_argument('-p','--homeopolymer_rate', dest='hp_rate', action='store',type=float, default=0.0, help='sets the homeopolymer problem rate')
    parser.add_argument('-q','--quality_dropoff', dest='qd_rate', action='store', type=float, default=0.0, help='sets the quality dropoff per base')

    parser.add_argument('-m','--substitution_matrix', dest='matrixfile', action='store', type=str, default='', help='allows setting specific error rates from a matrix file')

    parser.add_argument('--export_substitution_matrix', dest='matrixexportfile', action='store', type=str, default='', help='will export the substitution matrix generated to a file')

    args = parser.parse_args()
    return args

    ## END OF FUNCTION

def read_input(infile):
    sequence_data = {}
    for seq_record in SeqIO.parse(infile, "fasta"):
        sequence_data[seq_record.id] = seq_record.seq
    return sequence_data

    ## END OF FUNCTION

def build_substitution_matrix(sub_rate,ins_rate,del_rate,hp_rate,qd_rate):
    ## Format X->A,X->T,X->G,X->C,X->insA,X->insT,X->insG,X->insC,X->del,X->hp
    sm = {}
    sm["A"] = [0,sub_rate/3,sub_rate/3,sub_rate/3,ins_rate/4,ins_rate/4,ins_rate/4,ins_rate/4,del_rate,hp_rate,qd_rate]
    sm["T"] = [sub_rate/3,0,sub_rate/3,sub_rate/3,ins_rate/4,ins_rate/4,ins_rate/4,ins_rate/4,del_rate,hp_rate,qd_rate]
    sm["G"] = [sub_rate/3,sub_rate/3,0,sub_rate/3,ins_rate/4,ins_rate/4,ins_rate/4,ins_rate/4,del_rate,hp_rate,qd_rate]
    sm["C"] = [sub_rate/3,sub_rate/3,sub_rate/3,0,ins_rate/4,ins_rate/4,ins_rate/4,ins_rate/4,del_rate,hp_rate,qd_rate]
    return sm
    ## END OF FUNCTION

def load_substitution_matrix(matrixfile):
    ## Format X->A,X->T,X->G,X->C,X->insA,X->insT,X->insG,X->insC,X->del,X->hp
    sm = {}
    f = open(matrixfile, 'r')
    strvals = str.split(f.readline())
    subvals = [float(i) for i in strvals]
    sm["A"] = subvals
    strvals = str.split(f.readline())
    subvals = [float(i) for i in strvals]
    sm["T"] = subvals
    strvals = str.split(f.readline())
    subvals = [float(i) for i in strvals]
    sm["G"] = subvals
    strvals = str.split(f.readline())
    subvals = [float(i) for i in strvals]
    sm["C"] = subvals
    return sm
    ## END OF FUNCTION

def model_errors(fragments, subm, targetlen, lendistr):
    processed = {}
    qualities = {}
    nfrags = len(fragments.keys())
    i = 0
    Q = 1
    for fragID in fragments.keys():
        i = i + 1
        if (i / nfrags * 10 == int(i / nfrags * 10)):
            eprint(str(i / nfrags * 100) + "% ... ")
        procseq = []
        quality = []
        for bp in range(0,len(fragments[fragID])):
            nt = fragments[fragID][bp]
            Q = Q * (1 - subm[nt][10])
            r = random.random() * Q
            if (r <= subm[nt][0]):
                procseq.append("A")
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1]):
                procseq.append("T")
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2]):
                procseq.append("G")
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2] + subm[nt][3]):
                procseq.append("C")
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2] + subm[nt][3] + subm[nt][4]):
                procseq.append(nt)
                procseq.append("A")
                quality.append(random.gauss(Q,0.2))
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2] + subm[nt][3] + subm[nt][4] + subm[nt][5]):
                procseq.append(nt)
                procseq.append("T")
                quality.append(random.gauss(Q,0.2))
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2] + subm[nt][3] + subm[nt][4] + subm[nt][5] + subm[nt][6]):
                procseq.append(nt)
                procseq.append("G")
                quality.append(random.gauss(Q,0.2))
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2] + subm[nt][3] + subm[nt][4] + subm[nt][5] + subm[nt][6] + subm[nt][7]):
                procseq.append(nt)
                procseq.append("C")
                quality.append(random.gauss(Q,0.2))
                quality.append(random.random() * Q)
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2] + subm[nt][3] + subm[nt][4] + subm[nt][5] + subm[nt][6] + subm[nt][7] + subm[nt][8]):
                ## DO NOT APPEND ANYTHING, THIS IS A DELETION
                pass
            elif (r <= subm[nt][0] + subm[nt][1] + subm[nt][2] + subm[nt][3] + subm[nt][4] + subm[nt][5] + subm[nt][6] + subm[nt][7] + subm[nt][8] + subm[nt][9]):
                ntrep = 0
                for fw in range(bp+1,len(fragments[fragID])):
                    fwnt = fragments[fragID][fw]
                    if fwnt == nt:
                        ntrep = ntrep + 1
                    else:
                        break
                hr = random.random()
                if (hr < subm[nt][9] * ntrep):
                    if (random.random() > 0.5):
                        procseq.append(nt+nt)
                        quality.append(random.gauss(Q,0.2))
                        quality.append(random.gauss(Q,0.2))
                    else:
                        pass
                else:
                    procseq.append(nt)
                    quality.append(random.gauss(Q,0.2))
            else:
                procseq.append(nt)
                quality.append(random.gauss(Q,0.2))

        finalseq = ''.join(procseq)
        if (lendistr == 0):
            if (len(finalseq) > targetlen):
                finalseq = finalseq[0:targetlen]
                quality = quality[0:targetlen]
            while (len(finalseq) < targetlen):
                finalseq = finalseq + "N"
                quality.append(0)
                
        processed[fragID+".p"] = finalseq
        qualities[fragID+".p"] = quality

    eprint("\n")
    return processed, qualities
    ## END OF FUNCTION

def exportsubm(sm,matrixfile):
    fm = open(matrixfile, "w")
    p = re.compile('[^0-9. ]')
    print(p.sub("",str(sm["A"])), file = fm)
    print(p.sub("",str(sm["T"])), file = fm)
    print(p.sub("",str(sm["G"])), file = fm)
    print(p.sub("",str(sm["C"])), file = fm)
    ## END OF FUNCTION

def encode_quality(quality, format):
    qstring = ""
    maxQ = 40
    minQ = 3
    offset = 33
    if (format == "sanger"):
        maxQ = 60
        minQ = 0
        offset = 33
    if (format == "solexa"):
        maxQ = 40
        minQ = -5
        offset = 59
    if (format == "illumina64"):
        maxQ = 40
        minQ = 3
        offset = 64
    if (format == "illumina33"):
        maxQ = 40
        minQ = 3
        offset = 33
    if (format == "debug"):
        maxQ = 40
        minQ = 3
        offset = 0
        for qval in range(0,len(quality)):
            Q = quality[qval]
            Qs = str(Q) + " "
            qstring = qstring + Qs
    else:
        for qval in range(0,len(quality)):
            Q = int(minQ + (quality[qval] * maxQ))
            if (Q > maxQ):
                Q = maxQ
            Qs = str(chr(Q + offset))
            qstring = qstring + Qs
    return qstring
    ## END OF FUNCTION
    
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="", flush=True, **kwargs)

    ## END OF FUNCTION
    
if __name__ == "__main__":
    main()
