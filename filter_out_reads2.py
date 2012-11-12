#!/usr/bin/python
import sys
import ParseFastQ
import argparse
import itertools

def filterPEReads(inFilename1, inFilename2, outFilename1, outFilename2, ns):
    pa1 = ParseFastQ.ParseFastQ(inFilename1)
    pa2 = ParseFastQ.ParseFastQ(inFilename2)
    out1 = open(outFilename1, 'w')
    out2 = open(outFilename2, 'w')
    total = remained = 0
    for r1, r2 in itertools.izip(pa1, pa2):
        total += 1
        if 'N'*ns in r1[1] or 'N'*ns in r2[1]:
            continue
        out1.write("%s\n%s\n%s\n%s\n" % (r1[0], r1[1], r1[2], r1[3]))
        out2.write("%s\n%s\n%s\n%s\n" % (r2[0], r2[1], r2[2], r2[3]))
        remained += 1
    print "total #reads: %d" % total
    print "remaining #reads: %d" % remained
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='A tool that filters out paired-end reads containing at least a specified number of Ns')
    parser.add_argument('input1')
    parser.add_argument('input2')    
    parser.add_argument('output1')
    parser.add_argument('output2')
    parser.add_argument('-n','--number', type=int, default=1, help='number of N')
    args = parser.parse_args()
    filterPEReads(args.input1, args.input2, args.output1, args.output2, args.number)
