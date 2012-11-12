#!/usr/bin/python
import sys
import ParseFastQ
import argparse
import itertools

def filterReads(inFilename, outFilename, ns):
    parser = ParseFastQ.ParseFastQ(inFilename)  # optional arg: headerSymbols allows changing the header symbols
    output = open(outFilename, 'w')
    total = remained = 0
    for record in parser:
        total += 1
        if 'N'*ns in record[1]:
            continue
        output.write("%s\n%s\n%s\n%s\n" % (record[0], record[1], record[2], record[3]))
        remained += 1
    print "total #reads: %d" % total
    print "remaining #reads: %d" % remained

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
    parser = argparse.ArgumentParser(description='Reads filter')
    group = parser.add_mutually_exclusive_group(required=True)
    subgroup = group.add_argument_group()
    subgroup.add_argument('-1', '--read1')
    subgroup.add_argument('-2', '--read2')
    group.add_argument('-u', '--read')
    parser.add_argument('-o', '--output', help='prefix of output files', required=True)
    parser.add_argument('-n','--number', type=int, default=1, help='number of N')
    args = parser.parse_args()
    if args.read:
        filterReads(args.read, "%s.fq" % args.output, args.number)
    else:
        filterPEReads(args.read1, args.read2, "%s1.fq" % args.output, "%s2.fq" % args.output, args.number)
