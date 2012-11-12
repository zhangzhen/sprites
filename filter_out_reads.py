#!/usr/bin/python
import sys
import ParseFastQ
import argparse

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='A tool that filters out single-end reads containing at least a specified number of Ns')
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('-n','--number', type=int, default=1, help='number of N')
    args = parser.parse_args()
    filterReads(args.input, args.output, args.number)
