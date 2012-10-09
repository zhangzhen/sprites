#!/usr/bin/python
import sys
import ParseFastQ

def filterReads(inFilename, outFilename, numN):
    parser = ParseFastQ.ParseFastQ(inFilename)  # optional arg: headerSymbols allows changing the header symbols
    output = open(outFilename, 'w')
    total = 0
    rest = 0
    for record in parser:
        total += 1
        if 'N'*numN in record[1]:
            continue
        output.write("%s\n%s\n%s\n%s\n" % (record[0], record[1], record[2], record[3]))
        rest += 1
    print "total #reads: %d" % total
    print "remaining #reads: %d" % rest

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reads filter')
    parser.add_argument('-i','--input', help='input fastq file', required=True)
    parser.add_argument('-o','--output', help='output fastq file', required=True)
    parser.add_argument('-n','--number', type=int, default=1, help='number of N')
    filterReads(parser.input, parser.output, parser.number)
