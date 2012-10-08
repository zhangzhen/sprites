import sys
import fileinput

def countbases(filename):
    f = open(filename, "r")
    line = f.readline()
    if line[0] != ">":
        f.close()
        raise ValueError, "file not in proper format"
    cnt = 0
    while True:
        line = f.readline()
        if line == "":
            break
        cnt += len(line.rstrip('\n'))
    f.close()
    return cnt

def changeheaderline(filename):
    cnt = countbases(filename)
    for line in fileinput.input(filename, inplace=1):
        if fileinput.isfirstline():
            print ">chr0:1-%d" % cnt
        else:
            sys.stdout.write(line)

def main():
    changeheaderline(sys.argv[1])

if __name__ == '__main__':
    main()
