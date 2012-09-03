import pylab

def main(filename):
    f = open(filename, "r")
    line = f.readline()
    data = line.strip().split("\t")
    data = map(int, data)
    print data
    #pylab.hist(data)

if __name__ == "__main__":
    main("results.txt")
