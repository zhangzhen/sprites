#include "DeletionReader.h"
#include "SoftClipCounter.h"
#include <getopt.h>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

//
// Getopt
//
#define PROGRAM_NAME "countReads"
#define PROGRAM_VERSION "1.0"
#define PROGRAM_BUGREPORT "zhangz@csu.edu.cn"
#define DEFAULT_MIN_CLIP 5

static const char *COUNTREADS_VERSION_MESSAGE =
PROGRAM_NAME " Version " PROGRAM_VERSION "\n"
"Written by Zhen Zhang.\n"
"\n"
"Copyright 2013 netlab.csu.edu.cn\n";

static const char *COUNTREADS_USAGE_MESSAGE =
"Usage: " PROGRAM_NAME " [OPTION] ... BAMFILE\n"
"Count reads in BAMFILE for each deletion given by DELFILE\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -d, --del-file=FILE              all known deletions are stored in FILE\n"
"      -r, --radius=INT                 the radius for search area (default: 50)\n"
"      -c, --min-clip=SIZE              a soft-clip is defined as valid, when the clipped part is not less than SIZE (default: 5)\n"
"          --enhanced-mode              enable the enhanced mode, in which reads of type 2 are considered besides type 1\n"
"\nReport bugs to " PROGRAM_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string bamFile;
    static std::string delFile;
    static unsigned int radius = 50;
    static int mode = 0;

    static unsigned int minClip = DEFAULT_MIN_CLIP;
}

static const char* shortopts = "d:c:r:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_ENHANCED_MODE };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "min-clip",       required_argument, NULL, 'c' },
    { "radius",         required_argument, NULL, 'r' },
    { "del-file",       required_argument, NULL, 'd' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { "enhanced-mode",  no_argument,       NULL, OPT_ENHANCED_MODE },
    { NULL, 0, NULL, 0 }
};

void parseOptions(int argc, char** argv);

void count(const vector<Deletion>& dels);

void showDiff(const vector<int>& positions);

void fetchClippingPositions(vector<int>& positions);

//
// Main
//
int main(int argc, char *argv[]) {
    parseOptions(argc, argv);
    vector<int> positions;
    fetchClippingPositions(positions);
    showDiff(positions);

    // read deletion records from file
//    std::vector<Del> dels = DelReader::readDelsFromFile(opt::delFile);
//    count(dels);

//    std::vector<Del> input;

//    copy_if(dels.begin(), dels.end(), back_inserter(input), [](const Del& d) { return d.isHomogeneous(); });
//    count(input);

//    input.clear();
//    copy_if(dels.begin(), dels.end(), back_inserter(input), [](const Del& d) { return !d.isHomogeneous(); });
//    count(input);

//    input.clear();
//    copy_if(dels.begin(), dels.end(), back_inserter(input), [](const Del& d) { return d.hasInsertedSeq(); });
//    count(input);

//    input.clear();
//    copy_if(dels.begin(), dels.end(), back_inserter(input), [](const Del& d) { return !d.hasInsertedSeq(); });
//    count(input);

//    input.clear();
//    copy_if(dels.begin(), dels.end(), back_inserter(input), [](const Del& d) { return d.hasHomseq(); });
//    count(input);

//    input.clear();
//    copy_if(dels.begin(), dels.end(), back_inserter(input), [](const Del& d) { return !d.hasHomseq(); });
//    count(input);

    return 0;
}

//
// Handle command line arguments
//
void parseOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'c': arg >> opt::minClip; break;
            case 'r': arg >> opt::radius; break;
            case 'd': arg >> opt::delFile; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_ENHANCED_MODE: opt::mode = true; break;
            case OPT_HELP:
                std::cout << COUNTREADS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << COUNTREADS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1)
    {
        std::cerr << PROGRAM_NAME ": missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        std::cerr << PROGRAM_NAME ": too many arguments\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << COUNTREADS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filename
    opt::bamFile = argv[optind++];
}

void count(const vector<Deletion> &dels) {
    SoftClipReader *pReader = new SoftClipReader(opt::bamFile, opt::minClip, opt::mode);
    SoftClipCounter counter(pReader, opt::radius);

    int n1 = 0;
    int n2 = 0;
    int n3 = 0;

    vector<int> labels(dels.size(), -1);

    int i = 0;
    // for each deletion count soft-clips that are used to call it
    for (auto itr = dels.begin(); itr != dels.end(); ++itr) {
        int referenceId = pReader->getReferenceId((*itr).getReferenceName());
        if (referenceId == -1) {
            i++;
            continue;
        }
        cout << ">>>>>>>>>>>>>>>>>>>>>> " << (*itr).getStart2() << "\t[left]" << endl;
        int cntLeftBp = counter.countLeftBp(referenceId, (*itr).getStart2());
        cout << ">>>>>>>>>>>>>>>>>>>>>> " << (*itr).getEnd2() << "\t[right]" << endl;
        int cntRightBp = counter.countRightBp(referenceId, (*itr).getEnd2());
        cout << cntLeftBp << "\t" << cntRightBp << endl;
        int s = cntLeftBp + cntRightBp;
        if (s == 0) {
            labels[i] = 0;
            n1++;
        } else if (cntLeftBp > 0 && cntRightBp > 0) {
            labels[i] = 2;
            n3++;
        } else {
            labels[i] = 1;
            n2++;
        }
        i++;
    }

    // output the counting result to a file
    cout << n1 << "\t" << n2 << "\t" << n3 << endl;

    ofstream out("deletion-labels.txt");
    for(size_t i = 0; i < labels.size(); ++i) {
        out << i << "\t" << labels[i] << endl;
    }

    delete pReader;
}

void showDiff(const vector<int> &positions)
{
    for (size_t i=1; i < positions.size(); ++i) {
        int diff = positions[i] - positions[i-1];
        if (diff < 100)
            cout << positions[i-1] << "\t" << positions[i] << "\t" << diff << endl;
    }
}

void fetchClippingPositions(vector<int> &positions) {
    SoftClipReader reader(opt::bamFile, opt::minClip, opt::mode);
    SoftClip clip;
    while (reader.getSoftClip(clip)) {
        if (clip.isForRightBp()) positions.push_back(clip.getClipPosition());
    }
}
