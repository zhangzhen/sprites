#include "DelReader.h"
#include "SoftClipCounter.h"
#include <getopt.h>
#include <sstream>

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
"\nReport bugs to " PROGRAM_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string bamFile;
    static std::string delFile;
    static unsigned int radius = 50;

    static unsigned int minClip = DEFAULT_MIN_CLIP;
}

static const char* shortopts = "d:c:r:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "min-clip",      required_argument, NULL, 'c' },
    { "radius",        required_argument, NULL, 'r' },
    { "del-file",      required_argument, NULL, 'd' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void parseOptions(int argc, char** argv);

//
// Main
//
int main(int argc, char *argv[]) {
    parseOptions(argc, argv);

    // read deletion records from file
    std::vector<Del> dels = DelReader::readDelsFromFile(opt::delFile);

    SoftClipReader *pReader = new SoftClipReader(opt::bamFile, opt::minClip);
    SoftClipCounter counter(pReader, opt::radius);

    // for each deletion count soft-clips that are used to call it    
    for (auto itr = dels.begin(); itr != dels.end(); ++itr) {
        int referenceId = pReader->getReferenceId((*itr).referenceName);
        if (referenceId == -1) continue;
        int cntLeftBp = counter.count(referenceId, (*itr).leftBp);
        int cntRightBp = counter.count(referenceId, (*itr).rightBp);
        cout << cntLeftBp << "\t" << cntRightBp << endl;
    }

    // output the counting result to a file
    delete pReader;

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
