#include <cmath>
#include <ctime>
#include <algorithm>
#include <getopt.h>
#include <sstream>

#include "error.h"
#include "Deletion.h"
#include "SoftClipReader.h"
#include "BamStatCalculator.h"
#include "Caller.h"
#include "Helper.h"
#include "Parameters.h"
#include "clip.h"

//
// Getopt
//
#define PROGRAM_NAME "dfinder"
#define PROGRAM_VERSION "1.0"
#define PROGRAM_BUGREPORT "zhangz@csu.edu.cn"
#define DEFAULT_MIN_OVERLAP 21

static const char *DFINDER_VERSION_MESSAGE =
PROGRAM_NAME " Version " PROGRAM_VERSION "\n"
"Written by Zhen Zhang.\n"
"\n"
"Copyright 2013 netlab.csu.edu.cn\n";

static const char *DFINDER_USAGE_MESSAGE =
"Usage: " PROGRAM_NAME " [OPTION] ... BAMFILE\n"
"Find deletions from records in BAMFILE\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -r, --reffile=FILE           read the reference sequence from FILE"
"      -o, --outfile=FILE               write the deletion calls to FILE (default: BAMFILE.calls)\n"
"      -e, --error-rate=F               the maximum error rate allowed between two sequences to consider them overlapped (default: 0.04)\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads (default: 31)\n"
"      -n, --allowed-num=SIZE              a soft-clip is defined as valid, when the clipped part is not less than SIZE (default: 5)\n"
"\nThe following two option must appear together (if ommitted, attempt ot learn the mean and the standard deviation of insert size):\n"
"      -i, --insert-mean=N              the mean of insert size\n"
"          --enhanced-mode              enable the enhanced mode, in which reads of type 2 are considered besides type 1\n"
"      -s, --insert-sd=N                the standard deviation of insert size\n"
"\nReport bugs to " PROGRAM_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string bamFile;
    static std::string refFile;
    static std::string outFile;
    static double errorRate = 0.04;
    static int minOverlap = DEFAULT_MIN_OVERLAP;
    static int allowedNum = 5;
    static int mode = 0;

    static bool bLearnInsert = true;
    static int insertMean;
    static int insertSd;
}

static const char* shortopts = "o:r:e:m:n:i:s:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_ENHANCED_MODE };

static const struct option longopts[] = {
    { "verbose",        no_argument,       NULL, 'v' },
    { "min-overlap",    required_argument, NULL, 'm' },
    { "allowed-num",       required_argument, NULL, 'n' },
    {"reffile", required_argument, NULL, 'r'},
    { "outfile",        required_argument, NULL, 'o' },
    { "error-rate",     required_argument, NULL, 'e' },
    { "insert-mean",    required_argument, NULL, 'i' },
    { "insert-sd",      required_argument, NULL, 's' },
    { "help",           no_argument,       NULL, OPT_HELP },
    { "version",        no_argument,       NULL, OPT_VERSION },
    { "enhanced-mode",  no_argument,       NULL, OPT_ENHANCED_MODE },
    { NULL, 0, NULL, 0 }
};

void parseOptions(int argc, char** argv);

//
// Main
//
int main(int argc, char *argv[]) {
    parseOptions(argc, argv);

    if (opt::bLearnInsert) {
        std::cout << "Estimate the mean and standard deviation of insert size:" << std::endl;
        BamStatCalculator calc(opt::bamFile);
        opt::insertMean = calc.getInsertMean();
        opt::insertSd = calc.getInsertSd();
        std::cout << "Mean: " << opt::insertMean << std::endl;
        std::cout << "Sd: " << opt::insertSd << std::endl;
    }

    Parameters params = { opt::allowedNum,
                          opt::mode,
                          opt::minOverlap,
                          1.0f - opt::errorRate,
                          opt::insertMean,
                          opt::insertSd };

/*
    BamTools::BamReader bamReader;
    bamReader.Open(opt::bamFile);

    ClipReader creader(opt::bamFile, opt::allowedNum, opt::mode);
    std::vector<Deletion> deletions;

    AbstractClip *pClip;
    while (!(pClip = creader.nextClip())) {
        try {
            deletions.push_back(pClip->call(bamReader, faidx,
                                            opt::insertMean + 3 * opt::insertSd, opt::minOverlap, 1.0f - opt::errorRate));
        } catch (ErrorException& ex) {
        }
    }
*/

    Caller caller(opt::bamFile, params);

    std::vector<SoftClip> clips;
    caller.readClipsForRightBp(clips);

    std::vector<Deletion> dels;
    caller.call(clips, dels);

    caller.output(opt::outFile, dels);

    return 0;
}

//
// Handle command line arguments
//
void parseOptions(int argc, char** argv)
{
    bool bInsertMean = false;
    bool bInsertSd = false;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'n': arg >> opt::allowedNum; break;
            case 'm': arg >> opt::minOverlap; break;
            case 'r' : arg >> opt::refFile; break;
            case 'o': arg >> opt::outFile; break;
            case 'e': arg >> opt::errorRate; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case 'i': arg >> opt::insertMean; bInsertMean = true; break;
            case 's': arg >> opt::insertSd; bInsertSd = true; break;
            case OPT_ENHANCED_MODE: opt::mode = true; break;
            case OPT_HELP:
                std::cout << DFINDER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << DFINDER_VERSION_MESSAGE;
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

    if (bInsertMean & bInsertSd) {
        opt::bLearnInsert = false;
    }

    if (bInsertMean ^ bInsertSd) {
        std::cerr << PROGRAM_NAME ": the mean and standard deviation of insert size must be specified together\n";
        die = true;
    }

    if(opt::errorRate > 1.0f)
    {
        std::cerr << PROGRAM_NAME ": invalid error-rate parameter: " << opt::errorRate << "\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << DFINDER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Validate parameters
    if(opt::errorRate <= 0)
        opt::errorRate = 0.0f;

    // Parse the input filename
    opt::bamFile = argv[optind++];

    std::string out_prefix = stripFilename(opt::bamFile);
    if(opt::outFile.empty())
    {
        opt::outFile = out_prefix + ".calls";
    }

}
