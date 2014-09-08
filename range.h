#ifndef RANGE_H
#define RANGE_H

#include <vector>

struct IRange {
    int start;
    int end;

    int length() const;

    bool operator<(const IRange &other) const;
};

struct IRangeEndPoint {
    int position;
    std::size_t ownerId;
    bool isStart;

    bool operator<(const IRangeEndPoint &other) const;
};

typedef std::vector<std::size_t> IRangeIdCluster;

void processRanges(const std::vector<IRange> &in, std::vector<IRange> &out);
void clusterRanges(const std::vector<IRange> &ranges, std::vector<IRangeIdCluster> &clusters);
void refineRanges(const std::vector<IRange> &in, const std::vector<IRangeIdCluster> &clusters, std::vector<IRange> &out);

#endif // RANGE_H
