#ifndef RANGE_H
#define RANGE_H

#include <vector>
#include <string>

struct IRange {
    int start;
    int end;

    int length() const;

    bool operator<(const IRange &other) const;
    bool overlaps(const IRange& other) const;
};

struct IRangeEndPoint {
    int position;
    std::size_t ownerId;
    bool isStart;

    bool operator<(const IRangeEndPoint &other) const;
};

typedef std::vector<std::size_t> IdCluster;

void clusterRanges(const std::vector<IRange> &ranges, std::vector<IdCluster> &clusters);
void clusterRanges2(const std::vector<IRange> &ranges, std::vector<IdCluster> &clusters);

#endif // RANGE_H
