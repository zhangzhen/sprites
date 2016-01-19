#include "range.h"
#include <algorithm>
#include <set>
#include <queue>
#include <cassert>
#include <numeric>

using namespace std;

void clusterRanges(const vector<IRange> &ranges, std::vector<IdCluster> &clusters)
{

    vector<IRangeEndPoint> endPoints;
    for (size_t i = 0; i < ranges.size(); ++i) {
        endPoints.push_back({ranges[i].start, i, true});
        endPoints.push_back({ranges[i].end, i, false});
    }

    sort(endPoints.begin(), endPoints.end());
    set<size_t> usedIds;
    queue<size_t> buffer;

    for (auto it = endPoints.begin(); it != endPoints.end(); ++it) {
        if ((*it).isStart) buffer.push((*it).ownerId);
        else {
            if (usedIds.count((*it).ownerId)) continue;
            IdCluster clu;
            while (!buffer.empty()) {
                clu.push_back(buffer.front());
                usedIds.insert(buffer.front());
                buffer.pop();
            }
            if (!clu.empty()) clusters.push_back(clu);
        }
    }
    IdCluster clu;
    while (!buffer.empty()) {
        clu.push_back(buffer.front());
        usedIds.insert(buffer.front());
        buffer.pop();
    }
    if (!clu.empty()) clusters.push_back(clu);

}


void append(size_t startIndex, size_t endIndex, std::vector<IdCluster> &clusters)
{
    IdCluster buffer(endIndex - startIndex);
    std::iota(std::begin(buffer), std::end(buffer), startIndex);
    clusters.push_back(buffer);
}

void clusterRanges2(const vector<IRange> &ranges, std::vector<IdCluster> &clusters)
{
    size_t startIndex = 0;

    for (size_t i = 1; i < ranges.size(); ++i) {
        if (!ranges[i-1].overlaps(ranges[i])) {
            append(startIndex, i, clusters);
            startIndex = i;
        }
    }
    append(startIndex, ranges.size(), clusters);

}


int IRange::length() const
{
    return end - start + 1;
}

bool IRange::operator<(const IRange &other) const
{
    if (start != other.start) return start < other.start;
    return end < other.end;
}

bool IRange::overlaps(const IRange &other) const
{
    return (start >= other.start && start < other.end) ||
            (other.start >= start && other.start < end);
}

bool IRangeEndPoint::operator<(const IRangeEndPoint &other) const
{
    return position < other.position;
}
