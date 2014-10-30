#include "range.h"
#include <algorithm>
#include <set>
#include <queue>
#include <cassert>

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


int IRange::length() const
{
    return end - start + 1;
}

bool IRange::operator<(const IRange &other) const
{
    if (start != other.start) return start < other.start;
    return end < other.end;
}

bool IRangeEndPoint::operator<(const IRangeEndPoint &other) const
{
    return position < other.position;
}
