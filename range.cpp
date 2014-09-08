#include "range.h"
#include <algorithm>
#include <set>
#include <queue>

using namespace std;

void clusterRanges(const vector<IRange> &ranges, vector<IRangeIdCluster> &clusters)
{
    vector<IRangeEndPoint> endPoints;
    for (auto i =0; i < ranges.size(); ++i) {
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
            IRangeIdCluster clu;
            while (!buffer.empty()) {
                clu.push_back(buffer.front());
                usedIds.insert(buffer.front());
                buffer.pop();
            }
            if (!clu.empty()) clusters.push_back(clu);
        }
    }
    IRangeIdCluster clu;
    while (!buffer.empty()) {
        clu.push_back(buffer.front());
        usedIds.insert(buffer.front());
        buffer.pop();
    }
    if (!clu.empty()) clusters.push_back(clu);

}


void refineRanges(const std::vector<IRange> &in, const std::vector<IRangeIdCluster> &clusters, std::vector<IRange> &out)
{
    for (auto i = clusters.begin(); i != clusters.end(); ++i) {
        vector<int> rangeLengths;
        for (auto j = (*i).begin(); j != (*i).end(); ++j) {
            rangeLengths.push_back(in[*j].length());
        }
    }
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

void processRanges(const std::vector<IRange> &in, std::vector<IRange> &out)
{

}


bool IRangeEndPoint::operator<(const IRangeEndPoint &other) const
{
    return position < other.position;
}
