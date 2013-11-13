#include <algorithm>
#include <queue>
#include <set>
#include <iterator>
#include "DFinder.h"
#include "ChrRegionCluster.h"

DFinder::DFinder(const std::string& filename, int meanInsertSize, int stdInsertSize, int minOverlapLength, double maxMismatchRate, double discordant) :
    meanInsertSize(meanInsertSize), stdInsertSize(stdInsertSize), minOverlapLength(minOverlapLength), maxMismatchRate(maxMismatchRate), discordant(discordant) {
    if (!r1.Open(filename)) {
	std::cerr << "could not open BAM file" << std::endl;
    }

    // if (!r2.Open(filename)) {
    // 	std::cerr << "could not open BAM file" << std::endl;
    // }
    // if (!r2.LocateIndex()) {
    // 	if (!r2.CreateIndex())
    // 	    std::cerr << "could not find or create index" << std::endl;
    // }

    references = r1.GetReferenceData();
    size = r1.GetReferenceCount();

    leftClips.resize(size);
    leftParts.resize(size);

    rightClips.resize(size);
    rightParts.resize(size);

    intervals.resize(size);

    loadFrom();
    // assert(leftClips[19].size() > 0);
    // assert(rightClips[19].size() > 0);
    for (int i = 0; i < size; ++i) {
	sort(leftClips[i].begin(), leftClips[i].end(), SoftClip::compareL);
	sort(rightClips[i].begin(), rightClips[i].end(), SoftClip::compareR);
    }
}

DFinder::~DFinder() {
    for (int i = 0; i < size; ++i) {
	for_each(leftClips[i].begin(), leftClips[i].end(), DeletePtr<SoftClip>());
	for_each(leftParts[i].begin(), leftParts[i].end(), DeletePtr<SoftClip>());
	for_each(rightClips[i].begin(), rightClips[i].end(), DeletePtr<SoftClip>());
	for_each(rightParts[i].begin(), rightParts[i].end(), DeletePtr<SoftClip>());
	for_each(intervals[i].begin(), intervals[i].end(), DeletePtr<ChrRegion>());
    }
    r1.Close();
    // r2.Close();
}

bool DFinder::isLargeInsertSize(int insertSize) {
    return insertSize >= meanInsertSize + round(discordant * stdInsertSize);
}

void DFinder::loadFrom() {

    int cnt = 0;
    int n_read2 = 0;
    // int n_pairs = 0, found_pairs = 0;

    BamTools::BamAlignment ba1;
    BamTools::BamAlignment ba2;
    int xt1, mq;

    std::map<std::string, int> endPositions;

    while (r1.GetNextAlignment(ba1)) {

	if (!ba1.IsDuplicate() &&
	    ba1.IsPaired() &&
	    ba1.IsMapped() && ba1.IsMateMapped() &&
	    ba1.RefID == ba1.MateRefID &&
	    ba1.IsReverseStrand() != ba1.IsMateReverseStrand() &&
	    ba1.GetTag("XT", xt1)) {
	    std::vector<int> clipSizes, readPositions, genomePositions;
	    if (ba1.GetSoftClips(clipSizes, readPositions, genomePositions)) {
		if (ba1.Position == genomePositions.front() &&
		    ((xt1 == 'M' && ba1.IsProperPair() && !ba1.IsReverseStrand()) ||
		     (xt1 == 'M' && ba1.IsReverseStrand() && ba1.Position > ba1.MatePosition && -ba1.InsertSize > meanInsertSize + 2 * stdInsertSize) ||
		     (xt1 == 'U' && ba1.IsReverseStrand() && !ba1.IsProperPair() && ba1.Position > ba1.MatePosition))) {
		    leftClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
								genomePositions.front(),
								readPositions.front(),
								ba1.QueryBases,
								ba1.Qualities));
		}
		if (ba1.Position != genomePositions.back() &&
		    ((xt1 == 'M' && ba1.IsProperPair() && ba1.IsReverseStrand()) ||
		     (xt1 == 'M' && !ba1.IsReverseStrand() && ba1.Position < ba1.MatePosition && ba1.InsertSize > meanInsertSize + 2 * stdInsertSize) ||
		     (xt1 == 'U' && !ba1.IsReverseStrand() && !ba1.IsProperPair() && ba1.Position < ba1.MatePosition))) {
		    rightClips[ba1.RefID].push_back(new SoftClip(ba1.RefID,
								 genomePositions.back(),
								 ba1.Length - clipSizes.back(),
								 ba1.QueryBases,
								 ba1.Qualities));
		}
	    } else {
		if (ba1.MapQuality >= MapQualityThreshold &&
		    ba1.IsProperPair() && !ba1.IsReverseStrand()) {
		    leftParts[ba1.RefID].push_back(new SoftClip(ba1.RefID,
								ba1.Position,
								0,
								ba1.QueryBases,
								ba1.Qualities));
		} else if (ba1.MapQuality >= MapQualityThreshold &&
			   ba1.IsProperPair() && ba1.IsReverseStrand()) {
		    rightParts[ba1.RefID].push_back(new SoftClip(ba1.RefID,
								 ba1.GetEndPosition(),
								 ba1.Length,
								 ba1.QueryBases,
								 ba1.Qualities));
		}
	    }

	    if (endPositions.count(ba1.Name) && ba1.MapQuality >= MapQualityThreshold) {
		intervals[ba1.RefID].push_back(new ChrRegion(cnt,
							     ba1.Name,
							     ba1.RefID,
							     endPositions[ba1.Name],
							     ba1.Position,
							     -ba1.InsertSize,
							     ba1.Length));
		cnt++;
	    }

	    if (!ba1.IsReverseStrand() && ba1.Position < ba1.MatePosition &&
		ba1.MapQuality >= MapQualityThreshold &&
		isLargeInsertSize(ba1.InsertSize)) {
		endPositions[ba1.Name] = ba1.GetEndPosition();
	    }

	}
    }

    // std::cout << "Number of discordant pairs with one mate marked as 'M': " << n_read2 << std::endl;
}

void DFinder::callToFile(const std::string& filename) {
    std::vector<Deletion> calls;
    call(filename, calls);

    std::ofstream out(filename.c_str());
    out << "Chromosome\tType\tStart\tEnd\tLength" << std::endl;
    for(auto itr = calls.begin(); itr != calls.end(); ++itr) {
	out << references[(*itr).getReferenceId()].RefName << "\tDEL\t" << (*itr).getStart2() << "\t" << (*itr).getEnd2() << "\t" << (*itr).length()
	    << std::endl;
    }
}

void DFinder::callToVcf(const std::string& filename) {
    std::vector<Deletion> calls;
    call(filename, calls);

    std::ofstream out(filename.c_str());
    out << "##fileformat=VCFv4.1" << std::endl;
    out << "##reference=human_g1k_v37" << std::endl;
    out << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << std::endl;
    out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << std::endl;
    out << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << std::endl;
    out << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << std::endl;
    out << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << std::endl;
    out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << std::endl;
    out << "##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"List of samples\">" << std::endl;
    out << "##ALT=<ID=DEL,Description=\"Deletion\">" << std::endl;
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

    for(auto itr = calls.begin(); itr != calls.end(); ++itr) {
	out << references[(*itr).getReferenceId()].RefName << "\t" << (*itr).getStart2() << "\t.\tN\t<DEL>\t.\tPASS\t"
	    << "IMPRECISE;SVTYPE=DEL;END=" << (*itr).getEnd2() << ";SVLEN=-" << (*itr).length()
	    << std::endl;
    }
}

// int DFinder::numOfClipsIn(const TargetRegion& region, const std::vector<SoftClip*>& clips) {
//     std::map<int,std::vector<SoftClip*> > map;
//     for (auto itr = clips.begin(); itr != clips.end(); ++itr) {
//     // for (auto itr = lower_bound(clips.begin(), clips.end(), region.start, SoftClip::compare1);
//     // 	 itr != upper_bound(clips.begin(), clips.end(), region.end, SoftClip::compare2);
//     // 	 ++itr) {
// 	if ((region.start > (*itr)->position() && region.start <= (*itr)->endPosition()) ||
// 	    ((*itr)->position() >= region.start && (*itr)->position() <= region.end) ||
// 	    (region.end >= (*itr)->startPosition() && region.end < (*itr)->position()))
// 	    map[(*itr)->position()].push_back(*itr);
//     }
//     for (auto itr = map.begin(); itr != map.end(); ++itr) {
// 	std::transform(itr->second.begin(), itr->second.end(),
// 		       std::ostream_iterator<const SoftClip&>(std::cout, "\n"),
// 		       [](const SoftClip *sc) { return *sc; } );
//     }
//     std::cout << "######################" << std::endl;
//     return map.size();
// }

void DFinder::call(const std::string& filename, std::vector<Deletion>& calls) {
    for (int i = 0; i < size; i++) {

	std::vector<ChrRegionCluster> clusters;
	if (intervals[i].empty()) continue;
	clusterChrRegions(intervals[i], clusters);

	for (auto itr = clusters.begin(); itr != clusters.end(); ++itr) {
	    // std::cout << *itr << std::endl;
	    Deletion deletion;
	    if (callDeletionInCluster(*itr, deletion)) calls.push_back(deletion);
	}

    }
}

void DFinder::removeLargeChrRegions(std::vector<ChrRegion*>& regions) {
    if (regions.size() < 2) return;
    int l = (*min_element(regions.begin(), regions.end(), [](const ChrRegion* r1, const ChrRegion* r2) { return r1->length() < r2->length(); }))->length();
    regions.erase(remove_if(regions.begin(), regions.end(), [l, stdInsertSize](const ChrRegion *cr) { return cr->length() - l > 6 * stdInsertSize; }), regions.end());
}

void DFinder::clusterChrRegions(const std::vector<ChrRegion*>& remainder, std::vector<ChrRegionCluster>& clusters) {
    std::vector<EndPoint> endpoints;
    for (auto itr = remainder.begin(); itr != remainder.end(); ++itr) {
	endpoints.push_back((*itr)->getStart());
	endpoints.push_back((*itr)->getEnd());
    }
    sort(endpoints.begin(), endpoints.end());
    std::set<int> ids;
    std::queue<ChrRegion*> q;
    for (auto itr = endpoints.begin(); itr != endpoints.end(); ++itr) {
	if ((*itr).isStart()) {
	    q.push((*itr).getOwner());
	}
	else {
	    if (ids.count((*itr).ownerId())) continue;
	    std::vector<ChrRegion*> vec;
	    while (!q.empty()) {
		vec.push_back(q.front());
		ids.insert(q.front()->getId());
		q.pop();
	    }
	    removeLargeChrRegions(vec);
	    ChrRegionCluster clu;
	    for (auto itr2 = vec.begin(); itr2 != vec.end(); ++itr2) clu.add(*itr2);
	    clusters.push_back(clu);
	}
    }
}

bool DFinder::findReferenceId(const std::string& name, int& id) {
    for (int i = 0; i < size; ++i) {
	if (references[i].RefName == name) {
	    id = i;
	    return true;
	}
    }
    return false;
}

void DFinder::printOverlaps(const std::string& filename, int readlength) {
    std::vector<MyInterval> myIntervals;
    loadMyIntervals(filename, myIntervals);
    int cnt = 0, matched = 0;

    for (auto itr = myIntervals.begin(); itr != myIntervals.end(); ++itr) {

	std::cout << "#######################################" << std::endl;
	std::cout << (*itr).refname << ":" << (*itr).start << "-" << (*itr).end << "\t" << (*itr).length << std::endl;
	cnt++;

	int id;
	if (!findReferenceId((*itr).refname, id)) continue;
	auto first1 = upper_bound(leftClips[id].begin(), leftClips[id].end(), (*itr).end - readlength, SoftClip::compare2) - 1;
	auto last1 = upper_bound(leftClips[id].begin(), leftClips[id].end(), (*itr).end, SoftClip::compare2) - 1;
	auto first2 = lower_bound(rightClips[id].begin(), rightClips[id].end(), (*itr).start, SoftClip::compare1);
	auto last2 =  lower_bound(rightClips[id].begin(), rightClips[id].end(), (*itr).start + readlength, SoftClip::compare1);
	Overlap overlap;
	// if (overlaps(first1, last1, first2, last2, ((*itr).length < 150) ? std::max(0, (*itr).length - stdInsertSize) : std::max(0, (*itr).length - 3 * stdInsertSize), ((*itr).length < 150) ? (*itr).length + stdInsertSize : (*itr).length + 3 * stdInsertSize, overlap)) {
	//   std::cout << overlap << std::endl;
	//   matched++;
	// }
    }
    std::cout << "#Matched: " << matched << std::endl;
    std::cout << "Matching rate: " << float(matched) / cnt << std::endl;
}

void DFinder::checkAgainstGoldStandard(const std::string& filename) {
    std::vector<MyInterval> myIntervals;
    loadMyIntervals(filename, myIntervals);
    std::map<int, std::vector<ChrRegionCluster> > map;
    int cnt = 0;

    for (auto itr = myIntervals.begin(); itr != myIntervals.end(); ++itr) {
	int refid;
	if (!findReferenceId((*itr).refname, refid)) continue;
	if (!map.count(refid)) {
	    std::vector<ChrRegionCluster> clusters;
	    if (intervals[refid].empty()) continue;
	    clusterChrRegions(intervals[refid], clusters);
	    map[refid] = clusters;
	}
	if (checkMyInterval(*itr, refid, map[refid])) cnt++;
    }
    std::cout << std::endl << ">>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "Matched: " << cnt << "\t"
	      << "Missed: " << myIntervals.size() - cnt << std::endl;
}

bool DFinder::loadMyIntervals(const std::string& filename, std::vector<MyInterval>& out) {
    std::ifstream in(filename.c_str());

    int start, end, length;
    std::string refname;

    while(in >> refname >> start >> end >> length) {
	out.push_back({refname, start, end, length});
    }
    return true;
}

bool DFinder::checkMyInterval(const MyInterval& myInterval, int refId, const std::vector<ChrRegionCluster>& clusters) {
    bool res = false;
    std::cout << std::endl;
    std::cout << myInterval.refname << ":" << myInterval.start << "-" << myInterval.end << "\t" << myInterval.length << std::endl;
    for (auto itr = clusters.begin(); itr != clusters.end(); ++itr) {
	std::vector<const ChrRegion*> regions;
	(*itr).getOverlaps(myInterval.start, myInterval.end, regions);
	for (auto itr2 = regions.begin(); itr2 != regions.end(); ++itr2) {
	    std::cout << **itr2
		      << "\t" << (*itr2)->minDeletionLength(meanInsertSize, stdInsertSize)
		      << "\t" << (*itr2)->maxDeletionLength(meanInsertSize, stdInsertSize)
		      << std::endl;
	    res = true;
	}
    }
    return res;
}

bool DFinder::callDeletionInCluster(const ChrRegionCluster& cluster, Deletion& deletion) {
    for (auto itr = cluster.end() - 1; itr != cluster.begin() - 1; --itr) {
	Overlap overlap;
	if (getOverlapInRegion(**itr, overlap)) {
	    deletion = overlap.getDeletion();
	    return true;
	}
    }
    return false;
}

bool DFinder::getOverlapInRegion(const ChrRegion& region, Overlap& overlap) {

    std::vector<SoftClip*> part1s;
    getSoftClipsIn(getIntervalOfLeftClips(region), leftClips[region.getReferenceId()], part1s);
    std::map<std::pair<int,int>, std::vector<Overlap> > overlaps;

    for (auto ritr = part1s.rbegin(); ritr != part1s.rend(); ++ritr) {
	std::vector<SoftClip*> part2s;
	getSoftClipsIn(getIntervalOfRightClips(region),
		       rightClips[region.getReferenceId()],
		       part2s);
	for (auto itr = part2s.begin(); itr != part2s.end(); ++itr) {
	    Overlap ov;
	    if ((**itr).overlaps(**ritr, minOverlapLength, maxMismatchRate, ov) &&
		ov.deletionLength() >= lengthThreshold) {
		overlaps[std::make_pair(ov.start(), ov.end())].push_back(ov);
	    }
	}

	std::vector<SoftClip*> otherParts;
	getSoftClipsIn(getIntervalOfRightClips(region),
		       rightParts[region.getReferenceId()],
		       otherParts);
	for (auto itr = otherParts.begin(); itr != otherParts.end(); ++itr) {
	    Overlap ov;
	    if ((**itr).overlapWith(**ritr, minOverlapLength, maxMismatchRate, ov) &&
		ov.deletionLength() >= lengthThreshold) {
		overlaps[std::make_pair(ov.start(), ov.end())].push_back(ov);
		// std::cout << ov << std::endl;
	    }
	}
    }

    std::vector<SoftClip*> parts;
    getSoftClipsIn(getIntervalOfRightClips(region), rightClips[region.getReferenceId()], parts);

    for (auto ritr = parts.rbegin(); ritr != parts.rend(); ++ritr) {
    	std::vector<SoftClip*> otherParts;
    	getSoftClipsIn(getIntervalOfLeftClips(region),
    		       leftParts[region.getReferenceId()],
    		       otherParts);
    	for (auto itr = otherParts.begin(); itr != otherParts.end(); ++itr) {
    	    Overlap ov;
    	    if ((**ritr).overlapWith(**itr, minOverlapLength, maxMismatchRate, ov) &&
		ov.deletionLength() >= lengthThreshold) {
    		overlaps[std::make_pair(ov.start(), ov.end())].push_back(ov);
    		// std::cout << ov << std::endl;
    	    }
    	}
    }

    if (overlaps.empty()) { return false; }

    if (overlaps.size() == 1) {
	return Overlap::getHighScoreOverlap(overlaps.begin()->second, overlap);
    }

    auto itr = overlaps.begin();
    auto res = itr++;
    for (; itr != overlaps.end(); ++itr) {
	if (!contains(res->first, itr->first)) {
	    return Overlap::getHighScoreOverlap(res->second, overlap);
	}
	if ((res->first.second - res->first.first) < (itr->first.second - itr->first.first))
	    res = itr;
	/* std::cout << "(" << itr->first.first << ", " << itr->first.second << ")" << std::endl; */
    }
    /* copy(ovs.begin(), ovs.end(), std::ostream_iterator<Overlap>(std::cout, "\n")); */
    return Overlap::getHighScoreOverlap(res->second, overlap);
}

void DFinder::getSoftClipsIn(const Interval& interval, const std::vector<SoftClip*>& input, std::vector<SoftClip*>& output) {
    copy(lower_bound(input.begin(), input.end(), interval.start, SoftClip::compare1),
		 upper_bound(input.begin(), input.end(), interval.end, SoftClip::compare2),
		 back_inserter(output));
}

Interval DFinder::getIntervalOfLeftClips(const ChrRegion& regionOfInterest) {
    return { regionOfInterest.getEndPos() - meanInsertSize + 2 * regionOfInterest.getReadLength() - 3 * stdInsertSize,
	    regionOfInterest.getEndPos() };
}

Interval DFinder::getIntervalOfRightClips(const ChrRegion& regionOfInterest) {
    return { regionOfInterest.getStartPos(),
	    regionOfInterest.getStartPos() + meanInsertSize - 2 * regionOfInterest.getReadLength() + 3 * stdInsertSize };
}
