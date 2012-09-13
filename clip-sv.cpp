#include "clip-sv.h"
#include "error.h"

void countAlignments(BamTools::BamReader& reader) {
  int count = 0;
  int count2 = 0;
  int count3 = 0;

  BamTools::BamAlignment al;
  while (reader.GetNextAlignmentCore(al)) {
    std::vector<int> clipSizes, readPositions, genomePositions;
    count++;
    if (!al.IsMapped()) {
      count3++;
      continue;
    }
    if (al.GetSoftClips(clipSizes, readPositions, genomePositions))
      count2++;
    if (count <= 10) {
      std::cout << al.Position << std::endl;
    }
  }

  std::cout << "#alignments: " << count << std::endl;
  std::cout << "#alignments-soft clipping: " << count2 << " " << std::endl;
  std::cout << "#alignments-unmapped: " << count3 << " " << std::endl;

}

void getClips(BamTools::BamReader& reader, std::vector<Clip*>& leftClips, std::vector<Clip*>& rightClips) {
  BamTools::BamAlignment al;
  int report[5] = {0};
  while (reader.GetNextAlignment(al)) {
    std::vector<int> clipSizes, readPositions, genomePositions;
    if (!al.IsMapped()) {
      continue;
    }
    if (!al.GetSoftClips(clipSizes, readPositions, genomePositions)) {
      continue;
    }
    if (clipSizes.size() > 0 && clipSizes.size() <=5)
      report[clipSizes.size()-1]++;
    if (clipSizes.size() > 1) {
      continue;
    }
    if (al.Position == genomePositions[0]) { // left clip - or readPositions[i] == clipSizes[i]
      leftClips.push_back(new LeftClip(al.RefID, genomePositions[0], readPositions[0], clipSizes[0], al.QueryBases));
    } else if (al.Length == readPositions[0]+clipSizes[0]) { // right clip
      rightClips.push_back(new RightClip(al.RefID, genomePositions[0], readPositions[0], clipSizes[0], al.QueryBases));
    }
  }
  outputData(std::cout, report, 5);
}

void countClipsInLength(const std::vector<Clip*>& clips, int stats[], int binWidth) {
  for (int i = 0; i < MAX_BINS; i++) {
    stats[i] = 0;
  }
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int bin = (int)((*itr)->getSize() - 1) / binWidth;
    if (bin > 9) {
      std::cout << (*itr)->toString() << std::endl;
    }
    stats[bin]++;
  }
}

void outputData(std::ostream& output, int lenValues[], int nBins) {
  for (int i = 0; i < nBins; i++) {
    output << lenValues[i] << "\t";
  }
  output << std::endl;
}

Clip* findFirstClipInRange(const std::vector<Clip*>& clips, int min, int max) {
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int s = (*itr)->getSize();
    if (s >= min && s <= max) {
      return *itr;
    }
  }
  
  error("No such clip available.");
}

int countMismatches(const std::string& s1, const std::string& s2) {
  if (s1.size() != s2.size()) {
    error("Two strings must have the same size");
  }
  int count = 0;
  for (int i = 0; i < s1.size(); i++) {
    if (s1[i] != s2[i]) {
      count++;
    }
  }
  return count;
}

int countOverlappedReads(const std::vector<Clip*>& clips, Clip* cl) {
  int count = 0;
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    if ((*itr)->getType() == cl->getType()) {
      continue;
    }
    int len = cl->getSize() + (*itr)->getSize();
    if (len > READ_LENGTH) {
      continue;
    }
    std::string ls, rs;
    if (cl->getType() == Left) {
      ls = cl->getReadSeq().substr(0, len);
      rs = (*itr)->getReadSeq().substr((*itr)->getReadPosition()-cl->getSize());
      // std::cout << ls << std::endl;
      // std::cout << rs << std::endl;
    } else {
      ls = (*itr)->getReadSeq().substr(0, len);
      rs = cl->getReadSeq().substr(cl->getReadPosition()-(*itr)->getSize());
    }
    if (countMismatches(ls, rs) == 0) {
      count++;
    }
  }
  return count;
}

void matchClips(const std::vector<Clip*>& cls1, const std::vector<Clip*>& cls2, std::map<Clip*, std::vector<Clip*> >& matches) {
  for (std::vector<Clip*>::const_iterator itr = cls1.begin(); itr != cls1.end(); ++itr) {
    Clip* cl = *itr;
    for (std::vector<Clip*>::const_iterator itr2 = cls2.begin(); itr2 != cls2.end(); ++itr2) {
      if ((*itr2)->getType() == cl->getType()) {
        continue;
      }
      int len = cl->getSize() + (*itr2)->getSize();
      if (len > READ_LENGTH) {
        continue;
      }
      std::string ls, rs;
      if (cl->getType() == Left) {
        ls = cl->getReadSeq().substr(0, len);
        rs = (*itr2)->getReadSeq().substr((*itr2)->getReadPosition()-cl->getSize());
      } else {
        ls = (*itr2)->getReadSeq().substr(0, len);
        rs = cl->getReadSeq().substr(cl->getReadPosition()-(*itr2)->getSize());
      }
      if (ls == rs) {
        if (!matches.count(cl)) {
          matches[cl] = std::vector<Clip *>();
        }
        matches[cl].push_back(*itr2);
        if (!matches.count(*itr2)) {
          matches[*itr2] = std::vector<Clip *>();
        }
        matches[*itr2].push_back(cl);
      }
    }
  }  
}

void countClipsInLengthOneToFive(const std::vector<Clip*>& clips, int stats[]) {
  for (int i=0; i<5; i++) {
    stats[i] = 0;
  }
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    if ((*itr)->getSize() <= 5) {
      stats[(*itr)->getSize()-1]++;
    }
  }
}

void clusterClips(const std::vector<Clip*>& clips, std::map<int, std::set<Clip*> >& clusters) {
  for (std::vector<Clip*>::const_iterator itr = clips.begin(); itr != clips.end(); ++itr) {
    int pos = (*itr)->getPosition();
    if (!clusters.count(pos))
      clusters[pos] = std::set<Clip*>();
    clusters[pos].insert(*itr);
  }
}

int commonOverlap(const std::string& s1, const std::string& s2, int startLength) {
}
