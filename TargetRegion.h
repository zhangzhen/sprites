#ifndef _TARGETREGION_H_
#define _TARGETREGION_H_

#include <iostream>

struct TargetRegion {
  int start;
  int end;
  int minDeletionLength;
  int maxDeletionLength;
    int innermostStart;
    int innermostEnd;
    int numOfSupports;
    friend std::ostream& operator <<(std::ostream& os, const TargetRegion& tr) {
	return os << tr.start
		  << "\t" << tr.end
		  << "\t" << tr.minDeletionLength
		  << "\t" << tr.maxDeletionLength
		  << "\t" << tr.innermostStart
		  << "\t" << tr.innermostEnd
		  << "\t" << tr.numOfSupports;
    }
};

/* inline std::ostream& operator <<(std::ostream& os, const TargetRegion& tr) { */
/*     return os << tr.start */
/* 	      << "\t" << tr.end */
/* 	      << "\t" << tr.minDeletionLength */
/* 	      << "\t" << tr.maxDeletionLength */
/* 	      << "\t" << tr.innermostStart */
/* 	      << "\t" << tr.innermostEnd */
/* 	      << "\t" << tr.numOfSupports; */
/* } */

#endif /* _TARGETREGION_H_ */
