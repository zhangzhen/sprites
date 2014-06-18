#ifndef _SOFTCLIP_H_
#define _SOFTCLIP_H_

#include <string>

class SoftClip {
private:
    int referenceId;
    int position;
    int leftmostPosition;
    int clipPosition;
    int matePosition;
    bool is_reverse;
    bool is_mate_reverse;
    int clippedSize;
    int offset;
    std::string sequence;

public:
    SoftClip();
    SoftClip(int referenceId,
             int position,
             int leftmostPosition,
             int clipPosition,
             int matePosition,
             bool is_reverse,
             bool is_mate_reverse,
             int clippedSize,
             int offset,
             const std::string& sequence);

    int getReferenceId() const;
    int getPosition() const;
    int getClipPosition() const;
    int getAdjustedPosition() const;
    int getLeftmostPosition() const;

    bool isReverse() const;

    bool isTypeIForLeftBp() const;
    bool isTypeIIForLeftBp() const;
    bool isTypeIForRightBp() const;
    bool isTypeIIForRightBp() const;

    bool isForLeftBp() const;
    bool isForRightBp() const;

    const std::string& getSequence() const;
    int size() const;
    int getClippedSize() const;

    friend std::ostream& operator <<(std::ostream& stream, const SoftClip& o);

};

#endif /* _SOFTCLIP_H_ */
