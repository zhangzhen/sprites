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
    size_t clippedSize;
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
             size_t clippedSize,
             const std::string& sequence);

    int getReferenceId() const;
    int getPosition() const;
    int getClipPosition() const;
    int getClipPositionInRead() const;
    int getLeftmostPosition() const;

    bool isReverse() const;

    bool isForLeftBp() const;
    bool isForRightBp() const;

    const std::string& getSequence() const;
    size_t size() const;

    friend std::ostream& operator <<(std::ostream& stream, const SoftClip& o);

private:
    bool isTypeIForLeftBp() const;
    bool isTypeIIForLeftBp() const;
    bool isTypeIForRightBp() const;
    bool isTypeIIForRightBp() const;

};

#endif /* _SOFTCLIP_H_ */
