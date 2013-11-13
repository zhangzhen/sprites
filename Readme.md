DFinder: a deletion finder using discordant read pairs
  and overlaps between soft-clips

Methods of discovering overlaps between soft-clips
or between a soft-clip and a read have three advantages in comparison
to those of re-mapping unmapped segments of soft-clips:

First, the clip position of a read given by BWA where soft-clipping occurs is
usually inaccurate because of sequencing errors. For example, the clip
position might occur either before or after the real one. Methods
relying on
remapping need to take additional steps to adjust calls afterward in
order to obtain accurate variant identification. Those relying on
overlaps consider the fact while making calls, so they do not need such
post-processing.

Second, although the donor and reference genomes have huge similarity, they
differ in many places. Among those differences, big
ones are called structural variants, whereas small ones
includes indels ($\leq$ 50 bps) and SNPs. Thus, finding overlaps between reads
of the same genome
 leads to better result than re-mapping a small segment of reads on
 the reference.

Third, reads of length 35 or 50 tend to result in very short unmapped
segments, which causes difficulty in remapping them to the reference
because there are too many hits to decide which one is the
best. However, the length of an overlap is usually larger than the sum
of unmapped segments of two corresponding soft-clips, so they can overcome the
limitation and
are good at dealing with those reads.

我们看一个例子：
Ref: 	AGCATGTTAGATA[AGATAGCTGTGCTA]GTAGGCAGTCAGCGCCAT
Donor:	AGCATGTTAGATAGTAGGCAGTCAGCGCCAT
Donor是Ref方括号之间的部分被删除得到的。

有一种特殊的reads，他比对到Ref上，一部分比对上了，另一部分没有。read没有比对上的部分就叫clip,形象的看起来就是被剪掉的。因为发生剪切的位点很有可能是基因组结构变异的断裂点(breakpoint)，这样的reads对检测基因组结构变异是非常有用的，我们叫他split reads。ATAGTAGGCA和TAGATAGTAG是Donor的两个split reads。
一个read从左边被剪掉的部分就是left clip。An example for left clips: (ATA)GTAGGCA - 圆括号包含的是被剪掉的。这个clip的长度是3。相应的，一个read从右边被剪掉的部分就是right clip。An example for right clips: TAGATA(GTAG) - 这个clip的长度是4。
上面的一个left clip和一个right clip组成了支持Donor中那个删除的一对clip。观察得到，代表这两个clips的reads之间是有overlap。因为所有的读数都是等长的，发生overlap只有一种情况。


问题描述
输入：两个集合 - leftClips and rightClips。leftClips中的一个left clip可能与rightClips中的多个right clips配对。
目标：找出clip对，最终确定结构变异。（限制条件需要确定）
comments:
感觉有点想最大匹配问题。

All alignments (including unique, multiple, and unmapped) are kept in <NAME>.bam.

The multiple.bam actually is like a mapping distribution. Once we got multiple mappings for a mate/read, we record the chromosome ids and positions of them in multiple.bam. More precisely, if -om is not enabled, only chromosome ids and positions are shown in multiple.bam. bamtools (https://github.com/pezmaster31/bamtools) can handle that bam.

read position for a right clip: the first base position of the clipped sequence
read position for a left clip: the first base position of the aligned sequence
