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
