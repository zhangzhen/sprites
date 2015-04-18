#*sprites*

**Written by** Zhen Zhang (zhangz@csu.edu.cn)  
[Jianxin Wang Lab, Central South University](http://netlab.csu.edu.cn/)

**Please cite:**

---

**Current version:** 0.1.0

Support for Linux

##Summary
*sprites* is a sv caller that specializes in detecting deletion from low-coverage sequencing data. It works by identifying split reads from alignments based on soft-clipping information. By re-aligning a split read to one of its target sequences derived from paired-end reads that span it, a deletion is predicted and breakpoint ends are pinpointed with base-pair resolution. *sprites* uses alignments produced by BWA. Of course, it can also use those produced by other read aligners that support 5'- or 3'-end soft-clipping, like Bowtie2. It can also be extended to detect other types of sv.

##Installation

##Usage
