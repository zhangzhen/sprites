#*Sprites*

**Written by** Zhen Zhang (zhangz@csu.edu.cn)  
[Jianxin Wang Lab, Central South University](http://netlab.csu.edu.cn/)

**Please cite:**

---

**Current version:** 0.1.0

Support for Linux

##Summary
*Sprites* is a sv caller that specializes in detecting deletion from low-coverage sequencing data. It works by identifying split reads from alignments based on soft-clipping information. By re-aligning a split read to one of its target sequences derived from paired-end reads that span it, a deletion is predicted and breakpoint ends are pinpointed with base-pair resolution. *Sprites* uses alignments produced by BWA. Of course, it can also use those produced by other read aligners that support 5'- or 3'-end soft-clipping, like Bowtie2. It can also be extended to detect other types of sv.

##Installation

#### Requirements
- HTSlib ([http://www.htslib.org/](http://www.htslib.org/))
- BamTools ([https://github.com/pezmaster31/bamtools](https://github.com/pezmaster31/bamtools))
- CMake ([http://www.cmake.org](http://www.cmake.org))

##### Install
```
git clone git@github.com:zhangzhen/sprites.git
cd sprites
export BAMTOOLS_HOME=/path/to/bamtools
export HTSLIB_HOME=/path/to/htslib
mkdir build
cd build
cmake ..
make
cp build/sprites /usr/local/bin/
```
##Usage
```
sprites [options] sample.bam
```
The input bam file is required to be sorted.

**Options**
```
-r FILE 
```