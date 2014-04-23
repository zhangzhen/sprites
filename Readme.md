Here, we describe a new method for detecting deletions from
sequencing data, DFinder, which relies on overlaps between reads
across a deletion site to identify breakpoints of deletions. DFinder
works as follows: first, it identifies discordant read-pairs from the
input alignment file and then clusters these discordant read-pairs
to obtain an area that indicates roughly where a deletion might
occur. In the area, the search is performed to detect overlaps
between reads. A final call is then determined from the overlaps.
The minimum overlap size can play a role in preventing split-reads
with short clipped part from being discarded too early. Search in
the local area ensures time efficiency. Clustering makes our method
appliable to high-coverage data without losing time and memory
efficiency. We use both simulated and real data to evaluate our
caller compared to other detection tools with base-pair breakpoint
resolution. Experiments on simulated data show that DFinder can
deal with both low-coverage and high-coverage data. The result
of evaluation gives evidence that DFinder performs well on theses
dataset.
