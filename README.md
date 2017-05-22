# depth_track_window
A script for creating normalized chromosome-wise coverage tracks (in wiggle format) from an input BED file.

The script was used for data processing related to manuscript titled "Androgen receptor deregulation drives bromodomain-mediated chromatin alterations in prostate cancer" by Urbanucci et al.

The script was intended to be run as follows:

<code>python create_depth_track_window.py ifolder ifile ofolder win_size [-zig] [-hg18]</code>

The individual options and their usage:

* "ifolder": path to the folder containing the input BED file
* "ifile": name of the input file (a coordinate-sorted BED file)
* "ofolder": path to the folder for output files
* "win_size":  size of the sliding window used for smooth coverage calculation (even values will in effect be replaced by the closest higher odd value)
* "-zig":  an optional parameter, if "-zig" is present as an argument, bases with zero coverage will be ignored during the computation of mean genomic sequencing depth
* "-hg18": an optional parameter, if "-hg18" is present as an argument, the script will work with chromosome sizes derived from human genomic build hg18 (the defaule is to use chromosome sizes defined by build hg19)
