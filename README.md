submixbam
=========

A tool for subsampling and mixing reads from two BAM files over regions in a BED file.

### Installation
```
git clone https://github.com/DonFreed/submixbam.git
cd submixbam
make HTSDIR=/path/to/htslib/
```

### Usage

subsetbam [options] [-h inh.sam] <in1.bam> <in2.bam> <in.bed>

Options:
  -h FILE      copy the header in FILE to the output [in1.bam]
  -o FILE      write output to FILE [stdout]
  -m FLOAT     minimum faction of reads from <in1.bam> [0.0]
  -a FLOAT     maximum faction of reads from <in1.bam [1.0]
  -s VALUE     use VALUE to initalize the random seed [time]
  -d FLOAT     downsample from <in1.bam> by this fraction [1.0]
  -e FLOAT     downsample from <in2.bam> by this fraction [1.0]
  -b FILE      write the fractions of reads from <in1.bam> over each region to FILE
  -c INT       compression level for the output file [0]

### Examples

See platinum_200x_example.sh for an example of using submixbam on 200x data from the Illumina Platinum Genomes. 

