#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include "htslib/hts.h"
#include "htslib/sam.h"

#define IS_EMPTY UINT64_MAX

/* htslib/bam_sort.c license */

/*  bam_sort.c -- sorting and merging.

    Copyright (C) 2008-2014 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

typedef struct {
    int idx;
    uint64_t pos;
    bam1_t *b;
} heap2_t;

#define heap2_lt(a, b) ((a).pos < (b).pos || ((a).pos == (b).pos && ((a).idx < (b).idx || ((a).idx == (b).idx && rand() % 2))))

// bedidx declarations //
void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);
uint64_t bed_get_reg(const void *_h, const char *chr, int beg, int end);

void usage(FILE *fp)
{
    fprintf(fp,
"\n"
"submixbam\n"
"  Randomly subset from and merge sequence reads from two bam files over regions\n"
"\n"
"Usage: subsetbam [options] [-h inh.sam] <in1.bam> <in2.bam> <in.bed>\n"
"\n"
"Options:\n"
"  -h FILE      copy the header in FILE to the output [in1.bam]\n"
"  -o FILE      write output to FILE [stdout]\n"
"  -m FLOAT     minimum faction of reads from <in1.bam> [0.0]\n"
"  -a FLOAT     maximum faction of reads from <in1.bam [1.0]\n"
"  -s VALUE     use VALUE to initalize the random seed [time]\n"
"  -d FLOAT     downsample from <in1.bam> by this fraction [1.0]\n"
"  -e FLOAT     downsample from <in2.bam> by this fraction [1.0]\n"
"  -b FILE      write the fractions of reads from <in1.bam> over"
" each region to FILE\n"
"  -c INT       compression level for the output file [0]\n");
}

/*
 * Read a sam/bam file with downsampling
 */
inline void sam_downsample_next(samFile *fp, hts_itr_t *iter, heap2_t *h, float down, uint32_t max_tid)
{
    do {
        if (sam_itr_next(fp, iter, h->b) >= 0) {
            uint64_t tmp_pos = ((uint64_t)h->b->core.tid<<32) | (uint32_t)((uint32_t)h->b->core.pos+1)<<1 | bam_is_rev(h->b);
            if (tmp_pos == h->pos) {
                h->idx++;
            } else {
                h->idx = 0;
                h->pos = tmp_pos;
            }
        } else {
            h->pos = IS_EMPTY;
            return;
        }
    } while (drand48() > down);
    if (h->b->core.tid >= max_tid) { // Unaligned reads
        h->pos = IS_EMPTY;
        return;
    }
}

int main(int argc, char *argv[])
{
    int c, compression = 0;
    char *fn_in1 = 0, *fn_in2 = 0, *fn_bed, *fn_out = 0, *fn_hdr = 0, mode[16], *fn_bed_out = 0;
    long random_seed = (long)time(NULL);
    float min_frac = 0.0, max_frac = 1.0, down1 = 1.0, down2 = 1.0, diff_frac, fraction = 0.0;
    void *bed = 0;
    samFile *fp_out = 0, *fp_in1 = 0, *fp_in2 = 0;
    bam_hdr_t *hdr = 0, *hout = 0;
    hts_idx_t *idx = 0;
    hts_itr_t *iter1 = 0, *iter2 = 0;
    uint64_t reg = 0, last_reg = 1;
    FILE *fp_bed_out = 0;
    heap2_t h1, h2;

    while ((c = getopt(argc, argv, "b:h:o:m:a:s:d:e:c:")) >= 0) {
        if (c == 'h') fn_hdr = strdup(optarg);
        else if (c == 'o') fn_out = strdup(optarg);
        else if (c == 'm') min_frac = atof(optarg);
        else if (c == 'a') max_frac = atof(optarg);
        else if (c == 's') random_seed = atol(optarg);
        else if (c == 'd') down1 = atof(optarg);
        else if (c == 'e') down2 = atof(optarg);
        else if (c == 'c') compression = atoi(optarg);
        else if (c == 'b') fn_bed_out = strdup(optarg);
    }
    // User errors //
    if (argc - optind< 3 || argc - optind > 3) {
        usage(stderr);
        return 1;
    }
    if (min_frac < 0.0 || max_frac < 0.0) { 
        fprintf(stderr, "Warning: minimum fraction of reads can not be set below 0\n");
        min_frac = min_frac < 0.0 ? 0.0 : min_frac;
        max_frac = max_frac < 0.0 ? 0.0 : max_frac;
    }
    if (min_frac > 1.0 || max_frac > 1.0) {
        fprintf(stderr, "Warning: minimum fraction of reads can not be set above 1.0\n");
        min_frac = min_frac > 1.0 ? 1.0 : min_frac;
        min_frac = max_frac > 1.0 ? 1.0 : max_frac;
    }
    if (min_frac >= max_frac) {
        fprintf(stderr, "Error: the minimum fraction of reads must be below the maximum fraction of reads\n");
        return 1;
    }
    if (down1 > 1.0 || down1 < 0.0) {
        fprintf(stderr, "Error: '-d' must be between 0.0 and 1.0\n");
        return 1;
    }
    if (down2 > 1.0 || down2 < 0.0) {
        fprintf(stderr, "Error: '-e' must be between 0.0 and 1.0\n");
        return 1;
    }
    if (compression < 0)
        compression = 0;
    else if (compression > 9)
        compression = 9;

    srand48(random_seed);
    fn_in1 = strdup(argv[optind]);
    fn_in2 = strdup(argv[optind + 1]);
    fn_bed = strdup(argv[optind + 2]);
    if (!(bed = bed_read(fn_bed))) {
        fprintf(stderr, "Error: error reading bed file\n");
        return 1;
    }
    // The hard way //
    mode[0] = 'w'; mode[1] = 'b'; mode[2] = compression + '0'; mode[3] = 0;
    if (!(fp_out = sam_open(fn_out ? fn_out : "-", mode))) {
        fprintf(stderr, "Error: could not open output file %s\n", fn_out);
        return 1;
    }
    if (!(fp_in1 = sam_open(fn_in1, "rb"))) {
        fprintf(stderr, "Error: could not open input file %s\n", fn_in1);
        return 1;
    }
    if (!(fp_in2 = sam_open(fn_in2, "rb"))) {
        fprintf(stderr, "Error: could not open input file %s\n", fn_in2);
        return 1;
    }
    if (fn_bed_out) {
        if (!(fp_bed_out = fopen(fn_bed_out, "w"))) {
            fprintf(stderr, "Error: could not open output file %s\n", fn_bed_out);
            return 1;
        }
    }

    // Read the header and write to the output file //
    hdr = bam_hdr_init();
    hout = bam_hdr_init();
    if (fn_hdr) {
        htsFile *fp_hdr = hts_open(fn_hdr, "r");
        if (!fp_hdr) {
            fprintf(stderr, "Error: could not read header file\n");
            return 1;
        }
        hout = sam_hdr_read(fp_hdr);
        if (sam_hdr_write(fp_out, hout)) {
            fprintf(stderr, "Error: could not write header\n");
            return 1;
        }
        // Read other headers
        hdr = sam_hdr_read(fp_in1);
        bam_hdr_destroy(hdr);
        hdr = sam_hdr_read(fp_in2);
        bam_hdr_destroy(hdr);
        hts_close(fp_hdr);
    } else {
        hout = sam_hdr_read(fp_in1);
        if (sam_hdr_write(fp_out, hout)) {
            fprintf(stderr, "Error: could not wirte header\n");
            return 1;
        }
        // Read other header
        hdr = sam_hdr_read(fp_in2);
        bam_hdr_destroy(hdr);
    }

    // Initalize iterators //
    if (!(idx = sam_index_load(fp_in1, fn_in1))) {
        fprintf(stderr, "Could not load index of file %s\n", fn_in1);
        return 1;
    }
    if (!(iter1 = sam_itr_queryi(idx, HTS_IDX_START, 0, 0))) {
        fprintf(stderr, "Could not load iterator for file %s\n", fn_in1);
        return 1;
    }
    hts_idx_destroy(idx);
    if (!(idx = sam_index_load(fp_in2, fn_in2))) {
        fprintf(stderr, "Could not load index of file %s\n", fn_in2);
        return 1;
    }
    if (!(iter2 = sam_itr_queryi(idx, HTS_IDX_START, 0, 0))) {
        fprintf(stderr, "Could not load iterator for file %s\n", fn_in2);
        return 1;
    }

    hts_idx_destroy(idx);
    h1.b = bam_init1();
    h2.b = bam_init1();
    if ((c = sam_itr_next(fp_in1, iter1, h1.b)) >= 0) {
        h1.pos = ((uint64_t)h1.b->core.tid<<32) | (uint32_t)((uint32_t)h1.b->core.pos+1)<<1 | bam_is_rev(h1.b);
        h1.idx = 0;
    } else {
        h1.pos = IS_EMPTY;
    }
    if ((c = sam_itr_next(fp_in2, iter2, h2.b)) >= 0) {
        h2.pos = ((uint64_t)h2.b->core.tid<<32) | (uint32_t)((uint32_t)h2.b->core.pos+1)<<1 | bam_is_rev(h2.b);
    } else {
        h2.pos = IS_EMPTY;
    }
    // Merge data //
    diff_frac = max_frac - min_frac;
    //fraction = (float)(drand48() * diff_frac + min_frac); // fraction of reads from bam1
    while (h1.pos != IS_EMPTY || h2.pos != IS_EMPTY) {
        if (drand48() > fraction) { // draw from bam2
            reg = bed_get_reg(bed, hout->target_name[h2.b->core.tid], h2.b->core.pos, bam_endpos(h2.b));
            if (reg != last_reg) {
                // Move to a new region //
                while (reg == 0) {
                    sam_downsample_next(fp_in2, iter2, &h2, down2, hout->n_targets);
                    if (h2.pos == IS_EMPTY)
                        break;
                    reg = bed_get_reg(bed, hout->target_name[h2.b->core.tid], h2.b->core.pos, bam_endpos(h2.b));
                }
                if (h2.pos == IS_EMPTY) break;
                // Reads of different lengths or clipping might cause reg == last_reg //
                if (last_reg != reg) {
                    last_reg = reg;
                    // update bam1 //
                    reg = bed_get_reg(bed, hout->target_name[h1.b->core.tid], h1.b->core.pos, bam_endpos(h1.b));
                    while (heap2_lt(h1, h2) && reg != last_reg) {
                        sam_downsample_next(fp_in1, iter1, &h1, down1, hout->n_targets);
                        if (h1.pos == IS_EMPTY)
                            break;
                        reg = bed_get_reg(bed, hout->target_name[h1.b->core.tid], h1.b->core.pos, bam_endpos(h1.b));
                    }
                    fraction = (float)(drand48() * diff_frac + min_frac); // fraction of read from bam1
                    if (fp_bed_out)
                        fprintf(fp_bed_out, "%s\t%" PRIu32 "\t%" PRIu32 "\t%.3f\n", hout->target_name[h2.b->core.tid], (uint32_t)(last_reg >> 32), ((uint32_t)last_reg), fraction);
                    continue; // draw again for new region
                }
            }
            sam_write1(fp_out, hout, h2.b);
            // move bam1, if necessary
            while (heap2_lt(h1, h2)) {
                sam_downsample_next(fp_in1, iter1, &h1, down1, hout->n_targets);
            }
            // move bam2
            sam_downsample_next(fp_in2, iter2, &h2, down2, hout->n_targets);
        } else { // draw from bam1
            reg = bed_get_reg(bed, hout->target_name[h1.b->core.tid], h1.b->core.pos, bam_endpos(h1.b));
            if (reg != last_reg) {
                // Move to a new region //
                while (reg == 0) {
                    sam_downsample_next(fp_in1, iter1, &h1, down1, hout->n_targets);
                    if (h1.pos == IS_EMPTY)
                        break;
                    reg = bed_get_reg(bed, hout->target_name[h1.b->core.tid], h1.b->core.pos, bam_endpos(h1.b));
                }
                if (h1.pos == IS_EMPTY) break;
                // Reads of different lengths or clipping might cause reg == last_reg //
                if (last_reg != reg) {
                    last_reg = reg;
                    // update bam2 //
                    reg = bed_get_reg(bed, hout->target_name[h2.b->core.tid], h2.b->core.pos, bam_endpos(h2.b));
                    while (heap2_lt(h2, h1) && reg != last_reg) {
                        sam_downsample_next(fp_in2, iter2, &h2, down2, hout->n_targets);
                        if (h2.pos == IS_EMPTY)
                            break;
                        reg = bed_get_reg(bed, hout->target_name[h2.b->core.tid], h2.b->core.pos, bam_endpos(h2.b));
                    }
                    fraction = (float)(drand48() * diff_frac + min_frac); // fraction of read from bam1
                    if (fp_bed_out)
                        fprintf(fp_bed_out, "%s\t%" PRIu32 "\t%" PRIu32 "\t%.3f\n", hout->target_name[h1.b->core.tid], (uint32_t)(last_reg >> 32), ((uint32_t)last_reg), fraction);
                    continue; // draw again for new region
                }
            }
            sam_write1(fp_out, hout, h1.b);
            // move bam2, if necessary
            while (heap2_lt(h2, h1)) {
                sam_downsample_next(fp_in2, iter2, &h2, down2, hout->n_targets);
            }
            // move bam1
            sam_downsample_next(fp_in1, iter1, &h1, down1, hout->n_targets);
        }
    }

    free(fn_hdr);
    free(fn_out);
    free(fn_in1);
    free(fn_in2);
    free(fn_bed);
    free(fn_bed_out);
    bed_destroy(bed);
    hts_close(fp_out);
    hts_close(fp_in1);
    hts_close(fp_in2);
    fclose(fp_bed_out);
    bam_hdr_destroy(hout);
    hts_itr_destroy(iter1);
    hts_itr_destroy(iter2);
    bam_destroy1(h1.b);
    bam_destroy1(h2.b);
    return 0;
}
