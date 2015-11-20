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
inline uint64_t sam_downsample_next(samFile *fp, hts_itr_t *iter, bam1_t *b, float down)
{
    uint64_t ret;
    do {
        if (sam_itr_next(fp, iter, b) >= 0) {
            ret = ((uint64_t)b->core.tid<<32) | (uint32_t)((uint32_t)b->core.pos+1)<<1 | bam_is_rev(b);
        } else {
            ret = IS_EMPTY;
            break;
        }
    } while (drand48() > down);
    return ret;
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
    bam1_t *b1 = 0, *b2 = 0;
    uint64_t pos1 = 0, pos2 = 0, reg = 0, last_reg = 0;
    FILE *fp_bed_out = 0;

    while ((c = getopt(argc, argv, "h:o:m:a:s:d:e:c:")) >= 0) {
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

    fprintf(stderr, "1\n");
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
    fprintf(stderr, "%s\n", mode);
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

    fprintf(stderr, "2\n");

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
    fprintf(stderr, "3\n");
    fprintf(stderr, "hout has %" PRId32 " targets\n", hout->n_targets);

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

    fprintf(stderr, "iter1 has tid %d\n", iter1->tid);
    fprintf(stderr, "iter2 has tid %d\n", iter2->tid);

    hts_idx_destroy(idx);
    fprintf(stderr, "4\n");
    b1 = bam_init1();
    b2 = bam_init1();
    if ((c = sam_itr_next(fp_in1, iter1, b1)) >= 0) {
        pos1 = ((uint64_t)b1->core.tid<<32) | (uint32_t)((uint32_t)b1->core.pos+1)<<1 | bam_is_rev(b1);
        fprintf(stderr, "Read bam\n");
    } else {
        fprintf(stderr, "Did not read bam\n");
        pos1 = IS_EMPTY;
    }
    if ((c = sam_itr_next(fp_in2, iter2, b2)) >= 0) {
        pos2 = ((uint64_t)b2->core.tid<<32) | (uint32_t)((uint32_t)b2->core.pos+1)<<1 | bam_is_rev(b2);
        fprintf(stderr, "Read bam\n");
    } else {
        pos2 = IS_EMPTY;
        fprintf(stderr, "Did not read bam\n");
    }
    // Merge data //
    fprintf(stderr, "5\n");
    diff_frac = max_frac - min_frac;
    //fraction = (float)(drand48() * diff_frac + min_frac); // fraction of reads from bam1
    while (pos1 != IS_EMPTY || pos2 != IS_EMPTY) {
        if (drand48() > fraction) { // draw from bam2
            reg = bed_get_reg(bed, hout->target_name[b2->core.tid], b2->core.pos, bam_endpos(b2));
            if (reg != last_reg) {
                // Move to a new region //
                while (reg == UINT64_MAX) {
                    if ((pos2 = sam_downsample_next(fp_in2, iter2, b2, down2)) == IS_EMPTY) break;
                    reg = bed_get_reg(bed, hout->target_name[b2->core.tid], b2->core.pos, bam_endpos(b2));
                }
                last_reg = reg;
                // update bam1 //
                reg = bed_get_reg(bed, hout->target_name[b1->core.tid], b1->core.pos, bam_endpos(b1));
                while (pos1 < pos2 && reg != last_reg) {
                    if ((pos1 = sam_downsample_next(fp_in1, iter1, b1, down1)) == IS_EMPTY) break;
                    reg = bed_get_reg(bed, hout->target_name[b1->core.tid], b1->core.pos, bam_endpos(b1));
                }
                fraction = (float)(drand48() * diff_frac + min_frac); // fraction of read from bam1
                if (fp_bed_out)
                    fprintf(fp_bed_out, "%s\t%" PRIu32 "\t%" PRIu32 "\t%.3f", hout->target_name[b2->core.tid], (uint32_t)(last_reg >> 32), ((uint32_t)last_reg), fraction);
                continue; // draw again for new region
            }
            sam_write1(fp_out, hout, b2);
            // move bam1, if necessary
            while (pos1 < pos2) {
                pos1 = sam_downsample_next(fp_in1, iter1, b1, down1);
            }
            // move bam2
            pos2 = sam_downsample_next(fp_in2, iter2, b2, down2);
        } else { // draw from bam1
            reg = bed_get_reg(bed, hout->target_name[b1->core.tid], b1->core.pos, bam_endpos(b1));
            if (reg != last_reg) {
                // Move to a new region //
                while (reg == UINT64_MAX) {
                    if ((pos1 = sam_downsample_next(fp_in1, iter1, b1, down1)) == IS_EMPTY) break;
                    reg = bed_get_reg(bed, hout->target_name[b1->core.tid], b1->core.pos, bam_endpos(b1));
                }
                last_reg = reg;
                // update bam2 //
                reg = bed_get_reg(bed, hout->target_name[b2->core.tid], b2->core.pos, bam_endpos(b2));
                while (pos2 < pos1 && reg != last_reg) {
                    if ((pos2 = sam_downsample_next(fp_in2, iter2, b2, down2)) == IS_EMPTY) break;
                    reg = bed_get_reg(bed, hout->target_name[b2->core.tid], b2->core.pos, bam_endpos(b2));
                }
                fraction = (float)(drand48() * diff_frac + min_frac); // fraction of read from bam1
                if (fp_bed_out)
                    fprintf(fp_bed_out, "%s\t%" PRIu32 "\t%" PRIu32 "\t%.3f", hout->target_name[b1->core.tid], (uint32_t)(last_reg >> 32), ((uint32_t)last_reg), fraction);
                continue; // draw again for new region
            }
            sam_write1(fp_out, hout, b1);
            // move bam2, if necessary
            while (pos2 < pos1) {
                pos2 = sam_downsample_next(fp_in2, iter2, b2, down2);
            }
            // move bam1
            pos1 = sam_downsample_next(fp_in1, iter1, b1, down1);
        }
    }

    fprintf(stderr, "fraction is %f\n", fraction);
    fprintf(stderr, "min_frac is %f\n", min_frac);
    fprintf(stderr, "max_frac is %f\n", max_frac);
    fprintf(stderr, "down1 is %f\n", down1);
    fprintf(stderr,"down2 is %f\n", down2);
    fprintf(stderr, "seed is %ld\n", random_seed);
    fprintf(stderr, "fn1 is %s\n", fn_in1);
    fprintf(stderr, "fn2 is %s\n", fn_in2);
    fprintf(stderr, "fn_bed is %s\n", fn_bed);
    fprintf(stderr, "fn_hdr is %s\n", fn_hdr? fn_hdr : "not specified");
    fprintf(stderr, "fn_out is %s\n", fn_out? fn_out : "not specified");

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
    bam_destroy1(b1);
    bam_destroy1(b2);
    return 0;
}
