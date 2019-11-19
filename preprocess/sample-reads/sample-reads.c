#include <stdio.h>
#include <zlib.h>

#include "pcg_basic.h"
#include "kseq.h"
#include "ksort.h"

KSEQ_INIT(gzFile, gzread)
KSORT_INIT_GENERIC(uint64_t)

int main(int argc, char *argv[])
{
  gzFile fp1, fp2;
  kseq_t *seq1, *seq2;
  uint64_t total = 0, total_pairs = 0;

  if (argc < 2) {
    printf(" pass fastq1 fastq2 genome cov\n");
    return 1;
  }

  fp1 = gzopen(argv[1], "r");
  fp2 = gzopen(argv[2], "r");
  uint32_t g = atoi(argv[3]);
  uint32_t cov = atoi(argv[4]);

  seq1 = kseq_init(fp1);
  seq2 = kseq_init(fp2);
  while ((kseq_read(seq1) >= 0) && (kseq_read(seq2) >= 0 )) {
    total_pairs++;
    total += seq1->seq.l + seq2->seq.l;
  }
  
  int avg_len = total/total_pairs;
  printf(" total pairs %llu, average size %d, genome size %u\n", total_pairs, avg_len, g);
  printf("number of reads needed %u\n", (cov * g)/avg_len);
    
  kseq_destroy(seq1);
  gzclose(fp1);

  kseq_destroy(seq2);
  gzclose(fp2);

  pcg32_random_t rng;
  pcg32_srandom_r(&rng, 42u, 54u);
  

  /* Deal some reads */
  uint64_t sample_size = (cov * g) / avg_len;
  uint64_t *sample = (uint64_t *) malloc(sample_size * sizeof(uint64_t));
  uint64_t *reads = (uint64_t *) malloc(total_pairs * sizeof(uint64_t));
  uint64_t i = 0;

  for (i = 0; i < total_pairs; ++i) {
    reads[i] = i;
  }

  for (i = 0; i <  sample_size; ++i) {
    int chosen = pcg32_boundedrand_r(&rng, total_pairs - i);
    sample[i] = reads[chosen];
    reads[chosen] = reads[total_pairs - i - 1];
  }

  ks_introsort(uint64_t, sample_size, sample);

  
  fp1 = gzopen(argv[1], "r");
  char covfp1[800];
  sprintf(covfp1, "%s.cov%u.fastq", argv[1], cov);
  FILE *outfp1 = fopen(covfp1, "w");
  
  fp2 = gzopen(argv[2], "r");
  char covfp2[800];
  sprintf(covfp2, "%s.cov%u.fastq", argv[2], cov);
  FILE *outfp2 = fopen(covfp2, "w");
    
  seq1 = kseq_init(fp1);
  seq2 = kseq_init(fp2);
  int count, mark = 0;

  for (count = 0; (kseq_read(seq1) >= 0) && (kseq_read(seq2) >= 0 ); ++count) {
    if (count == sample[mark]) {
      mark++;
      // if (seq1->is_fastq)
      putc('@', outfp1);
      fprintf(outfp1, "%s", seq1->name.s);
      if (seq1->comment.l) {
	fprintf(outfp1, " %s\n", seq1->comment.s);
      }
      fprintf(outfp1, "%s\n", seq1->seq.s);
      fprintf(outfp1, "+\n");
      fprintf(outfp1, "%s\n", seq1->qual.s);

      putc('@', outfp2);
      fprintf(outfp2, "%s", seq2->name.s);
      if (seq2->comment.l) {
	fprintf(outfp2, " %s\n", seq2->comment.s);
      }
      fprintf(outfp2, "%s\n", seq2->seq.s);
      fprintf(outfp2, "+\n");
      fprintf(outfp2, "%s\n", seq2->qual.s);
      
    }
  }
    
  kseq_destroy(seq1);
  gzclose(fp1);
  fclose(outfp1);

  kseq_destroy(seq2);
  gzclose(fp2);
  fclose(outfp2);
    
  free(sample);
  free(reads);

  return 0;
}
