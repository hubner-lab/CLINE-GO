
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop("Usage: lfmm2vcf.R <imputed.lfmm> <vcf> <blocksize> <output.vcf>")
LFMM_imp = args[1]
VCF = args[2] # not imputed vcf for restoring information
BLOCKSIZE = as.numeric(args[3])
OUT_VCF = args[4] # the rule's {output}; was gsub('.lfmm', '.vcf', LFMM_imp), unanchored (953fac)
##########################

library(dplyr)
library(data.table)
library(WGCNA)
# Convert  LFMM_imp file to VCF_imp, !! implementation only for full homozygotes !!

# Transform lfmm imputed file into vcf file (WORKAROUND)
gt <-
	fread(LFMM_imp) %>%
	transposeBigData(blocksize = BLOCKSIZE)
# Replace with VCF notation
#TODO it's solution only for diploids (should look on .geno file for polyploids)
gt[gt == 0] <- "0/0"
gt[gt == 1] <- "0/1"
gt[gt == 2] <- "1/1"

gt <- gt %>% as.data.table
# Anchored: an unanchored CHROM first matched any earlier ##contig=<ID=CHROMOSOME_I> line (5f1e90).
colnames(gt) <- fread(cmd = paste("grep -m1 '^#CHROM'",
                                      VCF,
                                      "| cut -f10- "), header = T) %>% colnames
snps <- fread(cmd = paste("grep -v \'^##\'",
			  VCF,
			  "| cut -f1-9 "))

# Check the equal number of dim. Fatal: rows are paired by position, and every
# consumer (GAPIT GD, the EMMAX TPED, emmax.R's CHROM:POS labels) trusts the pairing.
if (nrow(gt) != nrow(snps)) {
	stop(sprintf("lfmm2vcf: %d SNPs in %s but %d in %s (LEA .removed sites?) - refusing to pair genotypes with the wrong sites",
	             nrow(gt), LFMM_imp, nrow(snps), VCF))
}

# Merge snps info and GT
vcf.restored <- cbind(snps, gt)
# Save restored VCF file

# Copy meta header: every leading ## line, however many (`head -1000` cut a long header,
# e.g. one ##contig line per scaffold, at line 1000; ab75b3's class)
system(paste("sed -n '/^##/!q;p'", VCF, ">", OUT_VCF))
fwrite(vcf.restored,
       OUT_VCF,
       append = T,
       col.names = T, row.names = F, sep = '\t', quote = F)
