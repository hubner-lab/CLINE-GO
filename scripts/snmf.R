library(LEA)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

# The project-wide seed. set.seed() governs R's RNG only; LEA::snmf() draws its own
# seed at random by default and never consults it, so without passing SEED through to
# snmf(seed=) the Q-matrices, cross-entropy curve (hence k_best), the sNMF-driven
# imputation and every structure-corrected p-value downstream differed between two
# runs of identical inputs (finding 8d07b1). Measured on SIMDATA: with seed=42 the
# cross-entropies are identical run to run, and the repetitions still differ from
# each other, so the multi-repetition best-run selection is unaffected.
SEED <- 42L
set.seed(SEED)

####################
GENO=args[1] ; LFMM=sub('\\.geno$', '.lfmm', GENO)
K_START=args[2] %>% as.numeric
K_END=args[3] %>% as.numeric
PLOIDY=args[4] %>% as.numeric
REPEAT=args[5] %>% as.numeric
CPU=args[6] %>% as.numeric
PROJECT=args[7] # new or continue
####################

# FUN SNMF analysis
FUN_snmf <- function(Ks,
                     Ke,
                     GENO,
		     ploidy = 2,
                     repetions = 10,
		     entropy = T,
                     project = 'new', # with force run TRUE
                     I = 10000,
                     CPU = 24,
                     seed = SEED){

  SNMF = paste0(gsub('.geno', '', GENO), '.snmfProject')

  

   message(paste0('INFO: sNMF seed = ', seed))
   project = snmf(GENO,
                    CPU = CPU,
                    K = Ks:Ke,
                    entropy = TRUE,
                    repetitions = repetions,
                    ploidy = ploidy,
                    I = I,
                    seed = seed,
                    project = project)                    

  return(project)
}


################### MAIN #########################

# Define the number of SNPs in GENO
nSNP = fread(GENO) %>% nrow
# Run SNMF
FUN_snmf(K_START,
	 K_END,
	 GENO,
	 ploidy = PLOIDY,
	 repetions = REPEAT,
	 project = PROJECT,
	 I = min(10000, nSNP),
	 CPU = CPU)
	
