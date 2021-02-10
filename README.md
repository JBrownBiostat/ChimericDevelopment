# ChimericDevelopment
Repository for code used perform analysis on 2020 (submitted) paper: Interspecies Chimeric Conditions Affect the Developmental Rate of Human Pluripotent Stem Cells, J. Brown, C. Barry, et al.

Gene sequencing data for this study is available from the Gene Expression Omnibus under accession number GSE157354.

"WorkDir" in the code indicates the directory which contains each of the subdirectories which are available on this repository.

If performing duplicate analysis, run the following scripts in order first:
1) PreprocessingAndNormalization/FormatRawData.R
2) PreprocessingAndNormalization/NormalizeData.R
3) Trendy/RunTrendy.R (long, consider reducing "nIter")
4) Trendy/ClassifyGenes.R
5) Trendy/goEnrichment.R
6) FACS_Data/PreprocessAndNormalize_ReAlign.R
7) FACS_Data/RunTrendy.R
8) FACS_Data/ClassifyGenes.R
9) FACS_Data/goEnrichment.R

In "BrainSpan" or "HumanProteinAtlas", similarly run the TypingAndTiming.R script first.

Some scripts (eg. in BrainSpan or HumanProteinAtlas) will require reference datasets. Please see our paper for details.
