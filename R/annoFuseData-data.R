# FilteredFusionAnnoFuse.tsv ----------------------------------------------
#' FilteredFusionAnnoFuse.tsv
#'
#' The exemplary output from the annoFuse workflow, containing a list of
#' filtered fusions, stored as tab-separated text file.
#' 22101 observations are included, with details about 25 different features (
#' see details below).
#'
#' This file is obtained from Star-Fusion, arriba and the rsem expression (fpkm)
#' files from release-v16-20200320 of
#' [OpenPBTA](https://github.com/AlexsLemonade/OpenPBTA-analysis)
#'
#' @format When read in via read.delim, the following information is available:
#' - LeftBreakpoint: Coordinates of breakpoint in Gene1A( if intergenic fusion the breakpoint here refers to breakpoint assigned to the nearest genes Gene1A/Gene2A)
#' - RightBreakpoint: Coordinates of breakpoint in Gene1B( if intergenic fusion the breakpoint here refers to breakpoint assigned to the nearest genes Gene1B/Gene2B)
#' - FusionName: Name of the fusion
#' - Sample: Sample ids for the input fusion calls
#' - Caller: Caller names for the input fusion calls
#' - Fusion_Type: Caller provided annotation for 3` gene in fusion is fused either as in-frame, frameshift or other
#' - JunctionReadCount: Number of supporting split read supporting the site of putative fusion
#' - SpanningFragCount: Number of spanning fragments suporting the site
#' - Confidence: Caller provided confidence of the putative fusion
#' - annots: Annotation from FusionAnnotator
#' - Gene1A: 5` end of fusion, if breakpoint is intergenic this column provides closest upstream gene
#' - Gene2A: 5` end of fusion, if breakpoint is intergenic this column provides the closest downstream gene
#' - Gene1B: 3` end of fusion, if breakpoint is intergenic this column provides closest upstream gene
#' - Gene2B: 3` end of fusion, if breakpoint is intergenic this column provides the closest downstream gene
#' - BreakpointLocation: Contains annotation of the location of breakpoint as intragenic or intergenic
#' - SpanningDelta: Provides difference between number of spanning fragments suporting the site mius number of supporting split read supporting the site
#' - reciprocal_exists: Annotates if a reciprocal fusion is also called in the input fusion calls
#' - Gene1A_anno: Annotates if Gene1A is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Gene1B_anno: Annotates if Gene1B is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Gene2A_anno: Annotates if Gene2A is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Gene2B_anno: Annotates if Gene2B is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Fusion_anno: Annotates if FusionName is in TCGA fusion calls
#' - integrated_diagnosis: Denotes the standardized 2016 WHO diagnosis
#' - broad_histology: Denotes the broad 2016 WHO classification for each tumor
#' - Kids_First_Participant_ID: Patient ids for the input fusion calls
#'
#' @name FilteredFusionAnnoFuse
#' @docType data
#'
#' @keywords datasets
NULL


# PutativeDriverAnnoFuse.tsv ----------------------------------------------
#' PutativeDriverAnnoFuse.tsv
#'
#' detailed desc
#' 4762 observations are included, with details about 29 different features (
#' see details below).
#'
#' @format When read in via read.delim, the following information is available:
#' - Sample: Sample ids for the input fusion calls
#' - LeftBreakpoint: Coordinates of breakpoint in Gene1A( if intergenic fusion the breakpoint here refers to breakpoint assigned to the nearest genes Gene1A/Gene2A)
#' - RightBreakpoint: Coordinates of breakpoint in Gene1B( if intergenic fusion the breakpoint here refers to breakpoint assigned to the nearest genes Gene1B/Gene2B)
#' - FusionName: Name of the fusion
#' - Caller: Caller names for the input fusion calls
#' - Fusion_Type: Caller provided annotation for 3` gene in fusion is fused either as in-frame, frameshift or other
#' - JunctionReadCount: Number of supporting split read supporting the site of putative fusion
#' - SpanningFragCount: Number of spanning fragments suporting the site
#' - Confidence: Caller provided confidence of the putative fusion
#' - annots: Annotation from FusionAnnotator
#' - Gene1A: 5` end of fusion, if breakpoint is intergenic this column provides closest upstream gene
#' - Gene2A: 5` end of fusion, if breakpoint is intergenic this column provides the closest downstream gene
#' - Gene1B: 3` end of fusion, if breakpoint is intergenic this column provides closest upstream gene
#' - Gene2B: 3` end of fusion, if breakpoint is intergenic this column provides the closest downstream gene
#' - BreakpointLocation: Contains annotation of the location of breakpoint as intragenic or intergenic
#' - SpanningDelta: Provides difference between number of spanning fragments suporting the site mius number of supporting split read supporting the site
#' - reciprocal_exists: Annotates if a reciprocal fusion is also called in the input fusion calls
#' - Gene1A_anno: Annotates if Gene1A is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Gene1B_anno: Annotates if Gene1B is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Gene2A_anno: Annotates if Gene2A is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Gene2B_anno: Annotates if Gene2B is TumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or/and Kinase
#' - Fusion_anno: Annotates if FusionName is in TCGA fusion calls
#' - kinaseDomainRetainedGene1A: Annotates if a kinase domain is retained in fused Gene1A
#' - kinaseDomainRetainedGene1B: Annotates if a kinase domain is retained in fused Gene1B
#' - CalledBy: Annotates all the callers that called the fusion
#' - caller_count: Number of callers that called the fusion
#' - integrated_diagnosis: Denotes the standardized 2016 WHO diagnosis
#' - broad_histology: Denotes the broad 2016 WHO classification for each tumor
#' - Kids_First_Participant_ID: Patient ids for the input fusion calls
#'
#' @references if any
#'
#' @name PutativeDriverAnnoFuse
#' @docType data
#'
#'
#' @keywords datasets
NULL

#' arriba_example.tsv
#'
#' detailed desc: Arriba fusion calls from v16 OpenPBTA
#' 25 rows are included, with details about 27 different features (see details below).
#'
#' @format When read in via read.delim(...,check.names=FALSE,stringsAsFactors = F), the following information is available:
#' - #gene1: 5` end of fusion, if breakpoint is intergenic downstream and upstream genes are provided as a comma separated
#' - gene2: 3` end of fusion,if breakpoint is intergenic downstream and upstream genes are provided as a comma separated
#' - strand1(gene/fusion): The strand before the slash reflects the strand of the gene according to the gene annotation supplied to arriba and the value after the slash reflects the strand that is transcribed. If breakpoint is intergenic the value is ".".
#' - strand2(gene/fusion): The strand before the slash reflects the strand of the gene according to the gene annotation supplied to arriba and the value after the slash reflects the strand that is transcribed.If breakpoint is intergenic the value is ".".
#' - breakpoint1: Genomic coordinate of 5` breakpoint
#' - breakpoint2: Genomic coordinate of 3` breakpoint
#' - site1: Provides information of location of 5` breakpoint values can be 5' UTR, 3' UTR, UTR (overlapping with a 5' UTR as well as a 3' UTR), CDS (coding sequence), exon, intron, and intergenic
#' - site2: Provides information of location of 5` breakpoint values can be 5' UTR, 3' UTR, UTR (overlapping with a 3' UTR as well as a 3' UTR), CDS (coding sequence), exon, intron, and intergenic
#' - type: Annotates the predicted event from the orientation of the reads supporting the fusion
#' - direction1: Denotes the orientation of the fusion
#' - direction2: Denotes the orientation of the fusion
#' - split_reads1: The number of supporting split fragments with an anchor in gene1
#' - split_reads2: The number of supporting split fragments with an anchor in gene2
#' - discordant_mates: Number of spanning fragments suporting the site
#' - coverage1: Denotes coverage over breakpoint1
#' - coverage2: Denotes coverage over breakpoint2
#' - confidence: Confidence annotation from arriba
#' - closest_genomic_breakpoint1: Only applicable when arriba is run with matched WGS sample
#' - closest_genomic_breakpoint2:  Only applicable when arriba is run with matched WGS sample
#' - filters: Denotes filters used to removed supporting reads
#' - fusion_transcript: Denotes the fusion transcript sequence
#' - reading_frame: Denotes if the 3` gene in fusion is fused either as in-frame, out-of-frame
#' - peptide_sequence: Denotes the peptide sequence
#' - read_identifiers: Comma separated names of reads that support the site
#' - tumor_id: Sample ids for the input fusion calls
#' - gene1..gene2: Name of the fusion
#' - annots: Annottaion from FusionAnnotator
#'
#' @references Arriba documentation: https://arriba.readthedocs.io/en/latest/output-files/
#' @name arriba_example
#' @docType data
#'
#' @keywords datasets
NULL


# example.rsem.genes.results.gz -------------------------------------------
#' example.rsem.genes.results.gz
#'
#' detailed desc
#' 58347 observations are included, with details about 7 different features (
#' see details below).
#'
#' @format When read in via read.delim, the following information is available:
#' - gene_id: Ensemble gene id and symbol
#' - transcript_id: Ensemble transcript id
#' - length: Denotes gene length
#' - effective_length: Denotes mappable gene regionlength
#' - expected_count: Number of reads mapping to transcripts of gene_id
#' - TPM: Transcripts Per Kilobase Million per gene_id
#' - FPKM: Fragments Per Kilobase Million per gene_id
#'
#'
#' @name example.rsem.genes.results
#' @docType data
#'
#' @keywords datasets
NULL


# exonsToPlot.RDS ---------------------------------------------------------
#' exonsToPlot.RDS
#'
#' detailed desc
#' 58347 observations are included, with details about 7 different features (
#' see details below).
#'
#' @format When read in via readRDS, the following information is available:
#' - contig: Denotes chromosome name
#' - type: Denotes if coordinates correspond to gene or exon
#' - start: Denotes genomic start position
#' - end: Denotes genomic end position
#' - strand: Denotes strand the feature is on
#' - attributes: Feature attributes in from gtf file
#' - geneName: Denotes gene symbol
#' - geneID: Denotes Ensemble gene id
#' - transcript: Denotes Ensemble transcript id
#' - exonNumber: Denotes exon rank in transcript
#' @references Feature information extracted from gencode.v27.primary_assembly.annotation.gtf.gz
#'
#' @name exonsToPlot
#' @docType data
#'
#' @keywords datasets
NULL


# expr_collapsed.rds ------------------------------------------------------
#' expr_collapsed.rds
#'
#' detailed desc
#' 58581 * 6 columns
#'
#' @references if any
#'
#' @format When read in via readRDS, the following information is available:
#' - gene_id: Denotes gene symbol
#' - BS_W97QQYKQ: Denotes sample id for example runs
#' @name expr_collapsed
#' @docType data
#'
#' @keywords datasets
NULL


# final_fusion_mmc2.txt ---------------------------------------------------
#' final_fusion_mmc2.txt
#'
#' detailed desc
#' 25664 * 7cols
#'
#' @references if any
#' @format When read in via read.delim, the following information is available:
#' - Cancer: Denotes TCGA cancer types
#' - Sample: Denotes TCGA sample id
#' - Fusion: Denotes TCGA fusion call
#' - Junction: Denotes Number of supporting split read supporting the site of putative fusion
#' - Spanning: Number of spanning fragments suporting the site
#' - Breakpoint1: Genomic coordinate of 5` breakpoint
#' - Breakpoint2: Genomic coordinate of 3` breakpoint
#'
#' @name final_fusion_mmc2
#' @docType data
#'
#' @keywords datasets
NULL


# fusionreference.txt -----------------------------------------------------
#' fusionreference.txt
#'
#' detailed desc
#' 17754 * 3
#'
#' @format When read in via read.delim, the following information is available:
#' - FusionName: Denotes Fusion name
#' - file: Source file used for annotation
#' - type: Denotes annotation type; values are TCGAFusion
#'
#' @name fusionreference
#' @docType data
#'
#' @keywords datasets
NULL


# genelistreference.txt ---------------------------------------------------
#' genelistreference.txt
#'
#' detailed desc
#' 4954 * 3 columns
#'
#' @format When read in via read.delim, the following information is available:
#' - Gene_Symbol: Denotes gene symbol
#' - type: Denotes annotation type, values could beTumorSuppressorGene, CosmicCensus, TranscriptionFactor.Oncogene or Kinase
#' - file: Source file used to gather annotation
#' @references Annotation	File	Source:
#' UCSC pfamID Description database:	http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz
#' UCSC pfamID Description database: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ucscGenePfam.txt.gz
#' TumorFusions: an integrative resource for cancer-associated transcript fusions PMID: 29099951	https://tumorfusions.org/PanCanFusV2/downloads/pancanfus.txt.gz
#' Oncogenes	http://www.bushmanlab.org/assets/doc/allOnco_Feb2017.tsv	www.bushmanlab.org
#' Tumor suppressor genes (TSGs)	https://bioinfo.uth.edu/TSGene/Human_TSGs.txt?csrt=5027697123997809089	Tumor Suppressor Gene Database 2.0 PMIDs: 23066107, 26590405
#' Kinases from The protein kinase complement of the human genome PMID: 12471243	http://kinase.com/human/kinome/tables/Kincat_Hsap.08.02.xls
#' TranscriptionFactor from https://www.sciencedirect.com/science/article/pii/S0092867418301065?via%3Dihub#mmc2 Table S1
#' COSMIC genes	https://cancer.sanger.ac.uk/census	Catalogue of Somatic Mutations in Cancer
#' Pediatric-specific oncogenes	MYBL1, SNCAIP, FOXR2, TTYH1, TERT	doi:10.1073/pnas.1300252110, doi:10.1038/nature11327, doi:10.1016/j.cell.2016.01.015, doi:10.1038/ng.2849, doi:10.1038/ng.3438, doi:10.1002/gcc.22110, doi:10.1016/j.canlet.2014.11.057, doi:10.1007/s11910-017-0722-5
#' Pediatric-specific TSGs	BCOR, QKI	doi:10.1016/j.cell.2016.01.015, doi:10.1038/ng.3500
#'
#' @name genelistreference
#' @docType data
#'
#' @keywords datasets
NULL


# gtex_collapsed.rds ------------------------------------------------------
#' gtex_collapsed.rds
#'
#' detailed desc
#' 58581*6
#' @format When read in via read.delim, the following information is available:
#' - gene_id: Denotes gene symbol
#' - SRR1069188: FPKM per gene for rna sample GTEX-XLM4-2926-SM-4AT59
#' - SRR1070986: FPKM per gene for rna sample GTEX-WZTO-0011-R6B-SM-4E3J6
#' - SRR1071289: FPKM per gene for rna sample GTEX-X4XY-0011-R8A-SM-46MVC
#' - SRR1071880: FPKM per gene for rna sample GTEX-X261-0011-R8A-SM-4E3I5
#' - SRR1072178: FPKM per gene for rna sample GTEX-XMD1-0011-R5A-SM-4AT47
#' @references GTEX: https://www.gtexportal.org/home/
#'
#' @name gtex_collapsed
#' @docType data
#'
#' @keywords datasets
NULL


# merged_arriba_tcga.tsv --------------------------------------------------
#' merged_arriba_tcga.tsv
#'
#' detailed desc
#' 8983 * 28 columns
#'
#' @format When read in via read.delim(...,check.names=FALSE,stringsAsFactors = F), the following information is available:
#' - #gene1: 5` end of fusion, if breakpoint is intergenic downstream and upstream genes are provided as a comma separated
#' - gene2: 3` end of fusion,if breakpoint is intergenic downstream and upstream genes are provided as a comma separated
#' - strand1(gene/fusion): The strand before the slash reflects the strand of the gene according to the gene annotation supplied to arriba and the value after the slash reflects the strand that is transcribed. If breakpoint is intergenic the value is ".".
#' - strand2(gene/fusion): The strand before the slash reflects the strand of the gene according to the gene annotation supplied to arriba and the value after the slash reflects the strand that is transcribed.If breakpoint is intergenic the value is ".".
#' - breakpoint1: Genomic coordinate of 5` breakpoint
#' - breakpoint2: Genomic coordinate of 3` breakpoint
#' - site1: Provides information of location of 5` breakpoint values can be 5' UTR, 3' UTR, UTR (overlapping with a 5' UTR as well as a 3' UTR), CDS (coding sequence), exon, intron, and intergenic
#' - site2: Provides information of location of 5` breakpoint values can be 5' UTR, 3' UTR, UTR (overlapping with a 3' UTR as well as a 3' UTR), CDS (coding sequence), exon, intron, and intergenic
#' - type: Annotates the predicted event from the orientation of the reads supporting the fusion
#' - direction1: Denotes the orientation of the fusion
#' - direction2: Denotes the orientation of the fusion
#' - split_reads1: The number of supporting split fragments with an anchor in gene1
#' - split_reads2: The number of supporting split fragments with an anchor in gene2
#' - discordant_mates: Number of spanning fragments suporting the site
#' - coverage1: Denotes coverage over breakpoint1
#' - coverage2: Denotes coverage over breakpoint2
#' - confidence: Confidence annotation from arriba
#' - closest_genomic_breakpoint1: Only applicable when arriba is run with matched WGS sample
#' - closest_genomic_breakpoint2:  Only applicable when arriba is run with matched WGS sample
#' - filters: Denotes filters used to removed supporting reads
#' - fusion_transcript: Denotes the fusion transcript sequence
#' - reading_frame: Denotes if the 3` gene in fusion is fused either as in-frame, out-of-frame
#' - peptide_sequence: Denotes the peptide sequence
#' - read_identifiers: Comma separated names of reads that support the site
#' - tumor_id: Sample ids for the input fusion calls
#' - gene1..gene2: Name of the fusion
#' - annots: Annottaion from FusionAnnotator
#' @references Arriba documentation: https://arriba.readthedocs.io/en/latest/output-files/.TCGA sample information from supplementary file from Driver Fusions and Their Implications in the Development and Treatment of Human Cancers available at https://www.sciencedirect.com/science/article/pii/S2211124718303954
#'
#' @name merged_arriba_tcga
#' @docType data
#'
#' @keywords datasets
NULL


# merged_starfusion_tcga.tsv ----------------------------------------------
#' merged_starfusion_tcga.tsv
#'
#' detailed desc
#' 1527 * 27
#'
#' @format When read in via read.delim(...,check.names=FALSE,stringsAsFactors = F), the following information is available:
#' - #FusionName: Denotes Fusion name
#' - JunctionReadCount: Number of supporting split read supporting the site of putative fusion
#' - SpanningFragCount: Number of spanning fragments suporting the site
#' - SpliceType: Denotes whether breakpoint is located on an exon junction
#' - LeftGene: 5` end of fusion
#' - LeftBreakpoint: Genomic coordinate of 5` breakpoint
#' - RightGene: 3` end of fusion
#' - RightBreakpoint: Genomic coordinate of 3` breakpoint
#' - LargeAnchorSupport: Denotes whether there are split reads that provide long anchor supporting read on both sides of breakpoint
#' - FFPM; Denotes exprression estimate
#' - LeftBreakDinuc: Denotes 2 bases flanking to the left of the breakpoint
#' - LeftBreakEntropy: Denotes Shannon Entropy for 15 exon bases flanking to the left of the breakpoint
#' - RightBreakDinuc: Denotes 2 bases flanking to the right of the breakpoint
#' - RightBreakEntropy:Denotes Shannon Entropy for 15 exon bases flanking to the right of the breakpoint
#' - annots: Annotation from FusionAnnotator
#' - CDS_LEFT_ID: Denotes cds_id flanking to the left of the breakpoint
#' - CDS_LEFT_RANGE: Denotes cds range to the left of the breakpoint
#' - CDS_RIGHT_ID: Denotes cds_id flanking to the right of the breakpoint
#' - CDS_RIGHT_RANGE: Denotes cds range to the right of the breakpoint
#' - PROT_FUSION_TYPE: Denotes if the 3` gene in fusion is fused either as INFRAME, FRAMESHIFT or other
#' - FUSION_MODEL: Denotes gene structure of the reconstructed fusion
#' - FUSION_CDS: Denotes cds sequence of the reconstructed fusion
#' - FUSION_TRANSL: Denotes translated sequence of the reconstructed fusion
#' - PFAM_LEFT: Denotes pfam id and location for reconstructed fusion transcript left to the breakpoint
#' - PFAM_RIGHT: Denotes pfam id and location for reconstructed fusion transcript right to the breakpoint
#' - Sample: Sample ids for the input fusion calls
#' - Caller
#' @references Star-Fusion documentation: https://github.com/STAR-Fusion/STAR-Fusion/wiki. TCGA sample information from supplementary file from Driver Fusions and Their Implications in the Development and Treatment of Human Cancers available at https://www.sciencedirect.com/science/article/pii/S2211124718303954
#'
#' @name merged_starfusion_tcga
#' @docType data
#'
#' @keywords datasets
NULL


# pfamDataBioMart.RDS -----------------------------------------------------
#' pfamDataBioMart.RDS
#'
#' detailed desc
#' 57406 * 11
#'
#' @format When read in via readRDS, the following information is available:
#' - hgnc_symbol: Denotes gene symbol
#' - pfam_id: Denotes pfam id
#' - chromosome_name: Chromosome name for the gene region
#' - gene_start: Genomic start coordinate of gene
#' - gene_end: Genomic end coordinate of gene
#' - strand: Strand information of the gene
#' - NAME: Pfam domian name
#' - DESC: Pfam domain description
#' - domain_chr: Pfam domain chromsome name
#' - domain_start: Pfam domain genomic start coordinate
#' - domain_end: Pfam domain genomic end coordinate
#' @references if any
#'
#' @name pfamDataBioMart
#' @docType data
#'
#' @keywords datasets
NULL


# sample_aliquot_tcga.tsv -------------------------------------------------
#' sample_aliquot_tcga.tsv
#'
#' detailed desc
#'
#' @format When read in via read.delim, the following information is available:
#' - name: TCGA sample name
#' - sample_id: TCGA sample id
#' - aliquot_id: TCGA aliquot id
#' @references Supplementary file from Driver Fusions and Their Implications in the Development and Treatment of Human Cancers available at https://www.sciencedirect.com/science/article/pii/S2211124718303954
#'
#' @name sample_aliquot_tcga
#' @docType data
#'
#' @keywords datasets
NULL


# starfusion_example.tsv --------------------------------------------------
#' starfusion_example.tsv
#'
#' detailed desc
#' One row of the typical output file from STAR-fusion (version TODO)
#' @format When read in via read.delim(...,check.names=FALSE,stringsAsFactors = F), the following information is available:
#' - #FusionName: Denotes Fusion name
#' - JunctionReadCount: Number of supporting split read supporting the site of putative fusion
#' - SpanningFragCount: Number of spanning fragments suporting the site
#' - SpliceType: Denotes whether breakpoint is located on an exon junction
#' - LeftGene: 5` end of fusion
#' - LeftBreakpoint: Genomic coordinate of 5` breakpoint
#' - RightGene: 3` end of fusion
#' - RightBreakpoint: Genomic coordinate of 3` breakpoint
#' - LargeAnchorSupport: Denotes whether there are split reads that provide long anchor supporting read on both sides of breakpoint
#' - FFPM: Denotes exprression estimate
#' - LeftBreakDinuc: Denotes 2 bases flanking to the left of the breakpoint
#' - LeftBreakEntropy: Denotes Shannon Entropy for 15 exon bases flanking to the left of the breakpoint
#' - RightBreakDinuc: Denotes 2 bases flanking to the right of the breakpoint
#' - RightBreakEntropy: Denotes Shannon Entropy for 15 exon bases flanking to the right of the breakpoint
#' - annots: Annotation from FusionAnnotator
#' - CDS_LEFT_ID: Denotes cds_id flanking to the left of the breakpoint
#' - CDS_LEFT_RANGE: Denotes cds range to the left of the breakpoint
#' - CDS_RIGHT_ID: Denotes cds_id flanking to the right of the breakpoint
#' - CDS_RIGHT_RANGE: Denotes cds range to the right of the breakpoint
#' - PROT_FUSION_TYPE: Denotes if the 3` gene in fusion is fused either as INFRAME, FRAMESHIFT or other
#' - FUSION_MODEL: Denotes gene structure of the reconstructed fusion
#' - FUSION_CDS: Denotes cds sequence of the reconstructed fusion
#' - FUSION_TRANSL: Denotes translated sequence of the reconstructed fusion
#' - PFAM_LEFT: Denotes pfam id and location for reconstructed fusion transcript left to the breakpoint
#' - PFAM_RIGHT: Denotes pfam id and location for reconstructed fusion transcript right to the breakpoint
#' - tumor_id: Sample ids for the input fusion calls
#' @references Star-Fusion documentation: https://github.com/STAR-Fusion/STAR-Fusion/wiki
#'
#' @name starfusion_example
#' @docType data
#'
#' @keywords datasets
NULL
