#!/usr/bin/env Rscript

# load dplyr for select function
library(dplyr)

# specify meta directory
meta_dir <- "novaseq/data/meta"

# specify (and create) generated_meta directory
gen_meta <- "novaseq/data/generated_meta"
if(!dir.exists(gen_meta)) dir.create(gen_meta)

# specify (and create) COI directory
coi_dir <- "novaseq/COI"
if(!dir.exists(coi_dir)) dir.create(coi_dir)

# specify (and create) COI/fastq_files directory
fastq_dir <- "novaseq/COI/fastq_files"
if(!dir.exists(fastq_dir)) dir.create(fastq_dir)


# specify directory where the novaseq files are
novaseq_files <- "/cfs/klemming/home/g/gusrutlu/TestProject/ARMS-scripts/novaseq"




# read ARMS_sequences_ARMS_EMOBON_batch3_4_5.csv that I created from ARMS_sequences_ARMS_EMOBON_batch3_4_5.xlsx
# only the first folder (ARMS sequences) was saved to csv
arms_sequences <- read.csv(file = file.path(novaseq_files, "COI_batch3.4.5.6.csv"),
                           sep = ";", header = TRUE)

# read Genoscope batch 3 file
#arms_batch3 <- read.csv(file = file.path(meta_dir, "Genoscope batch 3 COI - Sheet1.csv"),
                        #   sep = ",", header = TRUE)





# extract only koster samples
#koster_subset <- arms_sequences[grep("Koster", arms_sequences$MaterialSample.ID..EMOBON.), ]

# extract only COI samples from koster samples


# extract only columns of interest from koster_coi
arms_sequences2 <- select(arms_sequences, MaterialSampleID,
                     MaterialSample_ID_EMOBON., Sequence_batch, Code, 
                     PCR_negative_control_Code_1, PCR_negative_control_Code_2, 
                     Number_of_sequences)



# change name of MaterialSample.ID.for.DP2..ARMSMBON. to MaterialSampleID for merging columns
colnames(arms_sequences2)[colnames(arms_sequences2) == "MaterialSampleID"] <- "MaterialSampleID"

# merge the curated arms_sequences file with the arms_batch3 file by the MaterialSampleID (MaterialSample.ID.for.DP2..ARMSMBON.) column
# here, there are two sequences from batch 3.1 that are not included, and they are the sample blanks from the ARMS units processing that
# won't be included anyway

# HERE include code for extracting Batch 3 samples from wherever they are






# create summary file for sequencing batch 4
#koster_batch4 <- koster_coi[koster_coi$Sequence_batch == "Batch4", ]

# specify directory containing the sequencing files and only include those ending in .gz
seq_files <- list.files(novaseq_files, pattern = ".gz")

message("The specified directory contains the following files:")
print(seq_files)

# check if the correct files are in the seq_files directory and sort

pcr_negative1 <- c()
pcr_negative2 <- c()
project_files <- c()

for (row in 1:nrow(arms_sequences2)){
    pcr1_seqfile <- grep(arms_sequences2$PCR_negative_control_Code_1[row], seq_files)
    pcr_negative1 <- append(pcr_negative1, seq_files[pcr1_seqfile])

    pcr2_seqfile <- grep(arms_sequences2$PCR_negative_control_Code_2[row], seq_files)
    pcr_negative2 <- append(pcr_negative2, seq_files[pcr2_seqfile])

    filename_split <- strsplit(arms_sequences2$Code[row], "_")[[1]]
    filename_match <- grep(filename_split[2], seq_files)
    project_files <- append(project_files, seq_files[filename_match])

    if (identical(pcr1_seqfile, integer(0)) | identical(pcr2_seqfile, integer(0)) | identical(filename_match, integer(0)) ) {
        message(arms_sequences2$MaterialSample_ID_EMOBON.[row], " ", "does not have all associated sequence files.")
    }
}

# remove duplicates
pcr_negative1 <- unique(pcr_negative1)
pcr_negative2 <- unique(pcr_negative2)
project_files <- unique(project_files)



# create list of unused files
unused_files <- setdiff(seq_files, list.files(fastq_dir,  pattern = ".gz"))

# save to novaseq/meta as two different files
write.csv(arms_sequences2, file = file.path(gen_meta, "COI_2022_2025_LR.csv"),
          row.names = FALSE)
#write.csv(koster_batch4, file = file.path(gen_meta, "koster_coi_2022_2023.csv"),
 #         row.names = FALSE)
