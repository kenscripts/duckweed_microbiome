library(tidyr)
library(plyr)
library(reshape2)

# converted and formatted biom file from q2
# used custom script biom_2_tsv

# formatted taxonomy file
# sed -i '' 's/Feature ID/feature_id/' <taxonomy file name>

###############################################################################
# Common Tasks
###############################################################################

read_my_table <- function(
                          TABLE.PATH,
                          rows = NULL
                          ) {
   MY.TABLE <- read.table(
                          TABLE.PATH,
                          sep = "\t",
                          header = T,
                          check.names = FALSE,
                          row.names = rows
                          )
   return(MY.TABLE)
}
   

###############################################################################
# Rarefaction Code
###############################################################################

get_rarefy_out <- function(
                           COUNT.TABLE,
                           subsampling = 100
                           ){
   COUNT.MATRIX <- COUNT.TABLE[,-1]
   MATRIX.TRANSPOSED <- t(COUNT.MATRIX)

   RAREMAX <- min(
                  rowSums(
                          MATRIX.TRANSPOSED
                          )
                  )

   RAREFACTION.OUT <- rarecurve(
                                MATRIX.TRANSPOSED,
                                sample = RAREMAX,
        	                label = FALSE,
                                step = subsampling
                                )

   return(RAREFACTION.OUT)
}

generate_rarefy_df <- function(
                               RAREFACTION.OUT,
                               COUNT.TABLE
                               ){
   COUNT.MATRIX <- COUNT.TABLE[,-1]

   SUBSAMPLES <- lapply(
                        RAREFACTION.OUT,
                        function(x) attr(x,"names")
                        )
   SAMPLES <- lapply(
                     seq_along(RAREFACTION.OUT),
                     function(i) rep(
                                     colnames(COUNT.MATRIX)[i],
                                     length(RAREFACTION.OUT[[i]])
                                     )
                     )
   DATA <- lapply(
                  seq_along(RAREFACTION.OUT),
                  function(i) as.vector(RAREFACTION.OUT[[i]])
                  )

   RAREFY.DF <- data.frame(
                           unlist(SAMPLES),
                           unlist(SUBSAMPLES),
                           unlist(DATA)
                           )
   colnames(RAREFY.DF) <- c(
                            "sample_id",
                            "sample_size",
                            "no_species"
                            )
   RAREFY.DF$sample_size <- gsub(
                                 "N",
                                 "",
                                 RAREFY.DF$sample_size
                                 )
   RAREFY.DF$sample_size <- as.numeric(RAREFY.DF$sample_size)

   return(RAREFY.DF)
}


##############################################################################
# Beta Diversity Analysis
##############################################################################

q2files_into_phyloseq <- function(table_tsv,tax_tsv,meta_tsv) {
	# create phyloseq otu_table
	otu_file <- read.table(table_tsv,
			       sep="\t",
			       header = T,
			       row.names=1,
			       check.names=FALSE)
	otu_matrix = as.matrix(otu_file)
	otus = otu_table(otu_matrix,
			 taxa_are_rows=TRUE)
	# create phyloseq tax_table
	tax_file <- read.table(tax_tsv,
			       sep="\t",
			       header = T,
			       row.names=1)
	tax_filtered <- tax_file[row.names(tax_file) %in% row.names(otu_file),]
	tax_filtered <- separate(tax_filtered,
   				 Taxon,
				 c("Kingdom",
				   "Phylum",
			 	   "Class",
				   "Order",
			 	   "Family",
				   "Genus",
				   "Species"
				   ),
				   sep=";",
				   remove=TRUE)
	tax_matrix = as.matrix(tax_filtered)
	tax = tax_table(tax_matrix)
	# create phyloseq sample_data
	meta_table <- read.table(meta_tsv,
	                 	sep="\t",
				header=T,
				row.names=1)
	meta = sample_data(meta_table)
	
	# create phyloseq object
	physeq = phyloseq(otus,
	 		  tax,
	 		  meta)
	print(physeq)
	return(physeq)
}



convert_q2distance <- function(DISTANCE.FILE) {
   DISTANCE.TABLE <- read.table(
                                DISTANCE.FILE,
                                sep = "\t",
                                header = T,
                                row.names = 1
                                )
   DISTANCE.MATRIX <- as.dist(
                              as.matrix(DISTANCE.TABLE)
                             )
   return(DISTANCE.MATRIX)
}

# q2distance_ordination <- ordinate(physeq,distance=q2distance_matrix,method="PCoA")
# plot_ordination(physeq,q2distance_ordination


q2distance_2_pcoa <- function(DISTANCE_FILE,META_FILE){
   library(vegan)
   DISTANCE <- read.table(DISTANCE_FILE,
                          sep="\t",
                          header=1,
                          row.names=1)
   MATRIX <- as.matrix(DISTANCE)
   MDS <- cmdscale(MATRIX,
                   eig=TRUE,
                   x.ret=TRUE)
   MDS.PER <- round(MDS$eig/sum(MDS$eig)*100,1)
   MDS.VALUES <- MDS$points
   MDS.SAMPLES <- rownames(MDS.VALUES)
   MDS.DATA <- data.frame(sample_id = rownames(MDS.VALUES),
                          x = MDS.VALUES[,1],
                          y = MDS.VALUES[,2])
   META <- read.table(META_FILE,
                      sep="\t",
                      header=T)
   MDS.META <- merge(META,
                     MDS.DATA,
                     by = "sample_id")
   MDS.OBJECTS <- list("points" = MDS.META,
                       "variance" = MDS.PER)
   return(MDS.OBJECTS)
}

##############################################################################
# Alpha Diversity Analysis
##############################################################################

addmeta_2_q2vector <- function(VECTOR.FILE,METADATA) {
   VECTOR <- read.table(VECTOR.FILE,
                        sep = "\t", 
                        header = T)
   #NAME.INDEX <- gregexpr(pattern = "vector.tsv",
   #                       VECTOR.FILE)[[1]][1]
   #NAME <- substr(VECTOR.FILE,
   #               0,
   #               NAME.INDEX-2)
   #colnames(VECTOR) <- c("sample_id",
   #                      sprintf("%s",NAME))
   colnames(VECTOR)[1]<- "sample_id"

   META <- read.table(METADATA,
                      sep = "\t", 
                      header = T)
   VECTOR.DATA <- merge(VECTOR,
                        META,
                        by = "sample_id")
   return(VECTOR.DATA)
}

##############################################################################
# Taxonomic Analysis
##############################################################################

calculateRA_4_q2table <- function(RAREFIED_TABLE_PATH) {
   # import formatted rarefied_table
   RAREFIED_TABLE <- read.table(RAREFIED_TABLE_PATH,
				sep = "\t",
				header = T,
				check.names = FALSE)
   # remove feature_id column to convert to matrix
   # calculate RA for each otu in sample
   RAREFIED_TABLE_MATRIX <- RAREFIED_TABLE[,-1]
   RA_MATRIX <- apply(RAREFIED_TABLE_MATRIX,
                      2,
		      function(x) x / sum(x) * 100)
   # add back feature_id column
   RA_TABLE <- as.data.frame(RA_MATRIX)
   RA_TABLE$feature_id <- RAREFIED_TABLE$feature_id

   return(RA_TABLE)
   }


addtax_addmeta <- function(
                           TABLE,
                           TAX.PATH,
                           META.PATH,
                           sep = FALSE,
                           value = "RA"
                           ){
   # use formatted taxonomy tsv file
   # melt normalized table and change col names
   TABLE.MELT <- melt(
                      TABLE,
                      variable.name = "sample_id",
                      value.name = paste(value)
                      )
   print(colnames(TABLE.MELT))
	
   # upload taxonomy and format
   TAX <- read.table(
                     TAX.PATH,
                     sep="\t",
                     header=T,
		     )
   colnames(TAX)[1] <- "feature_id"
   print(colnames(TAX))

   # upload metadata and format					   
   META <- read.table(
                      META.PATH,
                      sep="\t",
                      header=T,
		      )
   colnames(META)[1] <- "sample_id"
   print(colnames(META))
		
   # add taxonomy to table
   # separate Taxon column into different levels
   TABLE.TAX <- merge(
                      TABLE.MELT,
                      TAX,
                      by = "feature_id"
                      )
   if (sep == TRUE){
      TABLE.TAX <- separate(
                            TABLE.TAX,
                            Taxon,
                            c(
                              "Kingdom",
                              "Phylum",
      		              "Class",
             	              "Order",
                              "Family",
                              "Genus",
			      "Species"
			      ),
                            sep=";\\s[k,p,c,o,f,g,s]__",
                            remove = TRUE
                            )
   }
   # add metadata to table
   TABLE.INFO <- merge(
                       TABLE.TAX,
                       META,
                       by = "sample_id"
                       )
   print(colnames(TABLE.INFO))

   return(TABLE.INFO)			       			 
}

##############################################################################
# OTU Analysis
##############################################################################

# takes otu table and returns OTU list
extract_otulist <- function(FILTERED_TABLE) {
   OTU_LIST <- unique(FILTERED_TABLE$feature_id)
   OTU_LIST <- as.character(OTU_LIST)

   return(OTU_LIST)
   }


# find RA for OTUs in filtered table
addmeanRA_2_otulist <- function(OTU_LIST,FILTERED_TABLE) {
   # find mean RA of otus
   # create index counter
   i <- 1
   MEAN_RAS = c()
   for (OTU in OTU_LIST) {
      OTU_TABLE <- subset(FILTERED_TABLE,
                          feature_id == OTU)
      OTU_MEANRA <- mean(OTU_TABLE$RA)
      MEAN_RAS[[i]] <- c(OTU_MEANRA)
	
      i <- i + 1
      }

   # create dataframe from vectors
   MEAN_RAS_DF <- data.frame(OTU_LIST,
                             MEAN_RAS)
   colnames(MEAN_RAS_DF) <- c("feature_id",
                              "mean_RA")

   return(MEAN_RAS_DF)
   }

# takes otu list and filtered table 
# returns otu list with type
# this function creates a vector of character vectors
# turn a vector into character using toString
# using indexing with for loop to append character to vector
addtype_2_otulist <- function(OTU_List,Filtered_Table) {
	i <- 1
	OTU_Type = c()
	for (OTU in OTU_List) {
		otu_Table <- subset(Filtered_Table,
				    feature_id == OTU 
				    )
		Type <- toString(sort(unique(otu_Table$type)))
		OTU_Type[[i]] = c(Type)	
		i <- i + 1
	}
	OTU_Type_df <- data.frame(OTU_List,
				  OTU_Type
				  )
	colnames(OTU_Type_df) <- c("feature_id",
				   "otu_type"
				   )
	return(OTU_Type_df)
}



addtax_2_otulist <- function(OTU_LIST,TAX_TSV) {
	# read in taxonomy tsv file
	TAX_FILE <- read.table(TAX_TSV,
			       sep="\t",
			       header = T)
	colnames(TAX_FILE)[1] <- "feature_id"

	# separate taxonomy classification into levels
	TAX_SEPARATE <- separate(TAX_FILE,
   			 	 Taxon,
			 	 c("Kingdom",
			   	   "Phylum",
				   "Class",
				   "Order",
				   "Family",
				   "Genus",
				   "Species"),
		                 sep=";\\s[k,p,c,o,f,g,s]__",
			    	 remove=TRUE)
	# get taxonomy for otus
	OTU_TAX <- subset(TAX_SEPARATE,
	    	          feature_id 
			  %in%
			  OTU_LIST)
	return(OTU_TAX)
}

addinfo_2_OTUlist <- function(otu_list,filtered_table) {
	# filtered table for observed otus
	observed_otu_table <- subset(filtered_table,
				     RA > 0
				     )

	# create RA vector
	# create sample count vector
	i <- 1
	otu_meanRA = c()
	otu_samplecount = c()
	for (otu in otu_list) {
		# filter table for observed outs
		single_otu_table <- subset(observed_otu_table,
					   feature_id == otu
					   )
		# add RA of otu to vector
		mean_RA <- mean(single_otu_table$RA)
		otu_meanRA[[i]] = c(mean_RA)
		
		# add sample count to vector
		samplecount <- length(single_otu_table$sample_id)
		otu_samplecount[[i]] = c(samplecountex)
		i <- i + 1
	}

	# create dataframe from vectors
	otu_info_df <- data.frame(otu_list,
				  otu_meanRA,
				  otu_samplecount
				  )
	colnames(otu_info_df) <- c("feature_id",
				   "mean_RA",
				   "sample_count"
				   )

	return(otu_info_df)
}
###############################################################################
# SynCom Scripts
###############################################################################

get_phylum_RA <- function(TABLE){
   INDEX <- 1
   SAMPLE_LIST <- unique(TABLE$sample_id)
   PHYLA_LIST <- unique(TABLE$Phylum)
   SAMPLE_ID = c()
   SAMPLE_LOCATION = c()
   PHYLUM_NAME = c()
   PHYLUM_RA = c()
   for (SAMPLE in SAMPLE_LIST){
      SAMPLE_TABLE <- filter(TABLE,
                             sample_id == SAMPLE)
      LOCATION <- toString(unique(SAMPLE_TABLE$location))
      for (PHYLUM in PHYLA_LIST){
         PHYLUM_TABLE <- filter(SAMPLE_TABLE,
                                Phylum == PHYLUM)
         RA <- sum(PHYLUM_TABLE$RA,
                   na.rm = TRUE)
         SAMPLE_ID[[INDEX]] = c(SAMPLE)
         SAMPLE_LOCATION[[INDEX]] = c(LOCATION)
         PHYLUM_NAME[[INDEX]] = c(PHYLUM)
         PHYLUM_RA[[INDEX]] = c(RA)
         INDEX <- INDEX + 1   
      }     
   }
   PHYLUM_DF <- data.frame(SAMPLE_ID,
                           SAMPLE_LOCATION,
                           PHYLUM_NAME,
                           PHYLUM_RA)
   colnames(PHYLUM_DF) <- c("Sample",
                            "Location",
                            "Phylum",
                            "RA")
   PHYLUM_DF <- filter(PHYLUM_DF,
                       RA > 0)

   return(PHYLUM_DF)
}


group_lowRA_phyla <- function(PHYLUM_DF){
   levels(PHYLUM_DF$Phylum) <- c(levels(PHYLUM_DF$Phylum),
                                 "low_abundance (< 1% RA)")
   LOW_RA_PHYLA <- PHYLUM_DF$RA < 1
   PHYLUM_DF$Phylum[LOW_RA_PHYLA] <- "low_abundance (< 1% RA)"
   PHYLUM_DF <- ddply(PHYLUM_DF,
                      ~Sample+Location+Phylum,
                      summarise,
                      RA = sum(RA))

   return(PHYLUM_DF)
}


replace_lowphyla <- function(PHYLUM_COUNT_DF,NOZERO_RA_PHYLUMTABLE,INITIAL_TABLE){
   # input: phylum count df
   # groups all low abundant phyla and determines total otus
   PHYLA_LIST <- unique(INITIAL_TABLE$Phylum)
   HIGH_RA_PHYLA_LIST <- unique(NOZERO_RA_PHYLUMTABLE$PHYLUM_LIST)
   LOW_RA_PHYLA_LIST <- setdiff(PHYLA_LIST,
                                HIGH_RA_PHYLA_LIST)
   levels(PHYLUM_COUNT_DF$PHYLUM_NAME) <- c(levels(PHYLUM_COUNT_DF$PHYLUM_NAME),
            	                                   "low_abundance (< 1% RA)")
   for (i in 1:nrow(PHYLUM_COUNT_DF)){
      if (PHYLUM_COUNT_DF$PHYLUM_NAME[i] %in% LOW_RA_PHYLA_LIST){
         PHYLUM_COUNT_DF$PHYLUM_NAME[i] <- "low_abundance (< 1% RA)"
      }
   }
   PHYLUM_NEWCOUNT_DF <- ddply(PHYLUM_COUNT_DF,"PHYLUM_NAME",numcolwise(sum))
   return(PHYLUM_NEWCOUNT_DF)
}


get_phylum_samplecount <- function(TABLE){
   # finds the total number of samples phylum is found in
   TOTAL_SAMPLES <- length(unique(TABLE$sample_id))
   PHYLA_LIST <- unique(TABLE$Phylum)

   INDEX <- 1
   PHYLUM_NAME = c()
   PHYLUM_SAMPLE_PERCENT = c()
   for (PHYLUM in PHYLA_LIST){
      PHYLUM_TABLE <- filter(TABLE,Phylum == PHYLUM)
      NO_OF_SAMPLES <- length(unique(PHYLUM_TABLE$sample_id))
      SAMPLE_PERCENTAGE <- (NO_OF_SAMPLES / TOTAL_SAMPLES) * 100

      PHYLUM_NAME[[INDEX]] = c(PHYLUM)
      PHYLUM_SAMPLE_PERCENT[[INDEX]] = c(SAMPLE_PERCENTAGE)

      INDEX <- INDEX + 1   
   }
   PHYLUM_SAMPLECOUNT_DF <- data.frame(PHYLUM_NAME,
                                       PHYLUM_SAMPLE_PERCENT)
   colnames(PHYLUM_SAMPLECOUNT_DF) <- c("Phylum",
                                        "Sample_Percent")
   return(PHYLUM_SAMPLECOUNT_DF)
}


get_family_samplecount <- function(TABLE){
   # finds the total number of samples phylum is found in
   TOTAL_SAMPLES <- length(unique(TABLE$sample_id))
   FAMILY_LIST <- unique(TABLE$Family)

   INDEX <- 1
   FAMILY_NAME = c()
   FAMILY_SAMPLE_PERCENT = c()
   for (FAMILY in FAMILY_LIST){
      FAMILY_TABLE <- filter(TABLE,Family == FAMILY)
      NO_OF_SAMPLES <- length(unique(FAMILY_TABLE$sample_id))
      SAMPLE_PERCENTAGE <- (NO_OF_SAMPLES / TOTAL_SAMPLES) * 100
      FAMILY_NAME[[INDEX]] = c(FAMILY)
      FAMILY_SAMPLE_PERCENT[[INDEX]] = c(SAMPLE_PERCENTAGE)

      INDEX <- INDEX + 1   
   }
   FAMILY_SAMPLECOUNT_DF <- data.frame(FAMILY_NAME,
                                       FAMILY_SAMPLE_PERCENT)
   colnames(FAMILY_SAMPLECOUNT_DF) <- c("Family",
                                        "Sample_Percent")
   return(FAMILY_SAMPLECOUNT_DF)
}


get_family_RA <- function(TABLE){
   # finds RA of each family
   # create RA dataframe vectors
   INDEX <- 1
   SAMPLE_LIST <- unique(TABLE$sample_id)
   FAMILIES <- unique(TABLE$Family)
   SAMPLE_ID = c()
   SAMPLE_LOCATION = c()
   FAMILY_NAME = c()
   FAMILY_RA = c()
   for (SAMPLE in SAMPLE_LIST){
      SAMPLE_TABLE <- filter(TABLE,sample_id == SAMPLE)
      LOCATION <- toString(unique(SAMPLE_TABLE$location))
      for (FAMILY in FAMILIES){
         FAMILY_TABLE <- filter(SAMPLE_TABLE,Family == FAMILY)
         RA <- sum(FAMILY_TABLE$RA,na.rm = TRUE)

         SAMPLE_ID[[INDEX]] = c(SAMPLE)
         SAMPLE_LOCATION[[INDEX]] = c(LOCATION)
         FAMILY_NAME[[INDEX]] = c(FAMILY)
         FAMILY_RA[[INDEX]] = c(RA)
         INDEX <- INDEX + 1   
      }     
   }
   # format vectors into dataframe
   RA_FAMILYTABLE <- data.frame(SAMPLE_ID,
                                SAMPLE_LOCATION,
                                FAMILY_NAME,
                                FAMILY_RA)
   # remove phylum not observed
   RA_FAMILYTABLE <- filter(RA_FAMILYTABLE,
                            FAMILY_RA > 0)
   colnames(RA_FAMILYTABLE) <- c("Sample",
                                "Location",
                                "Family",
                                "RA")

   return(RA_FAMILYTABLE)
}


get_genus_samplecount <- function(TABLE){
   # finds the total number of samples genus is found in
   TOTAL_SAMPLES <- length(unique(TABLE$sample_id))
   GENERA <- unique(TABLE$Genus)

   INDEX <- 1
   GENUS_NAME = c()
   GENUS_SAMPLE_PERCENT = c()
   for (GENUS in GENERA){
      GENUS_TABLE <- filter(TABLE,Genus == GENUS)
      NO_OF_SAMPLES <- length(unique(GENUS_TABLE$sample_id))
      SAMPLE_PERCENTAGE <- (NO_OF_SAMPLES / TOTAL_SAMPLES) * 100

      GENUS_NAME[[INDEX]] = c(GENUS)
      GENUS_SAMPLE_PERCENT[[INDEX]] = c(SAMPLE_PERCENTAGE)

      INDEX <- INDEX + 1   
   }
   GENUS_SAMPLECOUNT_DF <- data.frame(GENUS_NAME,
                                      GENUS_SAMPLE_PERCENT)
   colnames(GENUS_SAMPLECOUNT_DF) <- c("Genus",
                                       "Sample_Percent")

   return(GENUS_SAMPLECOUNT_DF)
}


get_genus_RA <- function(TABLE){
   # finds RA of each genus
   # create RA dataframe vectors
   INDEX <- 1
   SAMPLE_LIST <- unique(TABLE$sample_id)
   GENERA <- unique(TABLE$Genus)
   SAMPLE_ID = c()
   SAMPLE_LOCATION = c()
   GENUS_NAME = c()
   GENUS_RA = c()
   for (SAMPLE in SAMPLE_LIST){
      SAMPLE_TABLE <- filter(TABLE,sample_id == SAMPLE)
      LOCATION <- toString(unique(SAMPLE_TABLE$location))
      for (GENUS in GENERA){
         GENUS_TABLE <- filter(SAMPLE_TABLE,Genus == GENUS)
         RA <- sum(GENUS_TABLE$RA,na.rm = TRUE)

         SAMPLE_ID[[INDEX]] = c(SAMPLE)
         SAMPLE_LOCATION[[INDEX]] = c(LOCATION)
         GENUS_NAME[[INDEX]] = c(GENUS)
         GENUS_RA[[INDEX]] = c(RA)
         INDEX <- INDEX + 1   
      }     
   }
   # format vectors into dataframe
   RA_GENUSTABLE <- data.frame(SAMPLE_ID,
                               SAMPLE_LOCATION,
                               GENUS_NAME,
                               GENUS_RA)
   # remove phylum not observed
   RA_GENUSTABLE <- filter(RA_GENUSTABLE,
                           GENUS_RA > 0)
   colnames(RA_GENUSTABLE) <- c("Sample",
                                "Location",
				"Genus",
                                "RA")

   return(RA_GENUSTABLE)
}

add_genus_tax <- function(
                          GENUS.INFO,
                          TAX.FILE
                          ) {
   TAX <- read.table(
                     TAX.FILE,
                     sep="\t",
                     header=T
                     )
   TAX.SEPARATE <- separate(
                            TAX,
             		    Taxon,
			    c(
                              "Kingdom",
			      "Phylum",
			      "Class",
		              "Order",
			      "Family",
			      "Genus",
			      "Species"
                              ),
			    sep = ";\\s[k,p,c,o,f,g,s]__",
			    remove = TRUE
                            )
   DROP.COLS <- c(
                  "Feature.ID",
                  "Confidence",
                  "Species"
                  )
   TAX.SEPARATE <- TAX.SEPARATE[,!colnames(TAX.SEPARATE) %in% DROP.COLS]
   TAX.SEPARATE <- distinct(TAX.SEPARATE)
   TAX.SEPARATE <- filter(
                          TAX.SEPARATE,
                          Genus %in% GENUS.INFO$Genus
                          )
   GENUS.TAX <- merge(
                      GENUS.INFO,
                      TAX.SEPARATE,
                      by = "Genus"
                      )
   return(GENUS.TAX)
}


genus_complement_matrix <- function(TABLE){
   GENUS_LIST <- unique(TABLE$Genus)
   MATRIX <- matrix(NA,
                    nrow = length(GENUS_LIST),
                    ncol = length(GENUS_LIST) + 1)
   TOTAL_SAMPLES <- length(unique(TABLE$sample_id))

   ROW_INDEX = 1
   for (GENUS in GENUS_LIST){
      COLUMN_INDEX = 1
      NO_COMPLEMENT_GENERA = 0

      GENUS_TABLE <- filter(TABLE,
                     Genus == GENUS)
      GENUS_SAMPLES <- unique(GENUS_TABLE$sample_id)
      for (OTHER_GENUS in GENUS_LIST){
         OTHER_GENUS_TABLE <- filter(TABLE,
                                     Genus == OTHER_GENUS)
         OTHER_GENUS_SAMPLES <- unique(OTHER_GENUS_TABLE$sample_id)
         RESULT <- setdiff(OTHER_GENUS_SAMPLES,
                           GENUS_SAMPLES)
         if (length(RESULT) >  0){
            MATRIX[ROW_INDEX,COLUMN_INDEX] = length(RESULT)    
            NO_COMPLEMENT_GENERA <- NO_COMPLEMENT_GENERA + 1
         } else{
            MATRIX[ROW_INDEX,COLUMN_INDEX] = 0
         }
         COLUMN_INDEX <- COLUMN_INDEX + 1
         MATRIX[ROW_INDEX,COLUMN_INDEX] = NO_COMPLEMENT_GENERA
      }
      ROW_INDEX <- ROW_INDEX + 1   
   }
   GENUS_COMPLEMENT_DF <- as.data.frame(MATRIX)
   colnames(GENUS_COMPLEMENT_DF) <- c(GENUS_LIST,"no_complement_genera")
   row.names(GENUS_COMPLEMENT_DF) <- GENUS_LIST
   return(GENUS_COMPLEMENT_DF)
}


get_feature_info <- function(TABLE) {
   TABLE <- TABLE
   FEATURES <- unique(TABLE$feature_id)
   TOTAL_SAMPLES <- length(unique(TABLE$sample_id))

   INDEX <- 1
   FEATURE_RA = c()
   FEATURE_SAMPLE_PERCENT = c()
   for (FEATURE in FEATURES) {
      FEATURE_TABLE <- filter(TABLE,
                              feature_id == FEATURE)
      MEAN_RA <- mean(FEATURE_TABLE$RA,
                      na.rm = TRUE)
      FEATURE_RA[[INDEX]] <- c(MEAN_RA)

      NO_OF_SAMPLES <- length(unique(FEATURE_TABLE$sample_id))
      SAMPLE_PERCENTAGE <- (NO_OF_SAMPLES / TOTAL_SAMPLES) * 100
      FEATURE_SAMPLE_PERCENT[[INDEX]] = c(SAMPLE_PERCENTAGE)

      INDEX <- INDEX + 1
   }
   FEATURE_INFO <- data.frame(FEATURES,
                              FEATURE_RA,
                              FEATURE_SAMPLE_PERCENT)
   colnames(FEATURE_INFO) <- c("feature_id",
                               "mean_ra",
                               "sample_percent")
   return(FEATURE_INFO)
}


add_feature_tax <- function(FEATURE_INFO,TAX_FILE) {
   TAX <- read.table(TAX_FILE,
         	     sep="\t",
		     header=T)
   TAX_SEPARATE <- separate(TAX,
             		    Taxon,
			    c("Kingdom",
			      "Phylum",
			      "Class",
		              "Order",
			      "Family",
			      "Genus",
			      "Species"),
			    sep=";\\s[k,p,c,o,f,g,s]__",
			    remove = TRUE)
   colnames(TAX_SEPARATE)[1] <- "feature_id"
   FEATURE_INFO_TAX <- merge(FEATURE_INFO,
                             TAX_SEPARATE,
                             by = "feature_id")
   return(FEATURE_INFO_TAX)
}


feature_complement_matrix <- function(TABLE,FEATURE_LIST){
   FEATURES <- as.character(FEATURE_LIST)
   MATRIX <- matrix(NA,
                    nrow = length(FEATURES),
                    ncol = length(FEATURES) + 2)

   ROW_INDEX = 1
   for (FEATURE in FEATURES){
      COLUMN_INDEX = 1
      NON_COEXISTING_FEATURES = 0

      FEATURE_TABLE <- filter(TABLE,
                              feature_id == FEATURE)
      FEATURE_SAMPLES <- unique(FEATURE_TABLE$sample_id)
      NO_OF_SAMPLES <- length(FEATURE_SAMPLES)
      for (OTHER_FEATURE in FEATURES){
         OTHER_FEATURE_TABLE <- filter(TABLE,
                                       feature_id == OTHER_FEATURE)
         OTHER_FEATURE_SAMPLES <- unique(OTHER_FEATURE_TABLE$sample_id)
         RESULT <- setdiff(OTHER_FEATURE_SAMPLES,
                           FEATURE_SAMPLES)
         if (length(RESULT) >  0){
            MATRIX[ROW_INDEX,COLUMN_INDEX] = length(RESULT)    
            NON_COEXISTING_FEATURES <- NON_COEXISTING_FEATURES + 1
         } else{
            MATRIX[ROW_INDEX,COLUMN_INDEX] = 0
         }
         COLUMN_INDEX <- COLUMN_INDEX + 1
      }
      MATRIX[ROW_INDEX,COLUMN_INDEX] = NON_COEXISTING_FEATURES
      COLUMN_INDEX <- COLUMN_INDEX + 1
      MATRIX[ROW_INDEX,COLUMN_INDEX] = NO_OF_SAMPLES

      ROW_INDEX <- ROW_INDEX + 1   
   }
   FEATURE_COMPLEMENT_DF <- as.data.frame(MATRIX)
   colnames(FEATURE_COMPLEMENT_DF) <- c(FEATURES,
                                        "non_coexisting_features",
                                        "no_of_samples")
   row.names(FEATURE_COMPLEMENT_DF) <- c(FEATURES)

   return(FEATURE_COMPLEMENT_DF)
}


feature_sample_matrix <- function(TABLE,FEATURE_LIST){
   FEATURES <- as.character(FEATURE_LIST)
   SAMPLES <- as.character(TABLE$sample_id)
   MATRIX <- matrix(NA,
                    nrow = length(FEATURES),
                    ncol = length(SAMPLES))

   ROW_INDEX = 1
   for (FEATURE in FEATURES){
      COLUMN_INDEX = 1
      FEATURE_TABLE <- filter(TABLE,
                              feature_id == FEATURE)
      FEATURE_SAMPLES <- unique(FEATURE_TABLE$sample_id)
      for (SAMPLE in SAMPLES){
         if (SAMPLE %in% FEATURE_SAMPLES){
            MATRIX[ROW_INDEX,COLUMN_INDEX] = 1
         } else{
            MATRIX[ROW_INDEX,COLUMN_INDEX] = 0
         }
         COLUMN_INDEX <- COLUMN_INDEX + 1
      }
      ROW_INDEX <- ROW_INDEX + 1   
   }
   FEATURE_SAMPLE_DF <- as.data.frame(MATRIX)
   colnames(FEATURE_SAMPLE_DF) <- c(SAMPLES)
   row.names(FEATURE_SAMPLE_DF) <- c(FEATURES)

   return(FEATURE_SAMPLE_DF)
}

location_genus_samplepercent <- function(
                                         TABLE.FILTERED,
                                         TABLE
                                         ){
   GENERA <- as.character(unique(TABLE.FILTERED$Genus))
   LOCATIONS <- as.character(unique(TABLE$location))
   MATRIX <- matrix(
                    NA,
                    nrow = length(GENERA),
                    ncol = length(LOCATIONS)
                    ) 

   LOCATION.SAMPLES.COUNT = c()
   for (LOCATION in LOCATIONS){
      LOCATION.TABLE <- filter(
                               TABLE,
                               location == LOCATION
                               )
      SAMPLE_NO <- length(unique(LOCATION.TABLE$sample_id))
      LOCATION.SAMPLES.COUNT <- c(
                                  LOCATION.SAMPLES.COUNT,
                                  SAMPLE_NO
                                  )
   } 
   print(LOCATION.SAMPLES.COUNT)

   ROW.INDEX = 1
   for (GENUS in GENERA){
      COLUMN.INDEX = 1
      LOCATION.INDEX = 1
      GENUS.TABLE <- filter(
                            TABLE,
                            Genus == GENUS
                            )
      GENUS.SAMPLES <- unique(GENUS.TABLE$sample_id)
      for (LOCATION in LOCATIONS){
         LOCATION.TABLE <- filter(
                                  GENUS.TABLE,
                                  location == LOCATION
                                  )
         SAMPLE_NO <- length(unique(LOCATION.TABLE$sample_id))
         LOCATION.SAMPLE.TOTAL <- LOCATION.SAMPLES.COUNT[[LOCATION.INDEX]]
         LOCATION.SAMPLE_PERCENT <- SAMPLE_NO / LOCATION.SAMPLE.TOTAL
         LOCATION.SAMPLE_PERCENT <- format(
                                           round(
                                                 LOCATION.SAMPLE_PERCENT,
                                                 4
                                                 ),
                                           nsmall = 4
                                           )
         LOCATION.SAMPLE_PERCENT <- as.numeric(LOCATION.SAMPLE_PERCENT)
         
         MATRIX[ROW.INDEX,COLUMN.INDEX] = LOCATION.SAMPLE_PERCENT
         COLUMN.INDEX <- COLUMN.INDEX + 1

         LOCATION.INDEX <- LOCATION.INDEX + 1
      }

      ROW.INDEX <- ROW.INDEX + 1   
   }

   GENUS.LOCATION.SAMPLE_PERCENT <- as.data.frame(MATRIX)
   colnames(GENUS.LOCATION.SAMPLE_PERCENT) <- LOCATIONS
   row.names(GENUS.LOCATION.SAMPLE_PERCENT) <- GENERA

   return(GENUS.LOCATION.SAMPLE_PERCENT)
}


feature_location_matrix <- function(DATAFRAME,TT_TABLE){
   FEATURES <- as.character(unique(DATAFRAME$feature_id))
   LOCATIONS <- as.character(unique(TT_TABLE$location))
   MATRIX <- matrix(NA,
                    nrow = length(FEATURES),
                    ncol = length(LOCATIONS)) 

   LOCATION_SAMPLES = c()
   for (LOCATION in LOCATIONS){
      LOCATION_TABLE <- filter(TT_TABLE,
                               location == LOCATION)
      NO_OF_SAMPLES <- length(unique(LOCATION_TABLE$sample_id))
      LOCATION_SAMPLES <- c(LOCATION_SAMPLES,
                            NO_OF_SAMPLES)
   } 
   print(LOCATION_SAMPLES)
   
   ROW_INDEX = 1
   for (FEATURE in FEATURES){
      COLUMN_INDEX = 1
      LOCATION_INDEX = 1
      FEATURE_TABLE <- filter(TT_TABLE,
                              feature_id == FEATURE)
      FEATURE_SAMPLES <- unique(FEATURE_TABLE$sample_id)
      for (LOCATION in LOCATIONS){
         LOCATION_TABLE <- filter(FEATURE_TABLE,
                                  location == LOCATION)
         SAMPLE_NO <- length(unique(LOCATION_TABLE$sample_id))
         TOTAL_SAMPLES <- LOCATION_SAMPLES[[LOCATION_INDEX]]
         LOCATION_PERCENT <- SAMPLE_NO / TOTAL_SAMPLES
         LOCATION_PERCENT <- format(round(LOCATION_PERCENT,4),
                                    nsmall = 4)
         LOCATION_PERCENT <- as.numeric(LOCATION_PERCENT)
         

         MATRIX[ROW_INDEX,COLUMN_INDEX] = LOCATION_PERCENT
         COLUMN_INDEX <- COLUMN_INDEX + 1

         LOCATION_INDEX <- LOCATION_INDEX + 1
      }
      ROW_INDEX <- ROW_INDEX + 1   
   }
   FEATURE_LOCATION_DF <- as.data.frame(MATRIX)
   colnames(FEATURE_LOCATION_DF) <- LOCATIONS
   row.names(FEATURE_LOCATION_DF) <- FEATURES

   return(FEATURE_LOCATION_DF)
}

write_feature_files <- function(SYNCOM_TABLE){
   SYNCOM_TABLE <- SYNCOM_TABLE
   GENUS_LIST <- unique(SYNCOM_TABLE$Genus)
   for (GENUS in GENUS_LIST){
      GENUS_TABLE <- filter(SYNCOM_TABLE,
                            Genus == GENUS)
      GENUS_FEATURES <- extract_otulist(GENUS_TABLE)
      GENUS_FILE <- sprintf("feature_files/%s_features.txt",GENUS)
      write.table(GENUS_FEATURES,
                  GENUS_FILE,
                  row.names = FALSE)
      }
}


water_tissue_RA <- function(TABLE) {
   library(dplyr)
   library(reshape)
   GENUS_LIST <- unique(TABLE$Genus)
   TYPES <- unique(TABLE$type)

   INDEX <- 1
   GENUS = c()
   TYPE = c()
   GENUS_RA = c()
   for (GENUS_NAME in GENUS_LIST) {
      GENUS_TABLE <- filter(TABLE,
                            Genus == GENUS_NAME)
      for (TYPE_NAME in TYPES){
         TYPE_TABLE <- filter(GENUS_TABLE,
                              type == TYPE_NAME )
         SAMPLES <- unique(TYPE_TABLE$sample_id)
         SAMPLES_RA <- c()
         INDEX_2 <- 1
         for (SAMPLE in SAMPLES){
            SAMPLE_TABLE <- filter(TYPE_TABLE,
                                   sample_id == SAMPLE)
            SAMPLE_RA <- sum(SAMPLE_TABLE$RA,na.rm = TRUE)
            SAMPLES_RA[[INDEX_2]] <- c(SAMPLE_RA)
            INDEX_2 <- INDEX_2 + 1
         }
         MEAN_RA <- mean(SAMPLES_RA,
                         na.rm = TRUE)
         GENUS[[INDEX]] <- c(GENUS_NAME)
         TYPE[[INDEX]] <- c(TYPE_NAME)
         GENUS_RA[[INDEX]] <- c(MEAN_RA)
         INDEX <- INDEX + 1
      }
   }
   WATER_TISSUE_DF <- data.frame(GENUS,
                                 TYPE,
                                 GENUS_RA)
   WATER_TISSUE_DF <- cast(WATER_TISSUE_DF,
                           GENUS~TYPE)
   return(WATER_TISSUE_DF)
}

location_diffgenus_count <- function(
                                     TABLE
                                     ){

   LOCATIONS <- as.character(unique(TABLE$location))
   MATRIX <- matrix(
                    NA,
                    nrow = length(LOCATIONS),
                    ncol = length(LOCATIONS)
                    ) 

   ROW.INDEX = 1
   for (LOCATION in LOCATIONS){
      COLUMN.INDEX = 1
      LOCATION.INDEX = 1
      LOCATION.1 <- filter(
                           TABLE,
                           location == LOCATION
                           )
      LOCATION.1.GENERA <- unique(LOCATION.1$Genus)
      for (LOCATION in LOCATIONS){
         LOCATION.2 <- filter(
                              TABLE,
                              location == LOCATION
                              )
         LOCATION.2.GENERA <- unique(LOCATION.2$Genus)
         COMBINED.GENERA <- c(
                              LOCATION.1.GENERA,
                              LOCATION.2.GENERA
                              )
         TOTAL.GENERA <- length(unique(COMBINED.GENERA))
         COMMON.GENERA <- intersect(
                                    LOCATION.1.GENERA,
                                    LOCATION.2.GENERA
                                    )
         DIFF.GENERA <- TOTAL.GENERA - length(COMMON.GENERA)

         print("comparison")
         print(length(LOCATION.1.GENERA))
         print(length(LOCATION.2.GENERA))
         print(TOTAL.GENERA)
         print(length(COMMON.GENERA))
         print(DIFF.GENERA)
         
         MATRIX[ROW.INDEX,COLUMN.INDEX] = DIFF.GENERA

         COLUMN.INDEX <- COLUMN.INDEX + 1
         LOCATION.INDEX <- LOCATION.INDEX + 1
      }
      ROW.INDEX <- ROW.INDEX + 1   
   }

   LOCATION.GENUS.COUNT <- as.data.frame(MATRIX)
   colnames(LOCATION.GENUS.COUNT) <- LOCATIONS
   row.names(LOCATION.GENUS.COUNT) <- LOCATIONS

   return(LOCATION.GENUS.COUNT)
}

###############################################################################
# Caldwell House and Passion Puddle Scripts
###############################################################################

find_tax_RA_between_types <- function(TABLE,TAX_LEVEl){
   INDEX <- 1
   SAMPLE_LIST <- unique(TABLE$sample_id)
   TAX_LIST <- unique(TABLE[[TAX_LEVEL]])

   SAMPLE_ID = c()
   SAMPLE_LOCATION = c()
   SAMPLE_TYPE = c()
   TAX_NAME = c()
   TAX_RA = c()
   for (SAMPLE in SAMPLE_LIST){
      SAMPLE_TABLE <- filter(TABLE,
                             sample_id == SAMPLE)
      LOCATION <- toString(unique(SAMPLE_TABLE$location))
      TYPE <- toString(unique(SAMPLE_TABLE$type))
      print(SAMPLE)
      print(LOCATION)
      print(TYPE)
      for (TAX in TAX_LIST){
         TAX_TABLE <- filter(SAMPLE_TABLE,
                             UQ(sym(TAX_LEVEL)) == TAX)
         RA <- sum(TAX_TABLE$RA,
                   na.rm = TRUE)

         SAMPLE_ID[[INDEX]] = c(SAMPLE)
         SAMPLE_LOCATION[[INDEX]] = c(LOCATION)
         SAMPLE_TYPE[[INDEX]] = c(TYPE)
         TAX_NAME[[INDEX]] = c(TAX)
         TAX_RA[[INDEX]] = c(RA)
         INDEX <- INDEX + 1   
      }     
   }
   TAX_DF <- data.frame(SAMPLE_ID,
                        SAMPLE_LOCATION,
                        SAMPLE_TYPE,
                        TAX_NAME,
                        TAX_RA)
   TAX_DF <- filter(TAX_DF,
                    TAX_RA > 0)
   return(TAX_DF)
}


genus_featurecount_intypes <- function(TABLE){
   # finds the total number of samples genus is found in
   GENERA <- unique(TABLE$Genus)
   TYPES <- unique(TABLE$type)

   INDEX <- 1
   GENUS_NAME = c()
   SAMPLE_TYPE = c()
   SAMPLE_ID = c()
   GENUS_FEATURECOUNT = c()
   for (GENUS in GENERA){
      GENUS_TABLE <- filter(TABLE,
                            Genus == GENUS)
      SAMPLES <- unique(GENUS_TABLE$sample_id)
      for (SAMPLE in SAMPLES){
           SAMPLE_TABLE <- filter(GENUS_TABLE,
                                  sample_id == SAMPLE )
           NO_OF_FEATURES <- length(unique(SAMPLE_TABLE$feature_id))
           TYPE <- toString(unique(SAMPLE_TABLE$type))

           GENUS_NAME[[INDEX]] = c(GENUS)
           SAMPLE_TYPE[[INDEX]] = c(TYPE)
           SAMPLE_ID[[INDEX]] = c(SAMPLE)
           GENUS_FEATURECOUNT[[INDEX]] = c(NO_OF_FEATURES)

           INDEX <- INDEX + 1   
      }
   }
   GENUS_FEATURECOUNT_DF <- data.frame(GENUS_NAME,
                                       SAMPLE_ID,
                                       SAMPLE_TYPE,
                                       GENUS_FEATURECOUNT)

   return(GENUS_FEATURECOUNT_DF)
}


water_tissue_RA.v2 <- function(TABLE) {
   library(dplyr)
   library(reshape)
   GENUS_LIST <- unique(TABLE$Genus)
   TYPES <- unique(TABLE$type)

   INDEX <- 1
   GENUS = c()
   TYPE = c()
   GENUS_RA = c()
   for (GENUS_NAME in GENUS_LIST) {
      GENUS_TABLE <- filter(TABLE,
                            Genus == GENUS_NAME)
      for (TYPE_NAME in TYPES){
         TYPE_TABLE <- filter(GENUS_TABLE,
                              type == TYPE_NAME )
         SUM <- sum(TYPE_TABLE$RA)
         NO_OF_SAMPLES <- length(unique(TYPE_TABLE$sample_id))
         MEAN_RA <- SUM/NO_OF_SAMPLES

         GENUS[[INDEX]] <- c(GENUS_NAME)
         TYPE[[INDEX]] <- c(TYPE_NAME)
         GENUS_RA[[INDEX]] <- c(MEAN_RA)
         INDEX <- INDEX + 1
      }
   }
   WATER_TISSUE_DF <- data.frame(GENUS,
                                 TYPE,
                                 GENUS_RA)
   WATER_TISSUE_DF <- cast(WATER_TISSUE_DF,
                           GENUS~TYPE)
   return(WATER_TISSUE_DF)
}

get_genus_info.v2 <- function(TABLE) {
   GENUS_LIST <- unique(TABLE$Genus)
   TOTAL_SAMPLES <- length(unique(TABLE$sample_id))

   INDEX <- 1
   GENUS_MEAN_RA = c()
   GENUS_MEDIAN_RA = c()
   GENUS_SAMPLE_PERCENT = c()
   GENUS_OTUS = c()
   for (GENUS in GENUS_LIST) {
      GENUS_TABLE <- filter(
                            TABLE,
                            Genus == GENUS
                            )

      RA_SUM <- sum(GENUS_TABLE$RA)
      NO_OF_SAMPLES <- length(unique(GENUS_TABLE$sample_id))
      MEAN_RA <- RA_SUM/NO_OF_SAMPLES
      GENUS_MEAN_RA[[INDEX]] <- c(MEAN_RA)

      MEDIAN_RA <- median(
                          GENUS_TABLE$RA,
                          na.rm = TRUE
                          )
      GENUS_MEDIAN_RA[[INDEX]] <- c(MEDIAN_RA)

      SAMPLE_PERCENTAGE <- (NO_OF_SAMPLES / TOTAL_SAMPLES) * 100
      GENUS_SAMPLE_PERCENT[[INDEX]] = c(SAMPLE_PERCENTAGE)
 
      NO_OF_OTUS <- length(unique(GENUS_TABLE$feature_id))
      GENUS_OTUS[[INDEX]] = c(NO_OF_OTUS)

      INDEX <- INDEX + 1
   }
   GENUS_INFO <- data.frame(
                            GENUS_LIST,
                            GENUS_MEAN_RA,
                            GENUS_MEDIAN_RA,
                            GENUS_SAMPLE_PERCENT,
                            GENUS_OTUS
                            )
   GENUS_INFO <- GENUS_INFO[order(-GENUS_INFO$GENUS_MEAN_RA),]
   GENUS_INFO$Mean_RA_Rank <- seq.int(
                                      nrow
                                          (
                                           GENUS_INFO
                                           )
                                      )
   colnames(GENUS_INFO) <- c(
                             "Genus",
                             "Mean_RA",
                             "Median_RA",
                             "Sample_Percent",
                             "OTU_Count",
                             "Mean_RA_Rank"
                             )
   return(GENUS_INFO)
}


remove_unknown_genus <- function(TABLE){
   KNOWN.GENUS <- filter(
                         TABLE,
                         Genus != ""
                         )
   KNOWN.GENUS <- KNOWN.GENUS[!is.na(KNOWN.GENUS$Genus),]
  
   return(KNOWN.GENUS)
}
   

get_legend <- function(MY.GGPLOT){
  library(gridExtra)
  TMP <- ggplot_gtable(
                       ggplot_build(
                                    MY.GGPLOT
                                    )
                       )
  LEG <- which(
               sapply(
                      TMP$grobs,
                      function(x) x$name == "guide-box"
                      )
               )
  LEGEND <- TMP$grobs[[LEG]]

  return(LEGEND)
}

addtype_2_genus <- function(
                            GENUS.VECTOR,
                            TABLE
                            ) {
   library(dplyr)
   INDEX <- 1
   GENUS.TYPE = c()
   for (GENUS in GENUS.VECTOR) {
      GENUS.TABLE <- filter(
                            TABLE,
			    Genus == GENUS
			    )
      TYPE <- toString(
                       sort(
                            unique(
                                   GENUS.TABLE$type
                                   )
                            )
                       )
      GENUS.TYPE[[INDEX]] = c(TYPE)	
      INDEX <- INDEX + 1
   }
   TYPE.DF <- data.frame(
                         GENUS.VECTOR,
                         GENUS.TYPE
			 )
   colnames(TYPE.DF) <- c(
                          "Genus",
		          "sample_type"
			  )
   return(TYPE.DF)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Differential Abundance Testing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

sort_table_samples <- function(
                               TABLE.META,
                               TABLE.COUNT
                               ){
   SAMPLES <- as.vector(TABLE.META$sample_id)
   COLS <- c(
             "feature_id",
             SAMPLES
             )
   TABLE.SAMPLES <- TABLE.COUNT[,colnames(TABLE.COUNT) %in% COLS]
   TABLE.SORTED <- TABLE.SAMPLES[COLS]
   
   return(TABLE.SORTED)
}

create_genus_count_table <- function (
                                      TABLE,
                                      TAX.PATH,
                                      unknowns = "no"
                                      ){
   TABLE.MELT <- melt(
                      TABLE,
                      variable.name = "sample_id",
                      value.name = "count"
                      )

   TAX <- read.table(
                     TAX.PATH,
                     sep="\t",
                     header=T,
		     )
   colnames(TAX)[1] <- "feature_id"
   TAX.SEPARATE <- separate(
                            TAX,
                            Taxon,
                            c(
                              "Kingdom",
                              "Phylum",
      		              "Class",
             	              "Order",
                              "Family",
                              "Genus",
   			      "Species"
      			      ),
                            sep=";\\s[k,p,c,o,f,g,s]__",
                            remove = TRUE
                            )
   TABLE.INFO <- merge(
                       TABLE.MELT,
                       TAX.SEPARATE,
                       by = "feature_id"
                       )
   if (unknowns == "no"){
      TABLE.INFO <- remove_unknown_genus(TABLE.INFO)
   }

   TABLE.SELECT <- select(
                          TABLE.INFO,
                          feature_id,
                          sample_id,
                          count,
                          Genus,
                          )
   TABLE.COLLAPSE <- ddply(
                           TABLE.SELECT,
                           ~Genus+sample_id,
                           summarise,
                           genus_count = sum(count)
                           )
   TABLE.CAST <- dcast(
                       TABLE.COLLAPSE,
                       Genus~sample_id
                       )
   rownames(TABLE.CAST) <- TABLE.CAST[,1]
   TABLE.CAST <- TABLE.CAST[,-1]
   TABLE.GENUS <- data.frame(
                             TABLE.CAST,
                             check.names = FALSE
                             )

   return(TABLE.GENUS)
}

combine_aldex_output <- function(
                                 CLR.EFFECT,
                                 CLR.TEST
                                 ){

   CLR.EFFECT$Genus <- rownames(CLR.EFFECT)
   CLR.TEST$Genus <- rownames(CLR.TEST)
   GENERA.INFO <- merge(
                        CLR.EFFECT,
                        CLR.TEST,
                        by = "Genus"
                        )
   
   return(GENERA.INFO)
}

clr_transform_counttable <- function(
                                     COUNT.TABLE,
                                     PRIOR = 0.001
                                     ) {
   # works correctly; ran on small test set
   # add prior then transform
   COUNT.PRIOR <- COUNT.TABLE + PRIOR 
   COUNT.CLR <- apply(
                      COUNT.PRIOR,
                      2,
                      function(x){
                                  log2(x)- mean(log2(x))
                                  }
                      )
   CLR.TABLE <- as.data.frame(COUNT.CLR)

   return(CLR.TABLE)
}

add_clr_info <- function(
                         CLR.TABLE,
                         TAX.PATH,
                         META.PATH
                         ){
   # melt clr table
   CLR.MELT <- melt(
                    CLR.TABLE,
                    variable.name = "sample_id",
                    value.name = "clr"
                    )
	
   # get and format tax
   TAX <- read.table(
                     TAX.PATH,
                     sep="\t",
                     header=T,
		     )
   TAX.SEPARATE <- separate(
                            TAX,
                            Taxon,
                            c(
                              "Kingdom",
                              "Phylum",
      		              "Class",
             	              "Order",
                              "Family",
                              "Genus",
   			      "Species"
      			      ),
                            sep=";\\s[k,p,c,o,f,g,s]__",
                            remove = TRUE
                            )
   TAX.SELECT <- select(
                        TAX.SEPARATE,
                        Phylum,
                        Family,
                        Genus
                        )
   TAX.DISTINCT <- distinct(TAX.SELECT)
   # get and format meta
   META <- read.table(
                      META.PATH,
                      sep="\t",
                      header=T,
		      )
   colnames(META)[1] <- "sample_id"
		
   # merge dataframes
   CLR.TAX <- merge(
                    CLR.MELT,
                    TAX.DISTINCT,
                    by = "Genus"
                    )
   CLR.INFO <- merge(
                     CLR.TAX,
                     META,
                     by = "sample_id"
                     )

   return(CLR.INFO)			       			 
}

abs_log <- function(x){
   x[x == 0] <- 1
   SIGN <- sign(x)
   LOG2 <- SIGN * log2(SIGN * x)
}
