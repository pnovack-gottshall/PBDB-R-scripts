## Format Paleobiology Database (PBDB) genus/subgenus occurrences into compact
## taxonomic structure (using a parallel-computing environment), including FADs
## and LADs. Once done, then checks data for homonyms, problematic duplicate
## entries, and the frequency of genera with subgenera. Written for Phil
## Novack-Gottshall, Benedictine University, Lisle, IL
## <pnovack-gottshall@ben.edu>

## 1. Download data ------------------------------------------------------------

# Download data directly from Paleobiology Database and save to working #
# directory. (Note the downloaded data file will be several tens of megabytes.)
# DO NOT TEMPORALLY RESTRICT TO PHANEROZOIC! Restricting only includes taxa with
# fossil occurrences. Here we want ALL taxa with taxonomic opinions.


## Easier if paste URL link into browser and save manually as pbdb_data.csv
getwd()
# pbdb <- read.csv("www.paleobiodb.org/data1.2/taxa/list.csv?base_name=Metazoa&show=app&vocab=pbdb")
# If want forams too, use base_name=Metazoa,Retaria
# If want forams and plants, too, use base_name=Metazoa,Retaria,Plantae
# If want the entire database, use base_name=Life
pbdb <- read.csv("pbdb_data.csv")
head(pbdb, 1)
# As a quick check, the first row should include the most inclusive taxon
# downloaded above in the taxon_name column.

# Simple example to identify possibly problematic homonym genera (or duplicates
# or multiple listings of a genus, as occurs when there are subgenera)
which.gsg <- 
  which((pbdb$accepted_rank == "genus" | pbdb$accepted_rank == "subgenus") 
        & pbdb$difference == "")
sort(table(pbdb$accepted_name[which.gsg]), decreasing = TRUE)[1:20]
# Example (as of 9/1/2019, includes a homonym and likely duplicate entry):
pbdb[which(pbdb$accepted_name == "Lowenstamia"), ]




## 2. FUNCTIONS -----------------------------------------------------------------

## prep.PBDB: Function to add higher taxonomic names (phylum, class, order,
## etc.) for PBDB genus (and subgenus) names.
# g = Vector (sequence) of number of genus names to process.
# gen.order = Vector of ordered PBDB genus (and subgenus) names.
# which.gsg = Vector of indices for PBDB entries tagged as accepted genus or 
#   subgenus names.
# pbdb = data frame of all PBDB occurrences.
#
# Output is a list, with each item the taxonomic ranks for a single genus.
# Extends LAD to 'Recent' if genus is extant and splits subgenus names into
# genus and subgenus components.
prep.PBDB <- function(g = 1, gen.order, which.gsg, pbdb) {
  scales <- c("superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", 
              "subphylum", "superclass", "class", "subclass", "infraclass", 
              "superorder", "order", "suborder", "infraorder", "superfamily", 
              "family", "subfamily", "tribe", "subtribe", "genus", "subgenus")
  out <- data.frame(Superkingdom = character(1), Kingdom = character(1), 
                    Subkingdom = character(1), Superphylum = character(1), 
                    Phylum = character(1), Subphylum = character(1), 
                    Superclass = character(1), Class = character(1), 
                    Subclass = character(1), Infraclass = character(1), 
                    Superorder = character(1), Order = character(1), 
                    Suborder = character(1), Infraorder = character(1), 
                    Superfamily = character(1), Family = character(1), 
                    Subfamily = character(1), Tribe = character(1), 
                    Subtribe = character(1), Genus = character(1), 
                    Subgenus = character(1), Species = "sp.", 
                    stringsAsFactors = FALSE)
  out$Genus <- as.character(pbdb$accepted_name[which.gsg][gen.order[g]])
  wh <- which.gsg[gen.order[g]]
  out$max_ma <- as.numeric(pbdb$firstapp_max_ma[wh])
  out$min_ma <- as.numeric(pbdb$lastapp_min_ma[wh])
  # Implement 'Pull-of-the-Recent' extension:
  if (any(pbdb$is_extant[wh] == "extant"))
    out$min_ma <- 0
  # Properly assign subgenera and genera:
  if (pbdb$accepted_rank[wh] == "subgenus") {
    split.subgenus <- strsplit(out$Genus, " ")[[1]]
    out$Genus <- as.character(split.subgenus[1])
    out$Subgenus <- as.character(gsub("[()]", "", split.subgenus[2]))
  }
  parent <- pbdb[which(pbdb$accepted_no == pbdb$parent_no[wh]), ][1, ]
  repeat {
    if (parent$accepted_rank %in% scales)
      out[1, which(scales == parent$accepted_rank)] <-
        as.character(parent$accepted_name)
    parent <-
      pbdb[which(pbdb$accepted_no == parent$parent_no),][1,]
    if (all(is.na(parent)))
      break
  }
  return(out)
}



## 3. Format the PBDB data using a parallel-computing environment ---------------

# The following only runs in-parallel in a Windows OS. If using a fork-style OS
# (i.e., Mac OS or LINUX), will need to modify. The simplest (but not necessarily
# the fastest for your OS) is to replace the sfLapply() line with mclapply().

# Version using parallel computing:
library(data.table) # Required below for merging parallel lists into dataframe
library(snowfall)
(t.start0 <- Sys.time())
# Initialize
which.gsg <- 
  which((pbdb$accepted_rank == "genus" | pbdb$accepted_rank == "subgenus") 
        & pbdb$difference == "")
gen.order <- order(pbdb$accepted_name[which.gsg])
gen.seq <- seq_along(gen.order)
# gen.seq <- 1:1000 # Use if want a faster example for first 1000 genera
# Set up computer cluster
library(parallel)
cpus <- parallel::detectCores() # Number of CPUs to cluster together
# sfSetMaxCPUs(cpus)			      # Use if plan more than 32 CPUs
sfInit(parallel = TRUE, cpus = cpus, slaveOutfile = "initfile") # Initialize cluster
stopifnot(sfCpus() == cpus)		    # Confirm set up CPUs properly
stopifnot(sfParallel() == TRUE)		# Confirm now running in parallel
sfExportAll()				            # Export all libraries, files, & objects
# Execute the function
prep <- NA
prep <- sfLapply(x = gen.seq, fun = prep.PBDB, gen.order = gen.order, 
                 which.gsg = which.gsg, pbdb = pbdb) # Version without load-balancing
sfStop()
output <- data.table::rbindlist(prep)
(Sys.time() - t.start0)
head(output)

## Save output
write.csv(output, file = "PBDBformatted.csv", row.names = FALSE)



## 4. Check for homonyms and possibly duplicate names ---------------------------

# Most genera with multiple entries are legitimate, caused by listing the genus
# as a whole, plus each subgenus separately. Saves the list to file specified
# below.

# Do you want to return the list of genera with subgenera? (DEFAULT = FALSE)
return.subgenera <- FALSE

mults <- sort(table(output$Genus), decreasing = TRUE)
mults <- mults[mults >= 2]
head(mults, 20)
file.name <- "multiGenera.txt"
sq <- 1:19     # Higher taxonomy columns
cat("The presence of subgenera, homonyms, and possible duplicates equals", 
    round(100 * length(mults) / nrow(output), 1), "% of the database\n", 
    file = file.name)
output[which(output$Genus == "Acanthopyge"), ] # Example of multiple subgenera

for(d in 1:length(mults)) {
  sus.gen <- names(mults[d])
  suspicious <- output[which(output$Genus == sus.gen), ]
  classes <- unique(suspicious$Class)
  if (length(classes) == 1L)
    # Identify likely subgenera:
    if (return.subgenera & all(sapply(sq, function(sq)
      length(unique(suspicious[[sq]])) == 1)) &
      suspicious$Subgenus[1] == "" & all(suspicious$Subgenus[-1] != ""))
      cat("OK: Genus", names(mults[d]), "has", nrow(suspicious) - 1, 
          "subgenera.\n", file = file.name, append = TRUE)
    # Identify likely problematic duplicated entries:
    if (any(sapply(sq, function(sq) 
      length(unique(suspicious[[sq]])) != 1)) & length(classes) == 1L)
      cat("WARNING: Genus", names(mults[d]), 
          "may be a duplicate genus entry. Investigate and override in PBDB if true.\n", 
          file = file.name, append = TRUE)
  # Identify likely legitimate homonyms:
  if (length(classes) == 2L)
    cat("OK: Genus", names(mults[d]), 
        "is a homonym for genera in difference classes:", classes, "\n", 
        file = file.name, append = TRUE)
}

# For any genera tagged as "WARNING", the best-practice is to add a new taxon to
# the PBDB that overrides the duplicate (and to re-classify their occurrences).
