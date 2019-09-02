## Code to quickly download and summarize stratigraphic intervals used in the
## Paleobiology Database (PBDB). Written for Phil Novack-Gottshall, Benedictine
## University, Lisle, IL <pnovack-gottshall@ben.edu>

# Download all PBDB stratigraphic intervals (requires internet access)
strat_names <-
  read.csv("https://www.paleobiodb.org/data1.2/intervals/list.csv?all_records&vocab=pbdb")
# Specify desired level: 1 = eons, 2 = eras, 3 = periods, 4 (the default) =
# subperiods, and 5 = epochs.
scale_level <- 4
ages <- strat_names[which(strat_names$scale_level == scale_level),]
# If wish to add in the Ediacaran, for example:
edia <- strat_names[which(strat_names$interval_name == "Ediacaran"), ]
ages <- rbind(ages, edia)
# Calculate and add midpoints:
mid_ma <- apply(ages[ ,9:10], 1, mean)
ages <- cbind(ages, mid_ma = mid_ma)
tail(ages)
