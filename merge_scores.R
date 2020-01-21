#!R

#################################################
# Merge all Polygenic Risk Scores into one file #
#################################################

library(magrittr)
library(data.table)

get_folder <- function(folder) {
	data <- list.files(path=folder, pattern="profile", full.names=T) %>%
        lapply(fread) %>%
        rbindlist()
	scores <- data[, .(score = sum (SCORESUM)), by=IID]
	scores[, model := folder]
	scores
}

scores <- list.dirs(path="scores", full.names=TRUE, recursive=FALSE) %>%
    lapply(get_folder) %>%
    rbindlist()

# The "model" column contains the names of the folders in the scores directory.
# These are not nice to work with.
# Let's make renaming easier
scores[, model := as.factor(model)]

# Now you can name the models something that makes sense
# E.g.
levels(scores$model) <- c("p1e1", "p1e2", "p1e3", "p1e0", "p3e1", "p3e2", "p3e3", "inf")

# From here on, you can merge your phenotype data into the data.table

# Sometimes, it might be more convenient to have each model in a separate column:
wide <- dcast(scores, ... ~ model, value.var="score")