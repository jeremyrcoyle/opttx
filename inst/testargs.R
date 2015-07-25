Wnodes <- grep("^W", names(data), value = TRUE)
Anode <- "A"
Ynode <- "Y"
Vnodes <- Wnodes
stratifyAY <- TRUE


verbose <- 4
parallel <- FALSE
perf_tmle <- TRUE
perf_dripcw <- TRUE
perf_cv <- TRUE
perf_full <- TRUE
maximize <- TRUE
SL.library <- opt_tmle.SL.library 
