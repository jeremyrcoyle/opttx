Wnodes <- grep("^W", names(data), value = T)
Anode <- "A"
Ynode <- "Y"
Vnodes <- Wnodes
stratifyAY <- TRUE


verbose <- 4
parallel <- F
perf_tmle <- T
perf_dripcw <- T
perf_cv <- T
perf_full <- T 
