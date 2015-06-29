Wnodes=grep("^W", names(data), value = T)
Anode="A"
Ynode="Y"
Vnodes="W3"
stratifyAY=TRUE
Q_library=c("SL.glm", "SL.glmem", "SL.glmnet", "SL.step.forward", "SL.gam", "SL.rpart", "SL.rpartPrune", 
            "SL.mean")
g_library=c("SL.glm", "SL.glmnet", "SL.step.forward", "SL.gam", "SL.rpart", "SL.rpartPrune", "SL.mean")
blip_library=c("SL.glm", "SL.glmnet", "SL.step.forward", "SL.gam", 
               "SL.rpart", "SL.rpartPrune", "SL.mean")
maximize=T
verbose=2
parallel=F