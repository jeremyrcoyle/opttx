library(origami)
library(reshape2)
setwd("~/opttx/inst/simulations/")
files <- dir(path = "simresults", pattern = "nestingsim", full.names = T)
B <- length(files)
allresults <- lapply(files, function(file) {
  load(file)
  results
})
#load("simresults/nesting_sim_all.rdata")

results <- apply(do.call(rbind, allresults), 2, as.list)
combined <- combine_results(results)

alist <- combined$times[[1]]
listtodf <- function(alist) {
  df <- data.frame(do.call(rbind, alist))
  df$condition <- rownames(df)
  df
}

# run times
conditions <- c("splitrule", "fullrule", "nestrule", "tmle", "cvtmle", "nestcvtmle")
goodconditions <- c("SplitSequential Rule", "FullSequential Rule", "Nested Rule", "(FullSequential) TMLE", "SplitSequential CV-TMLE", 
                    "Nested CV-TMLE")

times <- ldply(combined$times, listtodf)
times$iteration <- rep(seq_len(B), each = nrow(times)/B)
times$goodcondition <- goodconditions[match(times$condition, conditions)]
times$goodcondition <- factor(times$goodcondition, goodconditions)


pdf("runtimes.pdf")
ggplot(times, aes(x = goodcondition, y = elapsed)) + geom_boxplot() + theme_bw() + xlab("Method") + ylab("Computation Time (s)") + 
  coord_flip()
dev.off()

# rule performance
combined$test_perf$iteration=1:nrow(combined$test_perf)
perfs<-melt(combined$test_perf,id=c("iteration"),variable_name="condition")
head(perfs)
perfs$condition=gsub("test_","",perfs$condition)
table(perfs$condition)

conditions <- c("EYd0", "split_EYdn", "full_EYdn", "nested_EYdn", "test_EYd0", "test_split_EYdn", "test_full_EYdn", "test_nested_EYdn")
goodconditions <- c("d0", "SplitSequential dn", "FullSequential dn", "Nested dn", "test set d0", "test set SplitSequential dn", 
                    "test set FullSequential dn", "test set Nested dn")
#perfs$iteration <- rep(seq_len(B), each = nrow(perfs)/B)
perfs$goodcondition <- goodconditions[match(perfs$condition, conditions)]

pdf("ruleperformance.pdf")
ggplot(perfs[!is.na(perfs$goodcondition),], aes(x = goodcondition, y = value)) + geom_boxplot() + theme_bw() + ylab("EYd") + xlab("Method") + coord_flip()
dev.off()

ggplot(perfs[perfs$condition%in%c("full_QaV_mse","nested_QaV_mse","split_QaV_mse"),], aes(x = condition, y = value)) + geom_boxplot() + theme_bw() + ylab("EYd") + xlab("Method") + coord_flip()

aggregate(value ~ condition, perfs, summary)

perfswide <- dcast(perfs, iteration ~ condition, value = "value")
mean(perfswide$split_EYdn - perfswide$full_EYdn)

mean(perfswide$split_EYdn - perfswide$nested_EYdn)


# estimated rule performance
ests <- combined$perf_ests
head(ests)
ests$iteration <- rep(seq_len(B), each = nrow(ests)/B)
aggregate(est ~ rule+estimator+rule_fit, ests, summary)
ests<-ests[ests$rule=="A=d(V)",]
# compare estimates to true EYdn
m <- merge(perfswide, ests)
ggplot(m, aes(x = estimator, y = est - split_EYdn)) + geom_boxplot() +coord_flip()+ theme_bw() + geom_hline(yintercept = 0, linetype = "dashed")+facet_wrap(~rule_fit)
m$rule_fit
estperf <- ddply(m, .(rule_fit,estimator), function(estdata) {
  if (estdata$rule_fit[1]=="split") {
    truth <- estdata$split_EYdn
  } else if (estdata$rule_fit[1]=="nested") {
    truth <- estdata$nested_EYdn
  } else {
    truth <- estdata$full_EYdn
  }
  #truth <- estdata$EYd0  #compare to true opt
  evar <- var(estdata$est)
  ebias <- mean(estdata$est - truth)
  emse <- mean((estdata$est - truth)^2)
  ecoverage <- mean((estdata$lower < truth) & (truth < estdata$upper))
  data.frame(var = evar, bias = ebias, mse = emse, coverage = ecoverage)
})

estperflong <- melt(estperf, id = c("rule_fit","estimator"))
estperflong$condition=sprintf("%s %s",estperflong$estimator,estperflong$rule_fit)
ggplot(estperflong,aes(y=condition,x=value))+geom_point()+facet_wrap(~variable,scales="free")+theme_bw()
estperflong$iscvest <- estperflong$condition %in% c("cvtmle", "cvtmlenest", "cvtmlenaive", "nestcvtmle")
conditions <- c("cvtmle_full", "cvtmle_nested", "cvtmle_split", "nestedcvtmle_full", "nestedcvtmle_split", "tmle_full", 
                "tmle_nested", "tmle_split")
goodconditions <- c("SplitSequential CV-TMLE FullSequential dn", "SplitSequential CV-TMLE Nested dn", "SplitSequential CV-TMLE SplitSequential dn", 
                    "Nested CV-TMLE FullSequential dn", "Nested CV-TMLE SplitSequential dn", "(FullSequential) TMLE FullSequential dn", 
                    "(FullSequential) TMLE Nested dn", "(FullSequential) TMLE SplitSequential dn")

ssconditions <- c("SplitSequential CV-TMLE SplitSequential dn", "Nested CV-TMLE SplitSequential dn", "(FullSequential) TMLE SplitSequential dn")


perft <- c("var", "bias", "mse", "coverage")
goodperft <- c("Variance", "Bias", "MSE", "Coverage")
estperflong$goodvariable <- goodperft[match(estperflong$variable, perft)]
estperflong$goodvariable <- factor(estperflong$goodvariable, levels = goodperft)
estperflong$goodcondition <- goodconditions[match(estperflong$condition, conditions)]
pdf("estperf.pdf")
ggplot(estperflong[estperflong$goodvariable != "Coverage" & estperflong$goodcondition %in% ssconditions, ], aes(y = goodcondition, 
                                                                                                                x = value)) + geom_point() + geom_vline(xintercept = 0, alpha = 0) + facet_wrap(~goodvariable, scales = "free", ncol = 1) + 
  theme_bw() + ylab("Estimator") + xlab("")
dev.off()

coverage <- estperflong[estperflong$goodvariable == "Coverage", ]
pdf("coverage.pdf")
ggplot(coverage, aes(y = condition, x = value)) + geom_point() + geom_vline(xintercept = 0.95, 
                                                                                                                            linetype = "dashed") + theme_bw() + ylab("Estimator") + xlab("")
dev.off()

coverage$value[coverage$condition=="CV-TMLE split"]
coverage$value[coverage$condition=="Nested CV-TMLE split"]
mean(perfs$value[perfs$condition=="split_EYdn"])
mean(ests$est[ests$rule_fit=="split"&ests$estimator=="CV-TMLE"])
mean(ests$est[ests$rule_fit=="split"&ests$estimator=="Nested CV-TMLE"])
TRUE