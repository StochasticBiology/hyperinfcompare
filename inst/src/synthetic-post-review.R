library(hyperinf)
library(hyperdags)
library(ggpubr)
library(ape)
library(dplyr)

sf = 2

run.1 = TRUE     # very different example
run.2 = TRUE     # very similar example
nsamp = 100

if(run.1 == TRUE) {

  # synthetic data
  set.seed(1)

  ########### very different pathways

  res = data.frame()
  for(L in c(3, 5, 7)) {
    # basic set of states from this process
    data.a <- 1*lower.tri(matrix(0, L, L), diag = TRUE)
    for(n in c(2, 5, 10, 20)) {
      for(sample in 1:10) {
        # sample from base set, its inverse, and base set again
        data.a.tmp = data.a[sample(1:L, n, replace=TRUE),]
        data.b.tmp = 1-data.a[sample(1:L, n, replace=TRUE),]
        data.c.tmp = data.a[sample(1:L, n, replace=TRUE),]
        # run HyperHMM inference on each sample
        fit.a.tmp = hyperinf(data.a.tmp, boot.parallel = 10)
        fit.b.tmp = hyperinf(data.b.tmp, boot.parallel = 10)
        fit.c.tmp = hyperinf(data.c.tmp, boot.parallel = 10)
        # loop through thresholds, using method to characterise relative and absolute ordering differences
        for(p in (1:7)/10) {
          for(q in (1:7)/10) {
            tmp.co = compare_orderings(fit.a.tmp, fit.b.tmp,
                                       threshold=p, percentile=q,
                                       type = "relative", n.samples = nsamp)
            res = rbind(res, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "diff", method = "relative", diffs = nrow(tmp.co)))
            tmp.co = compare_orderings(fit.a.tmp, fit.b.tmp,
                                       threshold=p, percentile=q,
                                       type = "absolute", n.samples = nsamp)
            res = rbind(res, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "diff", method = "absolute", diffs = nrow(tmp.co)))

            tmp.co = compare_orderings(fit.a.tmp, fit.c.tmp,
                                       threshold=p, percentile=q,
                                       type = "relative", n.samples = nsamp)
            res = rbind(res, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "same", method = "relative", diffs = nrow(tmp.co)))
            tmp.co = compare_orderings(fit.a.tmp, fit.c.tmp,
                                       threshold=p, percentile=q,
                                       type = "absolute", n.samples = nsamp)
            res = rbind(res, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "same", method = "absolute", diffs = nrow(tmp.co)))
          }
        }
      }
    }
  }
  save(res, file="res-synth.Rdata")
}

if(run.2 == TRUE) {
  ########### very similar pathways

  res.2 = data.frame()
  for(L in c(3, 5, 7)) {
    # two basic state sets, one a slight variant on the other
    data.a <- 1*lower.tri(matrix(0, L, L), diag = TRUE)
    data.b <- 1*lower.tri(matrix(0, L, L), diag = TRUE)
    data.b[1,1] = 0; data.b[1,2] = 1;
    for(n in c(2, 5, 10, 20)) {
      for(sample in 1:10) {
        # sample from both sets
        data.a.tmp = data.a[sample(1:L, n, replace=TRUE),]
        data.b.tmp = data.b[sample(1:L, n, replace=TRUE),]
        # use HyperHMM to perform inference
        fit.a.tmp = hyperinf(data.a.tmp, boot.parallel = 10)
        fit.b.tmp = hyperinf(data.b.tmp, boot.parallel = 10)
        # loop through thresholds, using method to characterise relative and absolute ordering differences
        for(p in (1:7)/10) {
          for(q in (1:7)/10) {
            tmp.co = compare_orderings(fit.a.tmp, fit.b.tmp,
                                       threshold=p, percentile=q,
                                       type = "relative", n.samples = nsamp)
            res.2 = rbind(res.2, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "diff", method = "relative", diffs = nrow(tmp.co)))
            tmp.co = compare_orderings(fit.a.tmp, fit.b.tmp,
                                       threshold=p, percentile=q,
                                       type = "absolute", n.samples = nsamp)
            res.2 = rbind(res.2, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "diff", method = "absolute", diffs = nrow(tmp.co)))
          }
        }
      }
    }
  }
  save(res.2, file="res-2-synth.Rdata")
}

#######
# analyse the outputs

load("res-synth.Rdata")
load("res-2-synth.Rdata")

# bindiffs will be binarised differences (was at least one difference found?)
res$bindiffs = 0
res$bindiffs[res$diffs > 0] = 1

# take means over samples for absolute and relative orderings
res.abs.mean = res[res$method=="absolute",] %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs), bindiffs = mean(bindiffs), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = c(bindiffs, diffs)
  )
res.rel.mean = res[res$method=="relative",] %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs), bindiffs = mean(bindiffs), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = c(bindiffs, diffs)
  )

# same treatment for small difference pair
res.2$bindiffs = 0
res.2$bindiffs[res.2$diffs > 0] = 1

# take means over samples for absolute and relative orderings
res.2.abs.mean = res.2[res.2$method == "absolute",] %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs), bindiffs = mean(bindiffs), .groups = "drop")
colnames(res.2.abs.mean)[6:7] = c("diffs_diff", "bindiffs_diff")
res.2.rel.mean = res.2[res.2$method == "relative",] %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs), bindiffs = mean(bindiffs), .groups = "drop")
colnames(res.2.rel.mean)[6:7] = c("diffs_diff", "bindiffs_diff")

# join the null comparison from the first batch to the second for FPR
res.2.abs.mean <- res.2.abs.mean |>
  dplyr::left_join(
    res.abs.mean |>
      dplyr::select(L, n, p, q, bindiffs_same, diffs_same),
    by = c("L", "n", "p", "q")
  )
res.2.rel.mean <- res.2.rel.mean |>
  dplyr::left_join(
    res.rel.mean |>
      dplyr::select(L, n, p, q, bindiffs_same, diffs_same),
    by = c("L", "n", "p", "q")
  )

##### binary differences (FWER)

# big set of ROC curves for absolute/relative and big/small differences
big.diag = ggarrange(
  ggplot(res.abs.mean,
         aes(x = bindiffs_same, y = bindiffs_diff, color = p, shape = factor(q))) +
    geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR", shape = "q", color="p") +
    facet_grid(factor(paste("L =",L),
                      levels=unique(paste("L =", L))) ~
                 factor(paste("n =", n),
                        levels=unique(paste("n =", n)))) +
    xlim(0,1) + ylim(0,1) + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust = 1))
  ,
  ggplot(res.2.abs.mean,
         aes(x = bindiffs_same, y = bindiffs_diff, color = p, shape = factor(q))) +
    geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR", shape = "q", color="p") +
    facet_grid(factor(paste("L =",L),
                      levels=unique(paste("L =", L))) ~
                 factor(paste("n =", n),
                        levels=unique(paste("n =", n)))) +
    xlim(0,1) + ylim(0,1)+ theme_minimal()+ theme(axis.text.x = element_text(angle=45, hjust = 1))
  ,
  ggplot(res.rel.mean,
         aes(x = bindiffs_same, y = bindiffs_diff, color = p, shape = factor(q))) +
    geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR", shape = "q", color="p") +
    facet_grid(factor(paste("L =",L),
                      levels=unique(paste("L =", L))) ~
                 factor(paste("n =", n),
                        levels=unique(paste("n =", n)))) +
    xlim(0,1) + ylim(0,1) + theme_minimal()+ theme(axis.text.x = element_text(angle=45, hjust = 1))
  ,
  ggplot(res.2.rel.mean,
         aes(x = bindiffs_same, y = bindiffs_diff, color = p, shape = factor(q))) +
    geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR", shape = "q", color="p") +
    facet_grid(factor(paste("L =",L),
                      levels=unique(paste("L =", L))) ~
                 factor(paste("n =", n),
                        levels=unique(paste("n =", n)))) +
    xlim(0,1) + ylim(0,1) + theme_minimal()+ theme(axis.text.x = element_text(angle=45, hjust = 1)),
  labels = c("Ai", "ii", "Bi", "ii"),
  common.legend = TRUE,
  legend = "right"
)

# FPR for different p, q
little.diag =
  ggplot(res.abs.mean, aes(x=p, y=q, fill=bindiffs_same)) + geom_tile() +
  facet_grid(factor(paste("L =",L),
                    levels=unique(paste("L =", L))) ~
               factor(paste("n =", n),
                      levels=unique(paste("n =", n)))) +
  labs(fill = "FPR") + theme_minimal()

sf = 2
png("all-diags-update.png", width=700*sf, height=1000*sf, res=72*sf)
ggarrange(big.diag, little.diag, nrow = 2, labels = c("", "C"), heights=c(2, 1))
dev.off()

###### main text figure

set.seed(1)
L = 5
n = 10

# as above, sample a set of states
data.a <- 1*lower.tri(matrix(0, L, L), diag = TRUE)
data.b = data.a
data.b[1,1] = 0; data.b[1,2] = 1;

data.a.tmp = data.a[sample(1:L, n, replace=TRUE),]
data.b.tmp = data.b[sample(1:L, n, replace=TRUE),]
data.c.tmp = 1-data.a[sample(1:L, n, replace=TRUE),]

# do the inference
fit.a.tmp = hyperinf(data.a.tmp, boot.parallel = 10)
fit.b.tmp = hyperinf(data.b.tmp, boot.parallel = 10)
fit.c.tmp = hyperinf(data.c.tmp, boot.parallel = 10)

# plot the ordering comparisons
comp.1.plot = plot_hyperinf_compare_orderings(fit.a.tmp, fit.c.tmp,expt.names = c("i", "ii"),
                                              thetastep=3, p.scale=0.3)
comp.2.plot = plot_hyperinf_compare_orderings(fit.a.tmp, fit.b.tmp,expt.names = c("i", "iii"),
                                              thetastep=3, p.scale=0.3)
recolor2 = scale_fill_manual(values=c("steelblue", "tomato", "white"))
nominorx = scale_x_continuous(breaks = scales::breaks_width(1), minor_breaks = NULL)
nominory = scale_y_continuous(breaks = scales::breaks_width(1), minor_breaks = NULL)

# combine data plot, inference plot, and ROC curves
mega.plot=    ggarrange(
  ggarrange(plot_hyperinf_data(data.a.tmp),
            plot_hyperinf_data(data.c.tmp),
            plot_hyperinf_data(data.b.tmp), nrow=1, labels=c("Ai", "ii", "iii")),
  ggarrange(comp.1.plot + recolor2 + nominorx + nominory + scale_x_continuous(breaks = 1:5) + labs(fill = "Model", y = "Feature"),
            comp.2.plot + recolor2 + nominorx + nominory + labs(fill = "Model", y = "Feature"), labels = c("B", "C")),
  ggarrange(
    ggplot(res.abs.mean[res.abs.mean$L == 5 & res.abs.mean$n == 10,],
           aes(x = bindiffs_same, y = bindiffs_diff, color = p, shape = factor(q))) +
      geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR", shape = "q", color="p") +
      xlim(0,1) + ylim(0,1) + theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust = 1))
    ,
    ggplot(res.2.abs.mean[res.2.abs.mean$L == 5 & res.2.abs.mean$n == 10,],
           aes(x = bindiffs_same, y = bindiffs_diff, color = p, shape = factor(q))) +
      geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR", shape = "q", color="p") +
      xlim(0,1) + ylim(0,1)+ theme_minimal()+ theme(axis.text.x = element_text(angle=45, hjust = 1))
    ,
    labels = c("D", "E"),
    common.legend = TRUE,
    legend = "bottom"
  ),
  nrow = 3,
  heights = c(1, 1, 1)
)

png("mega-plot-diags.png", width=600*sf, height=600*sf, res=72*sf)
print(mega.plot)
dev.off()
