library(hyperinf)
library(hyperdags)
library(ggpubr)
library(ape)
library(dplyr)

sf = 2

run.1 = TRUE     # very different example
run.2 = TRUE     # very similar example

if(run.1 == TRUE) {
  
  # synthetic data
  set.seed(1)
  
  ########### very different pathways
  
  res = data.frame()
  for(L in c(3, 5, 7)) {
    data.a <- 1*lower.tri(matrix(0, L, L), diag = TRUE)
    for(n in c(2, 4, 20)) {
      for(sample in 1:20) {
        data.a.tmp = data.a[sample(1:L, n, replace=TRUE),]
        data.b.tmp = 1-data.a[sample(1:L, n, replace=TRUE),]
        data.c.tmp = data.a[sample(1:L, n, replace=TRUE),]
        fit.a.tmp = hyperinf(data.a.tmp, boot.parallel = 30)
        fit.b.tmp = hyperinf(data.b.tmp, boot.parallel = 30)
        fit.c.tmp = hyperinf(data.c.tmp, boot.parallel = 30)
        for(p in (1:7)/10) {
          for(q in (1:7)/10) {
            tmp.co = compare_orderings(fit.a.tmp, fit.b.tmp, threshold=p, percentile=q)
            res = rbind(res, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "diff", diffs = nrow(tmp.co)))
            tmp.co = compare_orderings(fit.a.tmp, fit.c.tmp, threshold=p, percentile=q)
            res = rbind(res, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "same", diffs = nrow(tmp.co)))
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
    data.b <- 1*lower.tri(matrix(0, L, L), diag = TRUE)
    data.b[1,1] = 0; data.b[1,2] = 1;
    for(n in c(2, 4, 20)) {
      for(sample in 1:20) {
        data.a.tmp = data.b[sample(1:L, n, replace=TRUE),]
        data.b.tmp = 1-data.b[sample(1:L, n, replace=TRUE),]
        fit.a.tmp = hyperinf(data.a.tmp, boot.parallel = 30)
        fit.b.tmp = hyperinf(data.b.tmp, boot.parallel = 30)
        for(p in (1:7)/10) {
          for(q in (1:7)/10) {
            tmp.co = compare_orderings(fit.a.tmp, fit.b.tmp, threshold=p, percentile=q)
            res.2 = rbind(res.2, data.frame(L = L, n = n, sample=sample, p = p, q = q, type = "diff", diffs = nrow(tmp.co)))
          }
        }
      }
    }
  }
  save(res.2, file="res-2-synth.Rdata")
}

#######

load("res-synth.Rdata")
load("res-2-synth.Rdata")

ggarrange(
  ggplot(res[res$L == 3,], aes(x = factor(p), y = factor(q), fill = diffs)) + geom_tile() + facet_grid(n~type),
  ggplot(res[res$L == 5,], aes(x = factor(p), y = factor(q), fill = diffs)) + geom_tile() + facet_grid(n~type),
  ggplot(res[res$L == 10,], aes(x = factor(p), y = factor(q), fill = diffs)) + geom_tile() + facet_grid(n~type)
)

binary.class = TRUE
if(binary.class == TRUE) {
res$diffs[res$diffs > 0] = 1
res.mean = res %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = diffs
  )
res.2$diffs[res.2$diffs > 0] = 1
res.2.mean = res.2 %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = diffs
  )
} else {
res.mean = res %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs/(L*L)), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = diffs
  )
res.2.mean = res.2 %>%
  dplyr::group_by(L, n, p, q, type) %>%
  dplyr::summarise(diffs = mean(diffs/(L*L)), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = diffs
  )
}

res.2.mean <- res.2.mean |>
  dplyr::left_join(
    res.mean |>
      dplyr::select(L, n, p, q, same),
    by = c("L", "n", "p", "q")
  )

ggarrange(
ggplot(res.mean[res.mean$n <= 20,], 
          aes(x = same, y = diff, color = p, shape = factor(q))) + 
  geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR") +
  facet_grid(L ~ factor(n))
,
ggplot(res.2.mean[res.2.mean$n <= 20,], 
       aes(x = same, y = diff, color = p, shape = factor(q))) + 
  geom_line() + geom_point() + geom_abline() + labs(x = "FPR", y = "TPR") +
  facet_grid(L ~ factor(n))
)

ggarrange(
  ggplot(res.mean[res.mean$L == 3,], aes(x = factor(p), y = factor(q), fill = diffs)) + geom_tile() + facet_grid(n~type),
  ggplot(res.mean[res.mean$L == 5,], aes(x = factor(p), y = factor(q), fill = diffs)) + geom_tile() + facet_grid(n~type),
  ggplot(res.mean[res.mean$L == 10,], aes(x = factor(p), y = factor(q), fill = diffs)) + geom_tile() + facet_grid(n~type)
)

ggplot(res.mean[res.mean$L == 3,], aes(x = factor(p), y = factor(q), fill = diffs)) + geom_tile() + facet_grid(n~type),

#########

