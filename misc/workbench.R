# plots ------------------------------------------------------------------------

obj <- myc_gwn

rm(obj)

hist(myc_gwn@stats$similarity, breaks = 30)

hist(myc_gwn@permuted_embd[[1]], breaks = 30)

hist(myc_gwn@permuted_embd[[2]], breaks = 30)

hist(myc_gwn@permuted_embd[[3]], breaks = 30)

summary(myc_gwn@stats$similarity)

summary(myc_gwn@permuted_embd[[3]])
