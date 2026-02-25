# plots ------------------------------------------------------------------------

obj <- myc_gwn

rm(obj)

hist(myc_gwn@stats$similarity, breaks = 30)

hist(myc_gwn@permuted_embd[[1]], breaks = 30)

hist(myc_gwn@permuted_embd[[2]], breaks = 30)

hist(myc_gwn@permuted_embd[[3]], breaks = 30)

summary(myc_gwn@stats$similarity)

summary(myc_gwn@permuted_embd[[3]])


object <- myc_gwn


random_embeddings <- purrr::imap(object@permuted_embd, \(x, i) {
  data.table::data.table(similarity = x, sample = sprintf("perm_%i", i))
}) %>%
  append(
    .,
    list(data.table::data.table(
      similarity = object@stats$similarity,
      sample = "actual"
    ))
  ) %>%
  data.table::rbindlist()

ggplot(data = random_embeddings, mapping = aes(x = similarity)) +
  facet_wrap(~sample, scales = "free") +
  geom_histogram(aes(y = after_stat(density)), bins = 50) +
  xlab("Cosine similarity") +
  ylab("Density") +
  theme_bw()
