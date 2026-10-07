# node2vec ---------------------------------------------------------------------

spec_node2vec <- param_spec(
  name = "node2vec",
  title = "Wrapper function for the node2vec parameters",
  description = paste(
    "These are general node2vec parameters, designed for speed. They will",
    "default to `available cores - 2 (min 1)`."
  ),
  checker = "Node2Vec",
  label = "node2vec parameters",
  return_order = c(
    "p",
    "q",
    "walks_per_node",
    "walk_length",
    "num_workers",
    "batch_size",
    "n_epochs",
    "n_negatives",
    "window_size",
    "lr"
  ),
  fields = list(
    p = p_dbl(1.0, doc = "Return parameter for biased random walks."),
    q = p_dbl(1.0, doc = "In-out parameter for biased random walks."),
    walks_per_node = p_int(40L, doc = "Number of random walks per node."),
    walk_length = p_int(40L, doc = "Length of each random walk."),
    batch_size = p_int(256L, doc = "Batch size for training."),
    n_epochs = p_int(20L, doc = "Number of training epochs."),
    n_negatives = p_int(5L, doc = "Number of negative samples."),
    window_size = p_int(2L, doc = "Context window size."),
    lr = p_dbl(1e-2, doc = "Learning rate."),
    num_workers = p_int(
      NULL,
      null_ok = TRUE,
      check_as = "I1",
      doc = paste(
        "Number of worker threads. If kept to `NULL`, it resolves to",
        "`available cores - 2 (min 1)`."
      )
    )
  ),
  extra_ctor = quote(
    num_workers <- if (is.null(num_workers)) {
      max(1L, parallel::detectCores() - 2L)
    } else {
      num_workers
    }
  )
)

spec_genewalk <- param_spec(
  name = "genewalk",
  title = "Wrapper function for the GeneWalk-specific node2vec parameters",
  checker = "Node2Vec",
  label = "node2vec parameters",
  return_order = c(
    "p",
    "q",
    "walks_per_node",
    "walk_length",
    "num_workers",
    "batch_size",
    "n_epochs",
    "n_negatives",
    "window_size",
    "lr"
  ),
  fields = list(
    p = p_dbl(1.0, doc = "Return parameter for biased random walks."),
    q = p_dbl(1.0, doc = "In-out parameter for biased random walks."),
    walks_per_node = p_int(100L, doc = "Number of random walks per node."),
    walk_length = p_int(10L, doc = "Length of each random walk."),
    batch_size = p_int(256L, doc = "Batch size for training."),
    n_epochs = p_int(5L, doc = "Number of training epochs."),
    n_negatives = p_int(5L, doc = "Number of negative samples."),
    window_size = p_int(
      1L,
      doc = paste(
        "Context window size. For GeneWalk you want embeddings from direct",
        "neighbours."
      )
    ),
    lr = p_dbl(1e-2, doc = "Learning rate."),
    num_workers = p_int(
      1L,
      doc = paste(
        "Number of worker threads. Anything above `1L` makes the fitting much",
        "faster, at the cost of reproducibility."
      )
    )
  )
)

# metapath2vec -----------------------------------------------------------------

spec_metapath2vec <- param_spec(
  name = "metapath2vec",
  title = "Wrapper function for the metapath2vec parameters",
  description = paste(
    "Parameters for [metapath2vec()]. There is no `p` or `q`: metapath walks",
    "are first-order and follow the schema, not a biased return/explore rule."
  ),
  checker = "Metapath2Vec",
  label = "metapath2vec parameters",
  return_order = c(
    "walks_per_node",
    "walk_length",
    "num_workers",
    "n_epochs",
    "n_negatives",
    "window_size",
    "lr",
    "sample"
  ),
  fields = list(
    walks_per_node = p_int(
      40L,
      range = "[1,)",
      doc = paste(
        "Number of random walks per node of the metapath's starting type."
      )
    ),
    walk_length = p_int(
      40L,
      range = "[2,)",
      doc = paste(
        "Length of each random walk. Rounded up so the walk closes on a full",
        "number of schema cycles."
      )
    ),
    n_epochs = p_int(20L, range = "[1,)", doc = "Number of training epochs."),
    n_negatives = p_int(
      5L,
      range = "[1,)",
      doc = "Number of negative samples."
    ),
    window_size = p_int(2L, range = "[1,)", doc = "Context window size."),
    lr = p_dbl(1e-2, range = "(0,)", doc = "Learning rate."),
    sample = p_dbl(
      1e-3,
      range = "[0,)",
      doc = "Subsampling threshold for frequent nodes."
    ),
    num_workers = p_int(
      NULL,
      range = "[1,)",
      null_ok = TRUE,
      check_as = "I1[1,)",
      doc = paste(
        "Number of worker threads. If kept to `NULL`, it resolves to",
        "`available cores - 2 (min 1)`."
      )
    )
  ),
  extra_ctor = quote(
    num_workers <- if (is.null(num_workers)) {
      max(1L, parallel::detectCores() - 2L)
    } else {
      num_workers
    }
  )
)

# synthetic data ---------------------------------------------------------------

spec_genewalk_data <- param_spec(
  name = "genewalk_data",
  title = "GeneWalk synthetic data params",
  checker = "GeneWalkData",
  label = "GeneWalk data parameters",
  fields = list(
    n_signal_genes = p_int(150L, doc = "Number of signal genes."),
    n_noise_genes = p_int(150L, doc = "Number of noise genes."),
    n_roots = p_int(10L, doc = "Number of ontology root nodes."),
    depth = p_int(3L, doc = "Depth of each ontology subtree."),
    branching = p_int(3L, doc = "Base branching factor per ontology node."),
    p_lateral = p_dbl(
      0.1,
      doc = "Probability of a lateral edge between sibling nodes."
    ),
    p_ppi = p_dbl(
      0.05,
      doc = paste(
        "Probability of a PPI edge between any two genes in the same group."
      )
    ),
    min_annotations = p_int(
      10L,
      doc = "Minimum ontology term annotations per gene."
    ),
    max_annotations = p_int(
      25L,
      doc = "Maximum ontology term annotations per gene."
    ),
    min_noise_subtrees = p_int(
      15L,
      doc = "Minimum distinct ontology subtrees a noise gene must span."
    )
  )
)

# diffusion kernels ------------------------------------------------------------

spec_kernel <- param_spec(
  name = "kernel",
  title = "Wrapper function for the diffusion kernel parameters",
  description = paste(
    "These are parameters for the graph diffusion kernels. Not all parameters",
    "are used by all kernels, please see Details."
  ),
  details = paste(
    "Parameter usage by kernel type:",
    "\\itemize{",
    "  \\item \\code{regularised_laplacian}: \\code{sigma2}, \\code{add_diag}",
    "  \\item \\code{commute_time}: none",
    "  \\item \\code{diffusion}: \\code{sigma2}",
    "  \\item \\code{inverse_cosine}: none",
    "  \\item \\code{pstep}: \\code{a}, \\code{p}",
    "}",
    sep = "\n"
  ),
  checker = "Kernel",
  label = "kernel parameters",
  fields = list(
    sigma2 = p_dbl(
      1.0,
      doc = paste(
        "Bandwidth parameter controlling the extent of diffusion. Used by",
        "`regularised_laplacian` and `diffusion` kernels."
      )
    ),
    add_diag = p_dbl(
      1.0,
      doc = paste(
        "Regularisation term added to the diagonal. Used by",
        "`regularised_laplacian`. Must be positive for every connected",
        "component to be invertible."
      )
    ),
    a = p_dbl(
      3.0,
      doc = paste(
        "Regularisation parameter for the `pstep` kernel. Must be `>= 2`,",
        "since the eigenvalues of the normalised Laplacian sit in `[0, 2]`."
      )
    ),
    p = p_int(
      5L,
      doc = "Number of steps for the `pstep` kernel. Must be positive."
    )
  ),
  extra_ctor = quote({
    checkmate::assertTRUE(a >= 2, .var.name = "a >= 2")
    checkmate::assertTRUE(p > 0L, .var.name = "p > 0")
  }),
  extra_check = quote({
    if (x$a < 2) {
      return("Parameter `a` must be >= 2 for the pstep kernel.")
    }
    if (x$p < 1L) {
      return("Parameter `p` must be a positive integer.")
    }
  })
)

# diffusion profiles -----------------------------------------------------------

spec_diffusion_profiles <- param_spec(
  name = "diffusion_profiles",
  title = "Wrapper function for the diffusion profile parameters",
  description = paste(
    "Parameters for [generate_profiles()], i.e. the constrained personalised",
    "PageRank behind the diffusion profiles of Ruiz et al."
  ),
  checker = "DiffusionProfiles",
  label = "diffusion profile parameters",
  fields = list(
    alpha = p_dbl(
      0.85,
      range = "(0,1)",
      doc = "Probability of continuing the walk instead of restarting at the seed."
    ),
    max_iter = p_int(
      100L,
      range = "[1,)",
      doc = "Maximum number of power iterations."
    ),
    tol = p_dbl(
      1e-6,
      range = "(0,)",
      doc = "Convergence threshold on the L1 change between two iterations."
    )
  )
)
