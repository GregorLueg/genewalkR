# Gene Walk Network Generator

R6 class to build Gene Walk networks from internal database sources.
Supports building a master network once and subsetting to multiple gene
sets.

## Public fields

- `pathway_sources`:

  Character vector of pathway sources

- `pathway_namespaces`:

  List of namespace filters per pathway source

- `ppi_sources`:

  Character vector of PPI sources

- `ppi_params`:

  List of PPI filtering parameters

## Methods

### Public methods

- [`GeneWalkGenerator$new()`](#method-GeneWalkGenerator-new)

- [`GeneWalkGenerator$print()`](#method-GeneWalkGenerator-print)

- [`GeneWalkGenerator$add_pathways()`](#method-GeneWalkGenerator-add_pathways)

- [`GeneWalkGenerator$add_ppi()`](#method-GeneWalkGenerator-add_ppi)

- [`GeneWalkGenerator$build()`](#method-GeneWalkGenerator-build)

- [`GeneWalkGenerator$create_for_genes()`](#method-GeneWalkGenerator-create_for_genes)

- [`GeneWalkGenerator$reset_choices()`](#method-GeneWalkGenerator-reset_choices)

- [`GeneWalkGenerator$return_full_network_dt()`](#method-GeneWalkGenerator-return_full_network_dt)

- [`GeneWalkGenerator$clone()`](#method-GeneWalkGenerator-clone)

------------------------------------------------------------------------

### Method `new()`

Initialise the generator

#### Usage

    GeneWalkGenerator$new()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print for GeneWalkGenerator

#### Usage

    GeneWalkGenerator$print(...)

#### Arguments

- `...`:

  Additional parameters to forward to the print

------------------------------------------------------------------------

### Method `add_pathways()`

Add pathway sources to the network

#### Usage

    GeneWalkGenerator$add_pathways(
      source = c("go", "reactome"),
      go_namespace = c("biological_process", "molecular_function", "cellular_component")
    )

#### Arguments

- `source`:

  Character vector of sources: `"go"`, `"reactome"`. Defaults to `"go"`
  only by default.

- `go_namespace`:

  Character vector. Specifically for GO. Choices are
  `c("biological_process", "molecular_function", "cellular_component")`

------------------------------------------------------------------------

### Method `add_ppi()`

Add PPI sources to the network

#### Usage

    GeneWalkGenerator$add_ppi(
      source = c("combined", "string", "signor", "reactome", "intact"),
      string_threshold = NULL,
      intact_threshold = NULL,
      intact_physical_only = FALSE
    )

#### Arguments

- `source`:

  Character vector of sources

- `string_threshold`:

  Numeric threshold for STRING scores

- `intact_threshold`:

  Numeric threshold for Intact scores

- `intact_physical_only`:

  Logical, only physical interactions from Intact

------------------------------------------------------------------------

### Method `build()`

Build the full network from selected sources

#### Usage

    GeneWalkGenerator$build(.verbose = TRUE)

#### Arguments

- `.verbose`:

  Boolean. Controls the verbosity of the function

------------------------------------------------------------------------

### Method `create_for_genes()`

Create a gene-specific GeneWalk object

#### Usage

    GeneWalkGenerator$create_for_genes(genes)

#### Arguments

- `genes`:

  Character vector of gene symbols

#### Returns

Returns the initialised `GeneWalk`.

------------------------------------------------------------------------

### Method `reset_choices()`

Resets the internal choices and erases any stored network data.

#### Usage

    GeneWalkGenerator$reset_choices()

------------------------------------------------------------------------

### Method `return_full_network_dt()`

Returns the full network data.table.

#### Usage

    GeneWalkGenerator$return_full_network_dt()

#### Returns

data.table with the full internal network.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    GeneWalkGenerator$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
