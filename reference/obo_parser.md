# Parse OBO file to list of data.tables

This is a generalised helper function to read in .obo files (ontology
files) and generate some useful format out of it.

## Usage

``` r
obo_parser(file_path, verbose = TRUE)
```

## Arguments

- file_path:

  String. Path to the .obo file

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list the following elements:

- names - data.table with name and ID of the term.

- alt_ids - data.table with alternative ids.

- namespace - data.table with the found name-spaces.

- term_subsets - data.table with the found subsets.

- term_obsolete - character vector of terms that were marked obsolete.

- term_relations - data.table with the ontological relationships found
  in the file.

- synonyms - data.table with synonyms provided in the .obo file.

- externalRef - data.table with external references
