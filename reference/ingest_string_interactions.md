# Pull out the string data from OpenTargets

Pull out the string data from OpenTargets

## Usage

``` r
ingest_string_interactions(dir_data, threshold = 0.85)
```

## Arguments

- dir_data:

  String. The directory in which you downloaded the data

- threshold:

  Numeric. The string threshold from which to consider a STRING
  connection high quality.

## Value

The data.table with the STRING network data.
