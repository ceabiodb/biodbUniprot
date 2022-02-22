# biodbUniprot

[![Codecov test coverage](https://codecov.io/gh/pkrog/biodbUniprot/branch/master/graph/badge.svg)](https://codecov.io/gh/pkrog/biodbUniprot?branch=master)

An R package for accessing [Uniprot](https://www.uniprot.org/) online database,
based on R package/framework [biodb](https://github.com/pkrog/biodb/).

## Introduction

*biodbUniprot* is an extension package of the *biodb* package.
It allows to connect to UniProt for retrieving entries, searching for entries by
name or organism, and convert gene symbols to UniProt IDs.

## Examples


Getting entries:
```r
bdb <- boidb::newInst()
uniprot <- bdb$getFactory()$createConn('uniprot')
entries <- uniprot$getEntry(c('P01011', 'P09237'))
bdb$entriesToDataframe(entries)
```

Run a web service query:
```r
bdb <- boidb::newInst()
uniprot <- bdb$getFactory()$createConn('uniprot')
uniprot$wsQuery('reviewed:yes AND organism:9606', columns=c('id', 'entry name'),
    limit=2, retfmt='parsed')
```

## Installation

Install the latest stable version using Bioconductor:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('biodbUniprot')
```

## Documentation

See the introduction vignette:
```r
vignette('biodbUniprot', package='biodbUniprot')
```

## Citations

 * The UniProt Consortium. UniProt: the universal protein knowledgebase. Nucleic Acids Res. 45: D158-D169 (2017), <https://doi.org/10.1093/nar/gkw1099>.
 * Pundir S., Martin M.J., O’Donovan C. (2017) UniProt Protein Knowledgebase. In: Wu C., Arighi C., Ross K. (eds) Protein Bioinformatics. Methods in Molecular Biology, vol 1558. Humana Press, New York, NY. <https://doi.org/10.1007/978-1-4939-6783-4_2>.
