test.uniprot.wsSearch.empty <- function(db) {
    n <- 2
    result <- db$wsSearch(fields='accession', format='tsv', size=n)
    expect_true( ! is.null(result))
    expect_true( ! is.na(result))
    expect_true(nchar(result) > 0)
    readtc <- textConnection(result, "r", local=TRUE)
    df <- read.table(readtc, sep="\t", header=TRUE)
    expect_true(colnames(df) == 'Entry')
    expect_true(nrow(df) == 0)
}

test.uniprot.wsSearch.by.name <- function(db) {
    n <- 2
    result <- db$wsSearch(query='protein_name:"prion protein"',
        fields='accession', format='tsv', size=n)
    expect_true( ! is.null(result))
    expect_true( ! is.na(result))
    expect_true(nchar(result) > 0)
    readtc <- textConnection(result, "r", local=TRUE)
    df <- read.table(readtc, sep="\t", header=TRUE)
    expect_true(colnames(df) == 'Entry')
    expect_true(nrow(df) == n)
}

test.uniprot.wsSearch.multiple.columns <- function(db) {
    n <- 2
    results <- db$wsSearch(query='e', fields=c('accession', 'id'), format='tsv',
        size=n, retfmt='parsed')
    testthat::expect_is(results, 'data.frame')
    testthat::expect_true(all(c('Entry', 'Entry Name') %in% colnames(results)))
    testthat::expect_true(nrow(results) == n)
}

# Set context
biodb::testContext("Test web services")

# Instantiate Biodb
biodb <- biodb::createBiodbTestInstance(ack=TRUE)

# Load package definitions
file <- system.file("definitions.yml", package='biodbUniprot')
biodb$loadDefinitions(file)

# Create connector
conn <- biodb$getFactory()$createConn('uniprot')

# Run tests
biodb::testThat('Uniprot entries query works fine with an empty query.',
                test.uniprot.wsSearch.empty, conn=conn)
biodb::testThat('Uniprot entries query works fine with multiple columns',
                test.uniprot.wsSearch.multiple.columns, conn=conn)
biodb::testThat('Uniprot entries query works fine with a query by name.',
                test.uniprot.wsSearch.by.name, conn=conn)

# Terminate Biodb
biodb$terminate()
