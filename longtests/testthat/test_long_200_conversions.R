test_geneSymbolsToUniprotIds_null <- function(conn) {
    
    # Null list
    testthat::expect_equal(conn$geneSymbolToUniprotIds(NULL), list())
    testthat::expect_equal(conn$geneSymbolToUniprotIds(character()), list())
    testthat::expect_equal(conn$geneSymbolToUniprotIds(NA_character_), list())
}

test_geneSymbolsToUniprotIds_exact <- function(conn) {

    # Exact search
    expected_ids <- list('TGF-b1'='Q4LDM3')
    not_expected_ids <- list('TGF-b1'=c('W1I9X7', 'Q08FI9'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids))
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids)) {
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
    }
}

test_geneSymbolsToUniprotIds_exact_2_genes <- function(conn) {

    # Exact search with two genes
    expected_ids <- list('TGF-b1'='Q4LDM3', 'G-CSF'='Q9GJU0')
    not_expected_ids <- list('TGF-b1'=c('W1I9X7', 'Q08FI9'),
                             'G-CSF'=c('P09919', 'Q8N4W3', 'A0A3G2Y4F6'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids))
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids)) {
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
    }
}

test_geneSymbolsToUniprotIds_ignore_alphanum <- function(conn) {

    # Ignore non-alphanum chars (e.g.: "G-CSF" == "GCSF")
    expected_ids <- list('TGF-b1'='P01137', 'G-CSF'=c('Q9GJU0', 'P09919'))
    not_expected_ids <- list('TGF-b1'=c('W1I9X7', 'Q08FI9'),
        'G-CSF'=c('Q8N4W3', 'A0A3G2Y4F6'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids),
        ignore.nonalphanum=TRUE)
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids)) {
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
    }
}
    
# XXX Partial match does not work anymore in UniProt
# XXX New UniProt REST API returns only exact match.
#test_geneSymbolsToUniprotIds_partial_match <- function(conn) {
#
#    # Partial match (e.g.: "G-CSFb1" matches when searching for "G-CSF")
#    expected_ids <- list('TGF-b1'=c('Q4LDM3', 'W1I9X7'),
#        'G-CSF'=c('Q9GJU0', 'A0A3G2Y4F6', 'C0STS3'))
#    not_expected_ids <- list('TGF-b1'=c('Q08FI9'),
#        'G-CSF'=c('P09919', 'Q8N4W3'))
#    ids <- conn$geneSymbolToUniprotIds(names(expected_ids), partial.match=TRUE,
#        ignore.nonalphanum=TRUE)
#    testthat::expect_equal(names(ids), names(expected_ids))
#    for (gene in names(ids)) {
#        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
#        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
#    }
#}

# XXX Filtering option is deprecated.
# XXX New UniProt REST API returns only exact match.
#test_geneSymbolsToUniprotIds_no_filtering <- function(conn) {
#
#    # No filtering
#    expected_ids <- list('TGF-b1'=c('Q4LDM3', 'W1I9X7'),
#        'G-CSF'=c('C0STS3', 'A0A679AQ73', 'Q8MKE0'))
#    ids <- conn$geneSymbolToUniprotIds(names(expected_ids), filtering=FALSE)
#    testthat::expect_equal(names(ids), names(expected_ids))
#    for (gene in names(ids))
#        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
#}

# Set context
biodb::testContext("Test conversions")

# Instantiate Biodb
biodb <- biodb::createBiodbTestInstance(ack=TRUE)

# Load package definitions
file <- system.file("definitions.yml", package='biodbUniprot')
biodb$loadDefinitions(file)

# Create connector
conn <- biodb$getFactory()$createConn('uniprot')

# Run tests
biodb::testThat('We can convert gene symbols to UniProt IDs.',
                test_geneSymbolsToUniprotIds_null, conn=conn)
biodb::testThat('We can convert gene symbols to UniProt IDs.',
                test_geneSymbolsToUniprotIds_exact, conn=conn)
biodb::testThat('We can convert gene symbols to UniProt IDs.',
                test_geneSymbolsToUniprotIds_exact_2_genes, conn=conn)
biodb::testThat('We can convert gene symbols to UniProt IDs.',
                test_geneSymbolsToUniprotIds_ignore_alphanum, conn=conn)
#biodb::testThat('We can convert gene symbols to UniProt IDs.',
#                test_geneSymbolsToUniprotIds_partial_match, conn=conn)
#biodb::testThat('We can convert gene symbols to UniProt IDs.',
#                test_geneSymbolsToUniprotIds_no_filtering, conn=conn)

# Terminate Biodb
biodb$terminate()
