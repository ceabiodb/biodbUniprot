test.uniprot.wsQuery.empty <- function(db) {
	n <- 2
	result <- db$wsQuery(columns = 'id', format = 'tab', limit = n)
	expect_true( ! is.null(result))
	expect_true( ! is.na(result))
	expect_true(nchar(result) > 0)
	readtc <- textConnection(result, "r", local = TRUE)
	df <- read.table(readtc, sep = "\t", header = TRUE)
	expect_true(colnames(df) == 'Entry')
	expect_true(nrow(df) == n)
}

test.uniprot.wsQuery.by.name <- function(db) {
	n <- 2
	result <- db$wsQuery(query = 'name:"prion protein"', columns = 'id', format = 'tab', limit = n)
	expect_true( ! is.null(result))
	expect_true( ! is.na(result))
	expect_true(nchar(result) > 0)
	readtc <- textConnection(result, "r", local = TRUE)
	df <- read.table(readtc, sep = "\t", header = TRUE)
	expect_true(colnames(df) == 'Entry')
	expect_true(nrow(df) == n)
}

test.uniprot.wsQuery.multiple.columns <- function(db) {
	n <- 2
	results <- db$wsQuery(columns = c('id', 'entry name'), format = 'tab', limit = n, retfmt = 'parsed')
	testthat::expect_is(results, 'data.frame')
	testthat::expect_true(all(c('Entry', 'Entry name') %in% colnames(results)))
	testthat::expect_true(nrow(results) == n)
}

test.geneSymbolsToUniprotIds <- function(conn) {
    
    # Null list
    testthat::expect_equal(conn$geneSymbolToUniprotIds(NULL), list())
    testthat::expect_equal(conn$geneSymbolToUniprotIds(character()), list())
    testthat::expect_equal(conn$geneSymbolToUniprotIds(NA_character_), list())

    # Exact search
    expected_ids <- list('TGF-b1'='Q4LDM3')
    not_expected_ids <- list('TGF-b1'=c('W1I9X7', 'Q08FI9'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids))
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids)) {
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
    }

    # Exact search with two genes
    expected_ids <- list('TGF-b1'='Q4LDM3', 'G-CSF'='Q9GJU0')
    not_expected_ids <- list('TGF-b1'=c('W1I9X7', 'Q08FI9'),
                             'G-CSF'=c('P09919', 'Q8N4W3', 'A0A3G2Y4F6',
                                       'C0STS3'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids))
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids)) {
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
    }

    # Ignore non-alphanum chars (e.g.: "G-CSF" == "GCSF")
    expected_ids <- list('TGF-b1'='Q4LDM3', 'G-CSF'=c('Q9GJU0', 'P09919'))
    not_expected_ids <- list('TGF-b1'=c('W1I9X7', 'Q08FI9'),
                             'G-CSF'=c('Q8N4W3', 'A0A3G2Y4F6', 'C0STS3'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids),
                                       ignore.nonalphanum=TRUE)
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids)) {
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
    }
    
    # Partial match (e.g.: "G-CSFb1" matches when searching for "G-CSF")
    expected_ids <- list('TGF-b1'=c('Q4LDM3', 'W1I9X7'),
                         'G-CSF'=c('Q9GJU0', 'A0A3G2Y4F6', 'C0STS3'))
    not_expected_ids <- list('TGF-b1'=c('Q08FI9'),
                             'G-CSF'=c('P09919', 'Q8N4W3'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids), partial.match=TRUE)
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids)) {
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
        testthat::expect_false(any(not_expected_ids[[gene]] %in% ids[[gene]]))
    }

    # No filtering
    expected_ids <- list('TGF-b1'=c('Q4LDM3', 'W1I9X7'),
                         'G-CSF'=c('C0STS3', 'A0A679AQ73', 'Q8MKE0'))
    ids <- conn$geneSymbolToUniprotIds(names(expected_ids), filtering=FALSE)
    testthat::expect_equal(names(ids), names(expected_ids))
    for (gene in names(ids))
        testthat::expect_true(all(expected_ids[[gene]] %in% ids[[gene]]))
}

# Main
################################################################

# Instantiate Biodb
biodb <- biodb::createBiodbTestInstance(log='uniprot_test.log', ack=TRUE)

# Load package definitions
file <- system.file("definitions.yml", package='biodbUniprot')
biodb$loadDefinitions(file)

# Set context
biodb::setTestContext(biodb, "Test Uniprot connector.")

# Create connector
conn <- biodb$getFactory()$createConn('uniprot')

# Run tests
biodb::runGenericTests(conn)
biodb::testThat('Uniprot entries query works fine with an empty query.',
                test.uniprot.wsQuery.empty, conn=conn)
biodb::testThat('Uniprot entries query works fine with multiple columns',
                test.uniprot.wsQuery.multiple.columns, conn=conn)
biodb::testThat('Uniprot entries query works fine with a query by name.',
                test.uniprot.wsQuery.by.name, conn=conn)
biodb::testThat('We can convert gene symbols to UniProt IDs.',
                test.geneSymbolsToUniprotIds, conn=conn)

# Terminate Biodb
biodb$terminate()
