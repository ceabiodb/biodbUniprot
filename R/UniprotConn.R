#' The connector class to Uniprot database.
#'
#' This is a concrete connector class. It must never be instantiated directly,
#' but instead be instantiated through the factory BiodbFactory.
#' Only specific methods are described here. See super classes for the
#' description of inherited methods.
#'
#' @examples
#' # Create an instance with default settings:
#' mybiodb <- biodb::newInst()
#'
#' # Get Uniprot connector
#' uniprot <- mybiodb$getFactory()$createConn('uniprot')
#'
#' # Access web service query
#' result <- uniprot$wsQuery(query='name:"prion protein"',
#'                            columns=c('id', 'entry name'),
#'                            format='txt', limit=10)
#'
#' # Terminate instance.
#' mybiodb$terminate()
#'
#' @import R6
#' @export
UniprotConn <- R6::R6Class("UniprotConn",
inherit=biodb::BiodbConn,

public=list(

#' @description
#' Calls query to the database for searching for compounds.
#' See http //www.uniprot.org/help/api_queries for details.
#' @param query The query to send to the database.
#' @param columns The field columns to retrieve from the database (e.g.: 'id',
#' 'entry name', 'pathway', 'organism', 'sequence', etc).
#' @param format The return format (e.g.: 'tab').
#' @param limit The maximum number of entries to return.
#' @param retfmt Use to set the format of the returned value. 'plain' will
#' return the raw results from the server, as a character value. 'parsed' will
#' return the parsed results, as a JSON object. 'request' will return a
#' BiodbRequest object representing the request as it would have been sent.
#' 'ids' will return a character vector containing the IDs of the matching
#' entries.
#' @return Depending on `retfmt` parameter.
wsQuery=function(query='', columns=NULL, format=NULL, limit=NULL,
    retfmt=c('plain', 'parsed', 'ids', 'request')) {

    # Check parameters
    retfmt <- match.arg(retfmt)
    if (retfmt == 'ids') {
        columns <- 'id'
        format <- 'tab'
    }
    if (is.null(format) || is.na(format))
        format <- 'tab'
    if (retfmt != 'ids' && (is.null(columns) || all(is.na(columns))))
        columns <- c("citation", "clusters", "comments", "domains", "domain",
            "ec", "id", "entry name", "existence", "families", "features",
            "genes", "go", "go-id", "interactor", "keywords", "last-modified",
            "length", "organism", "organism-id", "pathway", "protein names",
            "reviewed", "sequence", "3d", "version", "virus hosts")
    columns <- paste(columns, collapse=',')

    # Build request
    params <- list(query=query, columns=columns, format=format)
    if ( ! is.null(limit) && ! is.na(limit))
        params[['limit']] <- limit
    url <- BiodbUrl$new(url=c(self$getPropValSlot('urls', 'base.url'), ''),
        params=params)
    request <- self$makeRequest(method='get', url=url)

    # Return request
    if (retfmt == 'request')
        return(request)

    # Send request
    results <- self$getBiodb()$getRequestScheduler()$sendRequest(request)
    results <- unlist(results)

    # Parse
    if (retfmt != 'plain') {

        # Parse data frame
        readtc <- textConnection(results, "r", local=TRUE)
        df <- read.table(readtc, sep="\t", header=TRUE, check.names=FALSE)
        close(readtc)
        results <- df

        # Get IDs
        if (retfmt == 'ids')
            results <- as.character(results[[1]])
    }

    return(results)
},

#' @description
#' Gets UniProt IDs associated with gene symbols.
#' @param genes A vector of gene symbols to convert to UniProt IDs.
#' @param ignore.nonalphanum If set to TRUE, do not take into account
#' non-alphanumeric characters when comparing gene symbols. 
#' @param partial.match If set to TRUE, a match will be valid even if the
#' provided gene symbol is only a substring of the found gene symbol.
#' @param filtering If set to FALSE, do not run any filtering and return all
#' the UniProt IDs given by UniProt Query web service.
#' @param max.results Maximum of UniProt IDs returned for each gene symbol.
#' @return A named list of vectors of UniProt IDs. The names are gene
#' symbols provided with the genes parameter. For each gene symbol, a vector
#' of found UniProt IDs is set.
geneSymbolToUniprotIds=function(genes, ignore.nonalphanum=FALSE,
    partial.match=FALSE, filtering=TRUE, max.results=0) {
    
    ids <- list()
    
    if ( ! is.null(genes) && length(genes) > 0) {
        
        for (gene in genes)
            if ( ! is.na(gene)) {
                
                # Set limit
                limit <- NULL
                if ( ! filtering && max.results > 0)
                    limit <- max.results
                
                # Prepare query
                wanted_genes <- gene
                if (filtering && ignore.nonalphanum)
                    wanted_genes <- c(wanted_genes, gsub('[^A-Za-z0-9]', '',
                        gene))
                if (filtering && partial.match)
                    wanted_genes <- c(wanted_genes, paste0('*', wanted_genes),
                        paste0(wanted_genes, '*'))
                query <- paste(paste0('gene:', wanted_genes), collapse=' OR ')
                
                # Run query
                x <- self$wsQuery(query, retfmt='ids', limit=limit)
                
                # Filtering
                # Needed since UniProt web service may return entries with
                # gene symbols like "TGF-b1a" when asking for "TGF-b1"
                if (filtering)
                    x <- private$filterResults(x, gene=gene,
                        ignore.nonalphanum=ignore.nonalphanum,
                        partial.match=partial.match, limit=max.results)

                # Cut
                if ( max.results > 0 && length(x) > max.results)
                    x <- x[seq_len(max.results)]

                # Set results
                ids[[gene]] <- x
            }
    }

    return(ids)
}
),

private=list(

doGetEntryPageUrl=function(id) {

    u <- c(self$getPropValSlot('urls', 'base.url'), id)
    f <- function(x) BiodbUrl$new(url=u)$toString()
    return(vapply(id, f, FUN.VALUE=''))
}

,doSearchForEntries=function(fields=NULL, max.results=0) {

    query <- ''

    # Search by name
    if ('name' %in% names(fields)) {
        name.query <- paste('name', paste0('"', fields$name, '"'), sep=':')
        mnemonic.query <- paste('mnemonic', paste0('"', fields$name, '"'),
                                sep=':')
        query <- paste(name.query, mnemonic.query, sep=' OR ')
    }

    # Search by mass
    if ('molecular.mass' %in% names(fields)) {

        rng <- do.call(Range$new, fields[['molecular.mass']])

        # Uniprot does not accept mass in floating numbers
        uniprot.mass.min <- as.integer(rng$getMin())
        uniprot.mass.max <- as.integer(rng$getMax())
#        if (uniprot.mass.min != mass.min || uniprot.mass.max != mass.max)
#            biodb::warn0('Uniprot requires integers for mass range.',
#                          ' Range [', mass.min, ', ', mass.max,
#                          '] will be converted into [', uniprot.mass.min,
#                          ', ', uniprot.mass.max, '].')

        mass.query <- paste0('mass:[', uniprot.mass.min, ' TO ',
            uniprot.mass.max, ']')

        if (nchar(query) > 0) {
            query <- paste0('(', query, ')')
            query <- paste(query, mass.query, sep=' AND ')
        }
        else
            query <- mass.query
    }

    # Send query
    ids <- self$wsQuery(query=query, limit=max.results, retfmt='ids')

    return(ids)
}

,doGetEntryContentRequest=function(id, concatenate=TRUE) {

    url <- paste0(self$getPropValSlot('urls', 'base.url'), id, '.xml')

    return(url)
}

,doGetEntryIds=function(max.results=0) {

    ids <- self$wsQuery(limit=max.results, retfmt='ids')

    return(ids)
}

,filterResults=function(x, gene, ignore.nonalphanum, partial.match, limit=0) {

    # Prepare gene to find
    gene.to.find <- gene
    if (ignore.nonalphanum)
        gene.to.find <- gsub('[^A-Za-z0-9]', '', gene.to.find)

    # Ignore case
    gene.to.find <- tolower(gene.to.find)

    # Filtering function
    gene_matches <- function(id) {
        
        matches <- FALSE
        
        # Get gene symbols and prepare them
        entry <- self$getEntry(id)

        if ( ! is.null(entry)) {

            gene.symbols <- entry$getFieldValue('gene.symbol')
            if (ignore.nonalphanum)
                gene.symbols <- gsub('[^A-Za-z0-9]', '', gene.symbols)
            # Ignore case
            gene.symbols <- tolower(gene.symbols)
            
            if (partial.match)
                matches <- length(grep(gene.to.find, gene.symbols,
                    fixed=TRUE)) > 0
            else
                matches <- gene.to.find %in% gene.symbols
        }
        
        return(matches)
    }
    
    # Create progress instance
    prg <- Progress$new(biodb=self$getBiodb(), msg='Filtering results.',
        total=length(x))
    
    # Run filtering
    y <- character()
    for (id in x) {

        if (gene_matches(id))
            y <- c(y, id)

        if (limit > 0 && length(y) == limit)
            break

        # Progress message
        prg$increment()
    }

    return(y)
}
))
