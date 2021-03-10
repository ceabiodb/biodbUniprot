#' The connector class to Uniprot database.
#'
#' This is a concrete connector class. It must never be instantiated directly,
#' but instead be instantiated through the factory BiodbFactory.
#' Only specific methods are described here. See super classes for the
#' description of inherited methods.
#'
#' @examples
#' # Create an instance with default settings:
#' mybiodb <- biodb::Biodb()
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
#' @import methods
#' @export UniprotConn
#' @exportClass UniprotConn
UniprotConn <- methods::setRefClass("UniprotConn",
    contains=c("BiodbRemotedbConn", "BiodbCompounddbConn"),

methods=list(

wsQuery=function(query='', columns=NULL, format=NULL, limit=NULL,
                 retfmt=c('plain', 'parsed', 'ids', 'request')) {
    ":\n\nCalls query to the database for searching for compounds.
    See http://www.uniprot.org/help/api_queries for details.
    \nquery: The query to send to the database.
    \ncolumns: The field columns to retrieve from the database (e.g.: 'id',
    'entry name', 'pathway', 'organism', 'sequence', etc).
    \nformat: The return format (e.g.: 'tab').
    \nlimit: The maximum number of entries to return.
    \nretfmt: Use to set the format of the returned value. 'plain' will return
    the raw results from the server, as a character value. 'parsed' will return
    the parsed results, as a JSON object. 'request' will return a BiodbRequest
    object representing the request as it would have been sent. 'ids' will
    return a character vector containing the IDs of the matching entries.
    \nReturned value: Depending on `retfmt` parameter.
    "

    retfmt <- match.arg(retfmt)

    # Set parameters for retrieving IDs
    if (retfmt == 'ids') {
        columns <- 'id'
        format <- 'tab'
    }

    # Set columns
    if (is.null(columns) || all(is.na(columns)))
        columns <- c("citation", "clusters", "comments", "domains", "domain",
                     "ec", "id", "entry name", "existence", "families",
                     "features", "genes", "go", "go-id", "interactor",
                     "keywords", "last-modified", "length", "organism",
                     "organism-id", "pathway", "protein names", "reviewed",
                     "sequence", "3d", "version", "virus hosts")
    columns <- paste(columns, collapse=',')

    # Set format
    if (is.null(format) || is.na(format))
        format <- 'tab'

    # Build request
    params <- list(query=query, columns=columns, format=format)
    if ( ! is.null(limit) && ! is.na(limit))
        params[['limit']] <- limit
    url <- BiodbUrl(url=c(.self$getPropValSlot('urls', 'base.url'), ''),
                    params=params)
    request <- .self$makeRequest(method='get', url=url)

    # Return request
    if (retfmt == 'request')
        return(request)

    # Send request
    results <- .self$getBiodb()$getRequestScheduler()$sendRequest(request)

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

geneSymbolToUniprotIds=function(genes, ignore.case=TRUE,
                                ignore.nonalphanum=FALSE, partial.match=FALSE,
                                filtering=TRUE, max.results=NA_integer_) {
    ":\n\nGets UniProt IDs associated with gene symbols.
    \ngenes: A vector of gene symbols to convert to UniProt IDs.
    \nignore.case: If set to TRUE, ignore character case when comparing gene
    symbols.
    \nignore.nonalphanum: If set to TRUE, do not take into account
    non-alphanumeric characters when comparing gene symbols. 
    \npartial.match: If set to TRUE, a match will be valid even if the provided
    gene symbol is only a substring of the found gene symbol.
    \nfiltering: If set to FALSE, do not run any filtering and return all the
    UniProt IDs given by UniProt Query web service.
    \nmax.results: Maximum of UniProt IDs returned for each gene symbol.
    \nReturned value: A named list of vectors of UniProt IDs. The names are gene
    symbols provided with the genes parameter. For each gene symbol, a vector
    of found UniProt IDs is set.
    "
    
    ids <- list()
    
    if ( ! is.null(genes) && length(genes) > 0) {
        
        for (gene in genes)
            if ( ! is.na(gene)) {
                
                # Get UniProt IDs
                limit <- NULL
                if ( ! filtering && ! is.na(max.results))
                    limit <- max.results
                x <- .self$wsQuery(gene, columns='genes', retfmt='ids',
                                   limit=limit)
                
                # Filters
                if (filtering) {
                    
                    # Prepare gene to find
                    gene.to.find <- gene
                    if (ignore.nonalphanum)
                        gene.to.find <- gsub('[^A-Za-z0-9]', '', gene.to.find)
                    if (ignore.case)
                        gene.to.find <- tolower(gene.to.find)
                    
                    # Filtering function
                    gene_matches <- function(id) {
                        
                        matches <- FALSE
                        
                        # Get gene symbols and prepare them
                        entry <- .self$getEntry(id)
                        gene.symbols <- entry$getFieldValue('gene.symbol')
                        if (ignore.nonalphanum)
                            gene.symbols <- gsub('[^A-Za-z0-9]', '',
                                                 gene.symbols)
                        if (ignore.case)
                            gene.symbols <- tolower(gene.symbols)
                        
                        if (partial.match)
                            matches <- length(grep(gene.to.find, gene.symbols,
                                                   fixed=TRUE)) > 0
                        else
                            matches <- gene.to.find %in% gene.symbols
                        
                        return(matches)
                    }
                    
                    # Run filtering
                    x <- Filter(gene_matches, x)
                }

                # Cut
                if ( ! is.na(max.results) && length(x) > max.results)
                    x <- x[1:max.results]

                # Set results
                ids[[gene]] <- x
            }
    }

    return(ids)
},

getEntryPageUrl=function(id) {
    # Overrides super class' method.

    u <- c(.self$getPropValSlot('urls', 'base.url'), id)
    f <- function(x) BiodbUrl(url=u)$toString()
    return(vapply(id, f, FUN.VALUE=''))
},

searchCompound=function(name=NULL, mass=NULL, mass.field=NULL, mass.tol=0.01,
                        mass.tol.unit='plain', max.results=NA_integer_) {
    # Overrides super class' method.

    .self$.checkMassField(mass=mass, mass.field=mass.field)

    query <- ''

    # Search for name
    if ( ! is.null(name) && ! is.na(name)) {
        name.query <- paste('name', paste0('"', name, '"'), sep=':')
        mnemonic.query <- paste('mnemonic', paste0('"', name, '"'), sep=':')
        query <- paste(name.query, mnemonic.query, sep=' OR ')
    }

    # Search for mass
    if ( ! is.null(mass) && ! is.null(mass.field)) {

        mass.field <- .self$getBiodb()$getEntryFields()$getRealName(mass.field)

        if (mass.field != 'molecular.mass')
            .self$caution('Mass field "', mass.field, '" is not handled.')

        else {

            if (mass.tol.unit == 'ppm') {
                mass.min <- mass * (1 - mass.tol * 1e-6)
                mass.max <- mass * (1 + mass.tol * 1e-6)
            } else {
                mass.min <- mass - mass.tol
                mass.max <- mass + mass.tol
            }

            # Uniprot does not accept mass in floating numbers
            uniprot.mass.min <- as.integer(mass.min)
            uniprot.mass.max <- as.integer(mass.max)
            if (uniprot.mass.min != mass.min || uniprot.mass.max != mass.max)
                .self$caution('Uniprot requires integers for mass range.',
                              ' Range [', mass.min, ', ', mass.max,
                              '] will be converted into [', uniprot.mass.min,
                              ', ', uniprot.mass.max, '].')

            mass.query <- paste0('mass:[', uniprot.mass.min, ' TO ',
                                 uniprot.mass.max, ']')

            if (nchar(query) > 0) {
                query <- paste0('(', query, ')')
                query <- paste(query, mass.query, sep=' AND ')
            }
            else
                query <- mass.query
        }
    }

    # Send query
    ids <- .self$wsQuery(query=query, limit=max.results, retfmt='ids')

    return(ids)
},

.doGetEntryContentRequest=function(id, concatenate=TRUE) {

    url <- paste0(.self$getPropValSlot('urls', 'base.url'), id, '.xml')

    return(url)
},

.doGetEntryIds=function(max.results=NA_integer_) {

    ids <- .self$wsQuery(limit=max.results, retfmt='ids')

    return(ids)
}

))
