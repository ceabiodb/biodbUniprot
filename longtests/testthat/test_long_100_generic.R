# Set test context
biodb::testContext("Generic long tests")

# Instantiate Biodb
biodb <- biodb::createBiodbTestInstance(ack=TRUE)

# Load package definitions
defFile <- system.file("definitions.yml", package='biodbUniprot')
biodb$loadDefinitions(defFile)

# Create connector
conn <- biodb$getFactory()$createConn('uniprot')

# Run generic tests
testRefFolder <- system.file("testref", package='biodbUniprot')
biodb::runGenericTests(conn, pkgName='biodbUniprot', short=FALSE, long=TRUE,
    testRefFolder=testRefFolder, opt=list(max.results=3))

# Terminate Biodb
biodb$terminate()
