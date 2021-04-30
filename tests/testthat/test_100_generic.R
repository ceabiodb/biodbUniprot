# Set test context
biodb::testContext("Generic tests")

# Instantiate Biodb
biodb <- biodb::createBiodbTestInstance(ack=TRUE)

# Load package definitions
defFile <- system.file("definitions.yml", package='biodbUniprot')
biodb$loadDefinitions(defFile)

# Create connector
conn <- biodb$getFactory()$createConn('uniprot')

# Run generic tests
biodb::runGenericTests(conn, list(max.results=3))

# Terminate Biodb
biodb$terminate()
