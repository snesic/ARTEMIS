#
Sys.setenv(RETICULATE_PYTHON = "/Users/snesic/miniconda3/envs/tsw-env/bin/python3")

devtools::load_all("../ARTEMIS/")

##### INPUT #####
# TODO: check if we want to keep this as an example:
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                server="server/database",
                                                                user="username",
                                                                port = "9999",
                                                                password="password",
                                                                pathToDriver = "./JBDC")

# Load SQLite test database within the package
db_path <- system.file("extdata", "testing_db.sqlite", package = "ARTEMIS")

connectionDetails <- DatabaseConnector::createConnectionDetails(
    dbms = "sqlite",
    server = db_path
)


df_json <- loadCohort()
name <- "lungcancer"

validdrugs <- loadDrugs()
regimens <- loadRegimens()
regGroups <- loadGroups()

cdmSchema      <- "main"
writeSchema    <- "main"

##### MAIN #####
con_df <- getConDF(connectionDetails = connectionDetails, 
                   json = df_json$json, 
                   name = name, 
                   cdmSchema = cdmSchema, 
                   writeSchema = writeSchema)


# Check if the dates are correctly written, and if not: 
con_df$drug_exposure_start_date
con_df$drug_exposure_start_date <- as.POSIXct(con_df$drug_exposure_start_date, 
                                              origin = "1970-01-01", tz = "UTC")

stringDF <- stringDF_from_cdm(con_df = con_df, writeOut = F, validDrugs = validdrugs)

## Alignment

output_all <- stringDF %>% generateRawAlignments(regimens = regimens,
                                                 g = 0.4,
                                                 Tfac = 0.4,
                                                 method = "PropDiff",
                                                 verbose = 0)

## Post-process Alignment

processedAll <- output_all %>% 
    processAlignments(regimenCombine = 28, regimens = regimens)

processedEras <- processedAll %>% calculateEras()

# Regimen stats does not work 
# regStats <- processedEras %>% generateRegimenStats()

##### OUTPUT #####

writeOutputs(output_all, processedAll = processedAll, processedEras = processedEras,
             connectionDetails = connectionDetails, cdmSchema = cdmSchema, regGroups = regGroups,
             regStats = regStats, stringDF = stringDF, con_df = con_df)
