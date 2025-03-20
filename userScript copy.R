#
library(ARTEMIS)
library(tidyverse)

# ##### INPUT #####
# connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
#                                                                 server="server/database",
#                                                                 user="username",
#                                                                 port = "9999",
#                                                                 password="password",
#                                                                 pathToDriver = "./JBDC")

# json <- loadCohort()
stringDF <- read_csv("/path/to/stringDF_examples.csv")


name <- "lungcancer"

validDrugs <- loadDrugs(
  absolute="/Users/home/prj/art-ohdsi/etl_clean/prod/results/validdrugs_nsclc_20mar2025.rda"
  )

# regimens <- loadRegimens(condition = "lungCancer") 
#
# **With Update**
#
regimens <- loadRegimens(
    condition = "lungCancer", 
    absolute="/Users/home/prj/art-ohdsi/etl_clean/prod/results/regimens_nsclc_20mar2025.rda"
    )


regGroups <- loadGroups()

cdmSchema      <- "schema_containing_data"
writeSchema    <- "schmea_with_write_access"


##### MAIN #####
full_run_artemis(stringDF, regimens)

output_all = read_csv("Output.csv")
processedEras = read_csv("Output_Processed.csv")

# con_df <- getConDF(connectionDetails = connectionDetails, json = json, name = name, cdmSchema = cdmSchema, writeSchema = writeSchema)

# # stringDF <- stringDF_from_cdm(con_df = con_df, writeOut = F, validDrugs = validdrugs)
# stringDF <- stringDF_from_cdm(con_df = con_df, writeOut = F, validDrugs = validdrugs)

# output_all <- stringDF %>% generateRawAlignments(regimens = regimens,
# #                                                  g = 0.4,
# #                                                  Tfac = 0.4,
# #                                                  method = "PropDiff",
# #                                                  verbose = 0)

# processedAll <- output_all %>% processAlignments(regimenCombine = 28, regimens = regimens)

# processedEras <- processedAll %>% calculateEras()

# regStats <- processedEras %>% generateRegimenStats()

# ##### OUTPUT #####

# writeOutputs(output_all, processedAll = processedAll, processedEras = processedEras,
#              connectionDetails = connectionDetails, cdmSchema = cdmSchema, regGroups = regGroups,
#              regStats = regStats, stringDF = stringDF, con_df = con_df)
