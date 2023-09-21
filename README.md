<img align="center" src="/img/logo_Artemis.png?" alt="logo" title="" width="150"/>

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

ARTEMIS provides an interface for utilizing a modified Temporal
Smith-Waterman (TSW) algorithm, derived from
[10.1109/DSAA.2015.7344785](https://www.researchgate.net/publication/292331949_Temporal_Needleman-Wunsch),
to summarize longitudinal EHR data into discrete regimen eras. Primarily
intended to be used for cancer patients, ARTEMIS utilizes data derived
from the [HemOnc](https://hemonc.org/wiki/Main_Page) oncology reference
to form the basic regimen data used in testing.

<figure>
<img src="/img/Workflow_Detailed.png?" alt="ARTEMIS Workflow" />
<figcaption aria-hidden="true">ARTEMIS Workflow</figcaption>
</figure>

## Features

ARTEMIS is primarily useful for stratifying patients based on their most
likely prescribed regimens, for use in cohort construction via the
Episode Era table of the [OMOP
CDM](https://www.ohdsi.org/data-standardization/).

ARTEMIS may also be used for providing summary statistics on the number
and distribution of regimens found within a specific cohort, as well as
their coverage and length, as well as providing summary graphics for
patient treatment trajectories.

<figure>
<img src="/img/Networks.png?" alt="Treatment Trajectories" />
<figcaption aria-hidden="true">Treatment Trajectories</figcaption>
</figure>

## Installation

ARTEMIS can presently be installed directly from GitHub:

    # install.packages("devtools")
    devtools::install_github("odyOSG/oncoRegimens")

ARTEMIS relies on a python back-end via reticulate and depending on your
reticulate settings, system and environment, you may need to run the
following commands before loading the package:

    #reticulate::py_install("numpy")
    #reticulate::py_install("pandas")

## Usage

### CDMConnector

A cdm\_reference object is created from any DBI connection, utilising
the [CDMConnector](https://darwin-eu.github.io/CDMConnector/) package.

    dbiconn <- DBI::dbConnect(RPostgres::Redshift(),
                              dbname = "dbName",
                              host = "hostName",
                              port = "9999",
                              user = "user",
                              password = "password")

    cdmSchema      <- "schema_containing_data"

    cdm <- CDMConnector::cdm_from_con(con = dbiconn,
                                      cdm_schema = cdmSchema,
                                      write_schema = "schema_with_write_access")

### Input

An input JSON containing a cohort specification is input by the user.
Information on OHDSI cohort creation and best practices can be found
[here](https://ohdsi.github.io/TheBookOfOhdsi/Cohorts.html).

    #loadCohort()
    json <- CDMConnector::readCohortSet(path = here::here("myCohort/"))
    name <- "myExampleCohort"

Regimen data may be read in from the provided package, or may be
submitted directly by the user. All of the provided regimens will be
tested against all patients within a given cohort.

    #loadRegimens()
    regimens <- read.csv(here::here("data/myRegimens.csv"))

A set of valid drugs may also be read in using the provided data, or may
be curated and submitted by the user. Only valid drugs will appear in
processed patient strings.

    #loadDrugs()
    validDrugs <- read.csv(here::here("data/myDrugs.csv"))

### Pipeline

The cdm connection is used to generate a dataframe containing the
relevant patient details for constructing regimen strings.

    con_df <- getCohortSet(cdm = cdm, json = json, name = name)

Regimen strings are then constructed, collated and filtered into a
stringDF dataframe containing all patients of interest.

    stringDF <- stringDF_from_cdm(con_df = con_df, writeOut = F, validDrugs = validDrugs)

    stringDF <- stringDF %>% filter_stringDF(min = 20)

The TSW algorithm is then run using user input settings and the provided
regimen and patient data. Detailed information on user inputs, such as
the gap penalty, g, can be found
[here](www.github.com/odyOSG/oncoRegimens)

    output_all <- stringDF %>% generateRawAlignments(regimens = regimens,
                                                     g = 0.4,
                                                     Tfac = 0.5,
                                                     verbose = 0,
                                                     mem = -1,
                                                     removeOverlap = 1,
                                                     method = "PropDiff")

Raw output alignments are then post-processed and may be visualised.
Post-processing steps include the handling and combination of
overlapping regimen alignments, as well as formatting output for
submission to an episode era table.

    processedAll <- output_all %>% processAlignments(regimenCombine = 28)

    personOfInterest <- output_all[output_all$personID == unique(output_all$personID)[1337],]

    plotOutput(personOfInterest, fontSize = 2.5)

Data may then be further explored via several graphics which indicate
various information, such as regimen frequency or the score/length
distributions of a given regimen.

    plotFrequency(output_processed)

    plotScoreDistribution(regimen1 = "Acetaminophen Monotherapy", regimen2 = "Ibuprofen Monotherapy", processedAll = output_processed)

    plotRegimenLengthDistribution(regimen1 = "Acetaminophen Monotherapy", regimen2 = "Ibuprofen Monotherapy", processedAll = output_processed)

Treatment trajectories, or regimen eras, can then be calculated, adding
further information about the relative sequencing order of different
regimens and regimen types.

    output_eras <- output_processed %>% calculateEras(discontinuationTime = 90)

    regStats <- output_eras %>% generateRegimenStats()

    regStats[,-c(4,7)]

And resulting graphics, such as a sankey indicating the overall patterns
of treatment trajectories can then be constructed.

    plotErasFrequency(output_eras)

    loadGroups()
    plotSankey(output_eras, regimen_Groups)

## DBI Drivers

CDMConnector is tested using the following DBI driver backends:

-   [RPostgres](https://rpostgres.r-dbi.org/reference/postgres) on
    Postgres and Redshift
-   [odbc](https://solutions.posit.co/connections/db/r-packages/odbc/)
    on Microsoft SQL Server, Oracle, and Databricks/Spark
-   [duckdb](https://duckdb.org/docs/api/r)

## Getting help

If you encounter a clear bug, please file an issue with a minimal
[reproducible example](https://reprex.tidyverse.org/) at the [GitHub
issues page](https://github.com/OdyOSG/oncoRegimens/issues).
