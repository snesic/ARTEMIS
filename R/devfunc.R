library(tidyverse)
library(RColorBrewer)

# Test run 
# Development tools functionality
#' @export
full_run_artemis <- function(stringDF, regimens, g = 0.4) {
    
    #### RUN ARTEMIS
    output_all <- generateRawAlignments(stringDF,
                                        regimens = regimens,
                                        g = g,
                                        Tfac = 0.5,
                                        verbose = 0,
                                        mem = -1,
                                        removeOverlap = 1,
                                        method = "PropDiff")
    
    processedAlignments <- processAlignments(output_all,
                                             regimenCombine = 28,
                                             regimens = regimens)
    
    processedAlignments
}


#### PATCHED ERAS FUNCTION, PREVIOUS WAS BROKEN
#' @export
calculateEras_patched <- function(processedAll, discontinuationTime = 120){
    IDs_All <- unique(processedAll$personID)
    processedEras <- processedAll[0, ]
    
    for (i in c(1:length(IDs_All))) {
        tempDF <- processedAll[processedAll$personID == IDs_All[i],]
        tempDF <- tempDF[order(tempDF$t_start), ]
        toRemove <- c()
        if (dim(tempDF)[1] > 1) {
            for (i in c(2:length(tempDF$component))) {
                if (tempDF[i, ]$t_start < tempDF[i - 1, ]$t_end) {
                    toRemove <- c(toRemove, i)
                }
            }
        }
        
        if (length(toRemove) > 0) {
            tempDF <- tempDF[-toRemove, ]
        }
        
        tempDF <- tempDF %>% dplyr::mutate(timeToNextRegimen = dplyr::lead(t_start) -
                                               t_end)
        tempDF <- tempDF %>% dplyr::mutate(lag = dplyr::lag(.data$timeToNextRegimen),
                                           delete = ifelse((dplyr::lag(.data$timeToNextRegimen) <
                                                                discontinuationTime & component == dplyr::lag(component)),
                                                           "Y", "N"))
        tempDF[1, ]$delete <- "N"
        
        tempDF <- tempDF %>%
            mutate(newLine = cumsum(delete == "N")) %>%
            summarise(t_start = min(t_start), t_end = max(t_end), timToEod = min(timeToEOD), .by = c(component, newLine, personID)) %>%
            mutate(regLength = t_end - t_start,
                   timeToNextRegimen = lag(t_start,1) - t_end,
                   First_Line = 1 * (row_number() == 1),
                   Second_Line = 1 * (row_number() == 2),
                   Other = 1 * (row_number() > 2))
        
        processedEras <- rbind(processedEras, tempDF)
        
        #Handle overlapping regimens
        processedEras$timeToNextRegimen[processedEras$timeToNextRegimen < 0] <- 0
        
        
    }
    
    return(processedEras)
    
}




plot_alignment_for_patient <- function(stringDF_subset, processedAlignments) {
    
    df = stringDF_subset %>% 
        separate_rows(seq) %>% 
        separate(seq, into = c("time", 'component')) %>% 
        filter(time != "") %>% 
        group_by(person_id) %>% 
        mutate(time_start = cumsum(as.integer(time)),
               time_end = time_start + 1,
               #time_cumsum = cumsum(time),
               case = "patient", 
               person_id = as.character(person_id)) %>% 
        arrange(time)
    
    df = processedAlignments %>% 
        select(person_id = personID, 
               component, 
               time_start = t_start, 
               time_end = t_end) %>% 
        mutate(case = "aligned") %>% 
        bind_rows(df) %>% 
        mutate(component = fct_reorder(component, time_start))
        
    
    
    
    # Get unique components for WT and CT
    patient_components <- unique(df$component[df$case == "patient"])
    aligned_components <- unique(df$component[df$case == "aligned"])
    
    patient_components = as.character(patient_components)
    aligned_components = as.character(aligned_components
                                      )
    # Generate dynamic color palettes
    patient_colors <- setNames(brewer.pal(length(patient_components), "Set1"), patient_components)
    aligned_colors <- setNames(brewer.pal(length(aligned_components), "Set2"), aligned_components)
    # Combine color mappings
    colors <- c(patient_colors, aligned_colors)
    # Create separate aesthetics for WT and CT
    df$patient_components_col <- ifelse(df$case == "patient", as.character(df$component), NA)
    df$aligned_components_col <- ifelse(df$case == "aligned", as.character(df$component), NA)
    
    
    p = df %>%    
        ggplot() +
        geom_segment(aes(x = time_start, xend = time_end, y = component, 
                         yend = component, color = patient_components_col), 
                     linewidth = 2, na.rm = TRUE) +
        geom_segment(aes(x = time_start, xend = time_end, y = component, 
                         yend = component, color = aligned_components_col), 
                     linewidth = 2, na.rm = TRUE) +    
        geom_point(aes(x = time_start, y = component, color = patient_components_col)) +
        facet_grid(cols = vars(person_id), rows = vars(case), scale = "free_y") +
        scale_color_manual(name = "patient", values = colors, na.translate = FALSE, 
                           guide = guide_legend(order = 1)) +
        scale_color_manual(name = "aligned", values = colors, na.translate = FALSE, 
                           guide = guide_legend(order = 3, size = 3)) +
        labs(x = "Time", y = "Component", title = "Time Intervals per Component") +
        theme_bw() +
        theme(legend.position = "none") + # Move legends below 
        ggtitle(label = paste("Patient", unique(df$person_id)))
    
    return(p)
}
