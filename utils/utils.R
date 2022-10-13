# Function to get plots for the power calculation
powerPlot = function(df, mafs) {
  
  a = ggplot(finaldf, aes(fill = beta, x = maf, y = power)) + 
    geom_bar(position = "dodge", stat = "identity") +
    facet_wrap(~SS, scales = "free") +
    theme_bw(base_size = 13) +
    geom_hline(yintercept=0.8, linetype='dashed', color=c('blue')) +
    scale_x_continuous(breaks= mafs )  +
    scale_y_continuous(expand = c(0.005,0.005), limits = c(0,1)) +
    scale_fill_discrete(name = "Effect Size") +
    ggtitle("Power to detect genetic associations in LMMs") +
    theme(panel.spacing = unit (2, "lines"),
          axis.title.x = element_text(margin = margin(t = 20), size = 40), 
          axis.title.y = element_text(margin = margin(r = 20), size = 40),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          plot.title = element_text(hjust = 0.5, size = 40, face = "bold", margin = margin(b = 20)),
          strip.text = element_text(size = 30),
          legend.text = element_text(size=40),
          legend.title = element_text(size = 30))
  
  
  return(a)
}


# Functions to get foresplots
myforestplot = function(df, snpname = "SNP", n_studies) {
  
  # Three studiesn = 3
  numbers = c(1, n_studies + 3, n_studies + 4) %>% as.character()
  mylist = setNames(vector(mode = "list", length = 3), numbers)
  
  lines = lapply(mylist, function(x) {
    x = gpar(lwd = 1, columns = 1:4, col = "#444444") 
  })
  
  df %>% 
    forestplot::forestplot(
      labeltext = c(study, effect, ci), 
      is.summary = summary,
      clip = c(-1, 2), 
      xlog = FALSE,
      graph.pos = 4,
      colgap = unit(0.03, "npc"),
      boxsize = 0.25,
      xticks=c(-1.5,-1, -0.5, 0, 0.5, 1, 1.5),
      #txt_gp = fpTxtGp(cex=0.75),
      hrzl_lines = lines,
      col = fpColors(box = "royalblue",
                     line = "darkblue",
                     summary = "royalblue"),
      xlab = "Effect size and CI",
      title = paste0(snpname, " forest plot"),
      txt_gp = fpTxtGp(ticks=gpar(cex=0.8),
                       xlab = gpar(col = "red", lty = "solid", lwd = 3, fontsize = 20)))
  
}





# Function to automate the imputation of missing data
impute_UPDRSIII <- function(df, ID = "ID", score_column = "score", visit_column = "visit_number",
                            UPDRSIII_column = "UPDRSIII_measure", UPDRS_type = "total",
                            missing_threshold = c(6,5,1), keep_vars = "Yes") {
  # libs loading
  library(hash)      # To create key-value pairs
  library(tidyverse) # For tidy evaluation and wrangling
  library(glue)      # To dynamically name variables
  
  # Safety checks
  if(!"score" %in% colnames(df)) stop("Warning: score column not found")
  #if(!"visit_number" %in% colnames(df)) stop("Warining: visit_number column not found")
  if(length(colnames) > 4) stop("Consider using tidyr::gather to get your data
                               in long format as input for this function")
  
  # vector with UPDRS_limb terms
  limb_terms <- c("III_3a", "III_3b", "III_3c", "III_3d", "III_3e", 
                  "III_15a", "III_15b", "III_16a","III_16b", "III_17a", "III_17b",
                  "III_17c", "III_17d", "III_17e", "III_4a", "III_4b",
                  "III_5a", "III_5b", "III_6a", "III_6b", "III_7a", "III_7b",
                  "III_8a", "III_8b", "III_18",
                  "3_3a", "3_3b", "3_3c", "3_3d", "3_3e", 
                  "3_15a", "3_15b", "3_16a","3_16b", "3_17a", "3_17b",
                  "3_17c", "3_17d", "3_17e", "3_4a", "3_4b",
                  "3_5a", "3_5b", "3_6a", "3_6b", "3_7a", "3_7b",
                  "3_8a", "3_8b", "3_18",
                  "NP3RIGN", "NP3RIGRU", "NP3RIGLU", "NP3RIGRL", "NP3RIGLL",
                  "NP3PTRMR", "NP3PTRML", "NP3KTRMR", "NP3KTRML", "NP3RTARU", "NP3RTALU",
                  "NP3RTARL", "NP3RTALL", "NP3RTALJ", "NP3FTAPR", "NP3FTAPL",
                  "NP3HMOVR", "NP3HMOVL", "NP3PRSPR", "NP3PRSPL", 
                  "NP3TTAPR", "NP3TTAPL", "NP3LGAGR", "NP3LGAGL",
                  "NP3RTCON",
                  "RigidityNeck", "RigidityRightUpper", "RigidityRightLower", "RigidityLeftUpper", "RigidityLeftLower",
                  "PosturalTremorHandsRight", "PosturalTremorHandsLeft",
                  "KineticTremorHandsRight", "KineticTremorHandsLeft",
                  "RestTremorAmplitudeRightUpper", "RestTremorAmplitudeRightLower", "RestTremorAmplitudeLeftUpper", 
                  "RestTremorAmplitudeLeftLower", "RestTremorAmplitudeLipJaw",
                  "FingerTappingRight","FingerTappingLeft",
                  "HandMovementsRight", "HandMovementsLeft", "PSHandsRight", "PSHandsLeft",
                  "ToeTappingRight", "ToeTappingLeft", "LegAgilityRight", "LegAgilityLeft",
                  "ConstancyOfRestTremor",
                  "MOTOREXAM3NECK", "MOTOREXAM3RUE", "MOTOREXAM3LUE", "MOTOREXAM3RLE", "MOTOREXAM3LLE",
                  "MOTOREXAM15R", "MOTOREXAM15L", "MOTOREXAM16R", "MOTOREXAM16L", 
                  "MOTOREXAM17LIPJ", "MOTOREXAM17RUE", "MOTOREXAM17LUE", "MOTOREXAM17RLE", "MOTOREXAM17LLE",
                  "MOTOREXAM4R", "MOTOREXAM4L", "MOTOREXAM5R", "MOTOREXAM5L", "MOTOREXAM6R", "MOTOREXAM6L",
                  "MOTOREXAM7R", "MOTOREXAM7L", "MOTOREXAM8R", "MOTOREXAM8L", "MOTOREXAM18")
  
  # Vector with UPDRS_axial terms
  axial_terms <- c("III_1", "III_2", "III_9", "III_10", "III_11", "III_12", "III_13", "III_14", 
                   "3_1", "3_2", "3_9", "3_10", "3_11", "3_12", "3_13", "3_14",
                   "NP3SPCH", "NP3FACXP", "NP3RISNG", "NP3GAIT", "NP3FRZGT", "NP3PSTBL", 
                   "NP3POSTR", "NP3BRADY",
                   "Speech", "FacialExpression", "ArisingFromChair", "Gait", "FreezingOfGait",
                   "PosturalStability", "Posture", "GSM",
                   "MOTOREXAM1", "MOTOREXAM2", "MOTOREXAM9", "MOTOREXAM10", "MOTOREXAM11", "MOTOREXAM12", "MOTOREXAM13", "MOTOREXAM14")
  
  
  # Get all vars to merge later
  df_allvars <- df %>% select(!c(.data[[score_column]], .data[[UPDRSIII_column]])) %>%
    distinct(.data[[ID]], .data[[visit_column]], .keep_all = TRUE)
  
  
  # Removing UPDRSIII_total measure if present
  df <- df %>%
    dplyr::select(!ends_with(c("UPDRS_III", "total")))
  
  final_df = df %>% dplyr::select(.data[[ID]], .data[[visit_column]]) %>% 
    distinct()
  
  for (type_index in 1:length(UPDRS_type)) {
    
    # To be used when naming the imputed columns
    term = UPDRS_type[type_index]
    if (term == "total") {
      mydf <- df %>% 
        dplyr::select(.data[[ID]], .data[[visit_column]], .data[[score_column]], contains(c("UPDRS")))
      
      mydf <- mydf %>%
        dplyr::group_by(.data[[ID]], .data[[visit_column]]) %>% 
        dplyr::mutate(missing_total = sum(is.na(.data[[score_column]])),
                      UPDRSIIItotal_imputed = ifelse(missing_total > missing_threshold[1], NA,
                                                     (sum(score, na.rm = TRUE)) / (33-missing_total) * 33 )) 
      
      # UPDRSIII_limb imputation
    } else if (term == "limb") {
      
      mydf <- df %>%
        dplyr::select(.data[[ID]], .data[[score_column]], .data[[visit_column]], contains(c("UPDRS"))) %>%
        dplyr::filter(grepl(paste(limb_terms, collapse = "|"), .data[[UPDRSIII_column]]))
      
      mydf <- mydf %>%
        dplyr::group_by(.data[[ID]], .data[[visit_column]]) %>% 
        dplyr::mutate(missing_total = sum(is.na(.data[[score_column]])),
                      UPDRSIIIlimb_imputed = ifelse(missing_total > missing_threshold[2], NA,
                                                    (sum(score, na.rm = TRUE)) / (25-missing_total) * 25))
      # UPDRSIII_axial imputation
    } else if (term == "axial") {
      
      mydf <- df %>%
        dplyr::select(.data[[ID]], .data[[score_column]], .data[[visit_column]], contains(c("UPDRS"))) %>%
        dplyr::filter(!grepl(paste(limb_terms, collapse = "|"), .data[[UPDRSIII_column]]))
      
      mydf <- mydf %>%
        dplyr::group_by(.data[[ID]], .data[[visit_column]]) %>% 
        dplyr::mutate(missing_total = sum(is.na(.data[[score_column]])),
                      UPDRSIIIaxial_imputed = ifelse(missing_total > missing_threshold[3], NA,
                                                     (sum(score, na.rm = TRUE)) / (8-missing_total) * 8)) 
      
    } else {
      stop("Only UPDRS_total, UPDRS_limb, or UPDRS_axial supported")
    }
    
    measure_name = paste0("UPDRSIII_measure_", term)
    mydf <- mydf %>%
      dplyr::arrange(.data[[ID]], .data[[visit_column]]) %>%
      dplyr::filter(row_number() == 1) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate("{measure_name}" := paste0("V", .data[[visit_column]], "_UPDRS_III_", term)) %>%
      dplyr::select(-c(UPDRSIII_measure, score, missing_total))
    
    # final_df <- final_df %>% 
    #   dplyr::left_join(mydf, by = c(.data[[ID]], .data[[visit_column]]))
    
    final_df <- final_df %>% 
      dplyr::left_join(mydf, by = c(ID, visit_column))
    
  }
  
  if (keep_vars == "Yes") {
    final_df = df_allvars %>% inner_join(final_df)
  }
  
  return(final_df)
}
