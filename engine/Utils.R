#' -----------------------------------------------------------------------------
#' Script to hold utils (common fuctions). 
#' by Ham - May 2022
#'
#' -----------------------------------------------------------------------------
#' Functions:
#'   1. getcbPal - returns colorblind palette
#'   2. parse_fname - returns fname_df
#'   3. org_Behav_names - Used to assign and correct behaviour names
#'   4. create_facet_design - returns manual facet design for use with ggh4x::facet_manual()
#'
#' -----------------------------------------------------------------------------
#' 
#'
#' -----------------------------------------------------------------------------

getcbPal <- function(set = "Ori", black = T, show = F) {
  #' Returns colorblind palette for specific sets
  #' 
  #' @param set string, Ori: Returns palette specific to categories - c("Ori", "USV", "Behav").
  #' @param black bool, True: Returns palette w black (True) or grey (False).
  #' @param show bool, False: Shows palette if True.
  #'
  #' @return A vector of n colour hex codes
  
  if (set == "Ori") {
    cat("Returning color palette for", set, fill = T)
    if (black) {
      # The palette with black:
      cbPal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    } else {
      # The palette with grey:
      cbPal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    }
    
    
  } else if (set == "USV") {
    cat("Returning color palette for", set, fill = T)
    # Color palette corresponding to USVs
    cbPal <- c("Short"="#E69F00", 
               "Flat"="#56B4E9", 
               "Upwards"="#009E73",
               "Downwards"= "#0072B2", 
               "Chevron" = "#BCBD22",
               "Composite"= "#C7C7C7",
               "Harmonics" = "#98DF8A",
               "Complex"= "#D55E00",
               "TwoSyllable" = "#CC79A7",
               "Freq. Steps"= "#9467BD") 
    
    
  } else if (set == "Behav") {
    cat("Returning color palette for", set, fill = T)
    # Color palette corresponding to Behaviors
    cbPal <- c("Play"="#000000", 
               "AG"="#E69F00", 
               "blah"="#56B4E9") 
    
    
  } else if (set == "Opto") {
    cat("Returning color palette for", set, fill = T)
    # Color palette corresponding to Opto Conditions
    cbPal <- c(
      "No Dam" = "grey50",
      "Dam" = "#D55E00",
      # "40 Hz - Sync" = "#CC79A7",
      # "40 Hz - Desync" = "#9467BD")
      # "Sync" = "#CC79A7",
      # "Desync" = "#9467BD")
      "Sync" = "#66C2A5",
      "Desync" = "#FC8D62")
    
  } else if (set == "STFP") {
    cat("Returning color palette for", set, fill = T)
    # Color palette corresponding to STFP Conditions
    cbPal <- c(
      "No Dam" = "grey50",
      "Cinn" = "#D55E00",
      "Coco" = "#56B4E9",
      "Cori" = "#98DF8A")
    
  } else if (set == "Test") {
    cat("Returning color palette for", set, fill = T)
    # Color palette corresponding to STFP Conditions
    cbPal <- c(
      "Demo" = "#CC79A7",
      "NonDemo" = "#56B4E9") 
    
  } else if (set == "Animal") {
    cat("Returning color palette for", set, fill = T)
    # Color palette corresponding to STFP Conditions
    cbPal <- c(
      "Pup" = "#E41A1C",
      "Dam" = "#377EB8")
  }
  
  if (show) {
    show_col(cbPal)
  }
  
  return (cbPal)
}

assign_palette_colors <- function(labels, palette = NULL) {
  library(RColorBrewer)
  labels <- as.character(labels)  # Ensure it's character
  
  # If no palette is given, use RColorBrewer 'Set3'
  if (is.null(palette)) {
    n <- length(labels)
    # Brewer palettes repeat if n > max colors, so handle that
    pal <- brewer.pal(min(n, 12), "Dark2")
    if (n > 12) pal <- rep(pal, length.out = n)
    names(pal) <- NULL
  } else {
    pal <- palette
    n <- length(labels)
    # If named, remove names (will reassign below)
    pal <- unname(pal)
    if (n > length(pal)) pal <- rep(pal, length.out = n)
  }
  # Make named vector for use with scale_*_manual
  out <- setNames(pal[seq_along(labels)], labels)
  return(out)
}


org_Behav_names <- function(df){
  #' Returns dataframe with new cols for organising behav (Dyad, BType, Beh_ID)
  #' 
  #' @param df dataframe: dataframe containing Behaviors column
  #' 
  #' @return dataframe with new cols for organising behav (Dyad, BType, Beh_ID)
  #' Dyad (Dyad, Dam, Pup): Dyadic, Dam, or Pup Behavior
  #' BType (General, Parts, other): General orgnaisation of behaviours
  #' Beh_ID (1-24): Index behaviors
  #' Variable: Behavior names without "Dam" or "Pup"
  
  # Add Label
  lookup_table <- tibble::tribble(
    ~Long_Behav, ~Behavior,
    "Close Proximity", "Cls_Prox",
    "Dam Approach", "Dam_Appr",
    "Dam Investigates Any Body Area", "Dam_Inv_Any",
    "Pup Chase", "Dam_Avoid",
    "Pup Approach", "Pup_Appr",
    "Pup Investigates Any Body Area", "Pup_Inv_Any",
    "Dam Chase", "Pup_Avoid",
    "Window Nose-to-Nose Investigation", "Win_N2N",
    "Nose-to-Nose Investigation", "N2N",
    "Dam Investigates Head", "Dam_Inv_Head",
    "Dam Investigates Back", "Dam_Inv_Body",
    "Dam Investigates Abdominal Area", "Dam_Inv_Abdo",
    "Dam Investigates Anogenital Area", "Dam_Inv_Ano",
    "Dam Investigates Tail", "Dam_Inv_Tail",
    "Pup Investigates Waste", "Pup_Inv_Waste",
    "Pup Investigates Head", "Pup_Inv_Head",
    "Pup Investigates Back", "Pup_Inv_Body",
    "Pup Investigates Abdominal Area", "Pup_Inv_Abdo",
    "Pup Investigates Anogenital Area", "Pup_Inv_Ano",
    "Pup Investigates Tail", "Pup_Inv_Tail",
    "Pup Social Success", "Pup_SocialSuc",
    "Dam Social Success", "Dam_SocialSuc",
    "Pup-Initiated Close Proximity", "Pup_Cls_Prox",
    "Dam-Initiated Close Proximity", "Dam_Cls_Prox",
    "Mutually-Initiated Close Proximity", "Dyad_Cls_Prox")
  
  # lookup_table <- tibble::tribble(
  #   ~Long_Behav, ~Behavior,
  #   "Close Proximity", "Cls_Prox",
  #   "Dam Approach", "Dam_Appr",
  #   "Dam Social Investigation", "Dam_Inv_Any",
  #   "Dam Avoid", "Dam_Avoid",
  #   "Pup Approach", "Pup_Appr",
  #   "Pup Social Investigation", "Pup_Inv_Any",
  #   "Pup Avoid", "Pup_Avoid",
  #   "Window Nose to Nose", "Win_N2N",
  #   "Nose to Nose", "N2N",
  #   "Dam Investigates Head", "Dam_Inv_Head",
  #   "Dam Investigates Body", "Dam_Inv_Body",
  #   "Dam Investigates Abdomen", "Dam_Inv_Abdo",
  #   "Dam Investigates Anogenital", "Dam_Inv_Ano",
  #   "Dam Investigates Tail", "Dam_Inv_Tail",
  #   "Pup Investigates Waste", "Pup_Inv_Waste",
  #   "Pup Investigates Head", "Pup_Inv_Head",
  #   "Pup Investigates Body", "Pup_Inv_Body",
  #   "Pup Investigates Abdomen", "Pup_Inv_Abdo",
  #   "Pup Investigates Anogenital", "Pup_Inv_Ano",
  #   "Pup Investigates Tail", "Pup_Inv_Tail"
  # )
  # 
  # lookup_table <- tibble::tribble(
  #   ~Long_Behav, ~Behavior,
  #   "Dur_Cls_Prox", "Cls_Prox",
  #   "Dur_Dam_Appr", "Dam_Appr",
  #   "Dur_Dam_Inv_Any", "Dam_Inv_Any",
  #   "Dur_Dam_Avoid", "Dam_Avoid",
  #   "Dur_Pup_Appr", "Pup_Appr",
  #   "Dur_Pup_Inv_Any", "Pup_Inv_Any",
  #   "Dur_Pup_Avoid", "Pup_Avoid",
  #   "Dur_Win_N2N", "Win_N2N",
  #   "Dur_N2N", "N2N",
  #   "Dur_Dam_Inv_Head", "Dam_Inv_Head",
  #   "Dur_Dam_Inv_Body", "Dam_Inv_Body",
  #   "Dur_Dam_Inv_Abdo", "Dam_Inv_Abdo",
  #   "Dur_Dam_Inv_Ano", "Dam_Inv_Ano",
  #   "Dur_Dam_Inv_Tail", "Dam_Inv_Tail",
  #   "Dur_Pup_Inv_Waste", "Pup_Inv_Waste",
  #   "Dur_Pup_Inv_Head", "Pup_Inv_Head",
  #   "Dur_Pup_Inv_Body", "Pup_Inv_Body",
  #   "Dur_Pup_Inv_Abdo", "Pup_Inv_Abdo",
  #   "Dur_Pup_Inv_Ano", "Pup_Inv_Ano",
  #   "Dur_Pup_Inv_Tail", "Pup_Inv_Tail"
  # )
  if ("Long_Behav" %in% colnames(df)) {
    lab_df <- df %>% left_join(lookup_table, by = "Long_Behav")
  } else {
    lab_df <- df %>% 
      mutate(Behavior = sub("^Dur_", "", Behavior)) %>% 
      left_join(lookup_table, by = "Behavior")
  }
  
  behav_df <- lab_df %>%
    # left_join(Beh_IDs) %>%
    mutate(Ori = Behavior) %>%
    mutate(Behavior = as.character(gsub("_", " ", Behavior, fixed=TRUE))) %>%
    mutate(Dyad = factor(ifelse(grepl("Dam |Pup ", Behavior),
                                ifelse(grepl("Dam ", Behavior), "Dam", "Pup"), "Dyadic"),
                         levels = c("Dyadic", "Dam", "Pup"))) %>%
    mutate(Behavior = gsub("Dam |Pup ", "", Behavior)) %>% # Remove Dam|Pup str from Behavior
    mutate(BType = case_when(Behavior %in% c("Cls Prox", "Appr", "Chase", "Inv Any", "SocialSuc") ~ "General",
                             Behavior %in% c("Win N2N", "N2N", "Inv Head",
                                             "Inv Body", "Inv Abdo", "Inv Ano", "Inv Tail") ~ "Parts",
                             .default = "Others")) %>%
    mutate(BType = ifelse(Dyad %in% c("Dam","Pup") & Behavior == "Cls Prox", "Others", BType)) %>%
    # Factorize Behaviors
    mutate(BType = factor(BType, levels = c("General", "Parts", "Others")),
           Dyad = factor(Dyad, levels = c("Dyadic", "Dam", "Pup")),
           Behavior = factor(Behavior, levels = c("Cls Prox","Appr", "Chase", "Inv Any", "SocialSuc",
                                                  "Win N2N", "N2N", "Inv Head",
                                                  "Inv Body", "Inv Abdo", "Inv Ano", "Inv Tail",
                                                  "Dyad Cls Prox", "Prox Attmpt","Avoid Lat")))
  
  
  Beh_IDs <- behav_df %>% group_by(BType, Dyad, Behavior) %>% 
    slice_sample(n=1) %>% 
    select(BType, Dyad, Behavior) %>% 
    ungroup() %>% 
    arrange(BType, Dyad, Behavior) %>% 
    mutate(Beh_ID = row_number())
  
  out_df <- behav_df %>% left_join(Beh_IDs) %>% arrange(Beh_ID)
  
  return(out_df)
}



parse_fname <- function(f_list, Phase = "Test", Stage = "Feeding") {
  #' Returns dataframe with filenames split by Family, Phase, date, time
  #' 
  #' @param f_list char vector: list of filenames after regex selection
  #' @param Phase string, NA: filter selected phase only. If NA, no filter. 
  #'                  Options include: "Hab", "Test"
  #' 
  #' @return dataframe with filenames split by Family, Phase, date, time
  if (Stage == "Feeding") {
    fname_df <- f_list %>% data.frame(fpath = ., base = tools::file_path_sans_ext(basename(.))) %>%
      separate(base, into = c("Family", "Phase", "DateThh", "mm", "ss"), sep = "_", remove = F, extra = "drop") %>%
      separate(DateThh, into = c("Date", "hh"), sep = "T") %>%
      mutate(fpath = gsub("/", "\\\\", fpath))
    
    # Filter by phase
    if (is.na(Phase)) {
      cat("Filenames loaded without Phase filtering", fill=T)
    } else {
      cat("Filenames filtered by Phase:", Phase, fill=T)
      fname_df <- fname_df %>% 
        dplyr::filter(., Phase == !!Phase)
    }
  }
  
  if (Stage == "STFP") {
    fname_df <- f_list %>% data.frame(fpath = ., base = tools::file_path_sans_ext(basename(.))) %>%
      separate(base, into = c("Family", "Animal", "Angle", "DateThh", "mm", "ss"), sep = "_", remove = F, extra = "drop") %>%
      separate(DateThh, into = c("Date", "hh"), sep = "T") %>%
      mutate(fpath = gsub("/", "\\\\", fpath))
    
  }
  
  return(fname_df)
}

############################
#### Statistics related ####
############################

run_lmer_contrasts <- function(p_data, response,        # e.g., "tot_calls"
                               group_var,       # e.g., "Group_ID"
                               random_vars = c("Female_ID"),   # default single random effect
                               id_check_df = NULL,             # dataframe for the 'any(n>1)' check
                               id_check_col = "n",             # column in above for duplicate check
                               planned_contrasts = list()
) {
  #' Returns dataframe with filenames split by Family, Phase, date, time
  #' 
  #' @param p_data pd.df: contains columns with variables needed
  #' @param response char, NA: # e.g., "tot_calls"
  #'                  Options include: "Hab", "Test"
  #' 
  #' @return list(model = model, contrast_df = contrast_df)
  #' 
  #' Example usage:
  # results <- run_lmer_contrasts(p_data, "entR", "Group_ID", random_vars = c("Female_ID"), id_check_df = ani_rep, id_check_col = "n",
  #                               planned_contrasts = planned_contrasts)
  # stat_df <- results$contrast_df %>% filter(p.label != "ns")
  
  # NSE for variable names (not strictly necessary, kept for flexibility)
  response_sym   <- rlang::sym(response)
  group_var_sym  <- rlang::sym(group_var)
  
  # determine if Mouse_ID needs to be included (duplicated membership)
  if (!is.null(id_check_df) && any(id_check_df[[id_check_col]] > 1)) {
    cat(
      "Some animals are in multiple conditions. Random variable Mouse_ID added for lmer.\n"
    )
    # make sure Mouse_ID is in random_vars only once
    random_vars <- unique(c(random_vars, "Mouse_ID"))
  }
  
  # Build formula as string and parse
  rand_formula <- paste0("(1 | ", random_vars, ")", collapse = " + ")
  formula_str <- paste(response, "~", group_var, "+", rand_formula)
  model_formula <- as.formula(formula_str)
  
  # Fit model
  model <- lme4::lmer(model_formula, data = p_data)
  
  # Emmeans and pairwise or custom contrasts
  emmeans_obj <- emmeans::emmeans(model, group_var)
  if (length(planned_contrasts) == 0) {
    contrast_obj <- pairs(emmeans_obj, adjust = "holm")
  } else {
    contrast_obj <- contrast(emmeans_obj, method = planned_contrasts, adjust = "holm")
  }
  
  contrast_df <- contrast_obj %>%
    as.data.frame() %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
    dplyr::mutate(
      group1  = stringr::str_trim(stringr::str_replace_all(group1, "^\\(|\\)$", "")),
      group2  = stringr::str_trim(stringr::str_replace_all(group2, "^\\(|\\)$", "")),
      p.label = symnum(
        p.value, corr = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      ),
      y.position = max(p_data[[response]], na.rm = TRUE) * 1.1
    )
  print(contrast_df)
  return(list(model = model, contrast_df = contrast_df))
}




##########################
#### Plotting related ####
##########################

create_facet_design <- function(t_data) {
  # Script to produce manual facet design for use with ggh4x::facet_manual()
  
  # Ensure input data has the correct columns
  if (!all(c("Group_ID", "Mouse_ID") %in% names(t_data))) {
    stop("Data frame must contain 'Group_ID' and 'Mouse_ID' columns.")
  }
  # Count number of Mouse_ID per Group_ID, and preserve original order
  group_counts <- t_data %>%
    group_by(Group_ID) %>%
    summarise(count = n(), .groups = "drop")
  
  # Determine max number of columns (Mouse_IDs in largest group)
  max_cols <- max(group_counts$count)
  
  # Prepare letters (repeat if more than 26 are needed)
  total_letters <- sum(group_counts$count)
  alphabet <- toupper(letters)  # A-Z
  letter_sequence <- rep(alphabet, length.out = total_letters)
  
  letter_index <- 1
  output_lines <- character()
  
  for (i in seq_len(nrow(group_counts))) {
    n_mice <- group_counts$count[i]
    # Get the subset of letters
    current_letters <- letter_sequence[letter_index:(letter_index + n_mice - 1)]
    letter_index <- letter_index + n_mice
    
    # Pad with '#' if needed
    if (n_mice < max_cols) {
      current_letters <- c(current_letters, rep("#", max_cols - n_mice))
    }
    
    # Form a row string
    row_str <- paste(current_letters, collapse = "")
    output_lines <- c(output_lines, row_str)
  }
  
  # Return the formatted design as single string
  return(paste(output_lines, collapse = "\n"))
}

create_facet_design_colwise <- function(t_data) {
  # Script to produce manual facet design for use with ggh4x::facet_manual()
  
  # Ensure input data has the correct columns
  if (!all(c("Group_ID", "Mouse_ID") %in% names(t_data))) {
    stop("Data frame must contain 'Group_ID' and 'Mouse_ID' columns.")
  }
  
  # Count how many Mouse_IDs per Group
  group_counts <- t_data %>%
    dplyr::group_by(Group_ID) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  # Sort group_counts in desired columnwise order
  group_counts <- group_counts[order(group_counts$Group_ID), ]  # Optional reordering (alphabetical)
  
  # Total number of individual Mouse_IDs (i.e., cells to fill with letters)
  total_cells <- sum(group_counts$count)
  
  # Determine total number of columns needed: fill columns top-down
  # max group count is the height of each column
  n_rows <- max(group_counts$count)
  n_cols <- ceiling(total_cells / n_rows)
  
  # Generate letter sequence
  alphabet <- toupper(letters)
  letter_sequence <- rep(alphabet, length.out = total_cells)
  
  # Build a vector that represents the column-wise stacking of group IDs
  group_vector <- unlist(apply(group_counts, 1, function(row) rep(row[["Group_ID"]], row[["count"]])))
  
  # Fill matrix column-first with letters and group info
  mat <- matrix("#", nrow = n_rows, ncol = n_cols)
  letter_index <- 1
  
  for (col in 1:n_cols) {
    for (row in 1:n_rows) {
      if (letter_index <= total_cells) {
        mat[row, col] <- letter_sequence[letter_index]
        letter_index <- letter_index + 1
      }
    }
  }
  
  # Transpose for row-wise printing
  row_strings <- apply(mat, 1, paste, collapse = "")
  
  return(paste(row_strings, collapse = "\n"))
}
