#' -----------------------------------------------------------------------------
#' Extract USV data from DeepSqueak exported excels
#' by KE - Jul 2023
#' modified by Gao - Jan-May 2025
#'
#' -----------------------------------------------------------------------------
#' Functions:
#'   1. 
#'   2. 
#'   3. 
#'   4. 
#'
#' -----------------------------------------------------------------------------
#' Pre-requisites:
#'   1. Exported DeepSqueak detection excels
#'   2. Excel files named in format: 
#'          (Expt)_(Group_ID)_(Mouse_ID)_(Female_ID)_(Date)_(Time)
#'
#' -----------------------------------------------------------------------------

getUSVdata <- function(experiment = "none", save = T, data_dir = NA, fname_org = list()) {
  #' Returns compiled dataframe with all excel logs in chosen data_dir.
  #' 
  #' @param experiment str, "none": Returns reminder to manually input fname_org
  #' @param save bool, True: Saves final data frame as .csv file.
  #' @param data_drive str, NA: Returns directory to choose folder with data. 
  #' @param fname_org list(str), NA: Input list of fname org divided by "_" or spaces
  #'            e.g.  c("Expt","Group_ID", "Mouse_ID","Female_ID", "Date","Time","Other")
  #' 
  #' @description
  #' Create compiled dataframe with filename as new column
  #' Create new columns - Group_ID, Mouse_ID, PD, Phase
  #' Select variables of interest (Group, Mouse, Age, Phase, Call Type, 
  #' start and End time, Call duration, Frequency (kHz) and Mean power)
  #' Remove Calls classified as Noise in DeepSqueak but accidentally accepted 
  #'      during manual pressing of keys
  #' @return list of two items
  #'            fname_df: broken down list of file names
  #'            raw_df: final data frame with selected variables (excludes expt with no calls)
  {
    #### LOAD PACKAGES -------------------------------------------------------------
    library(dplyr)
    library(readxl)
    library(tidyr)
    library(tidyverse)
    library(fs)
    source("Utils.R")
  }
  
  {
    # Presets for different data
    if (experiment == "none" && is_empty(fname_org)) return (stop("Input list for fname_org"))
    if (experiment == "FI") fname_org <- c("Group_ID", "PD", "FI_Phase", "Mouse_ID", "DateTime")
    if (experiment == "STFP") fname_org <- c("Group_ID", "Mouse_ID", "DateTime")
    if (experiment == "Hensch") fname_org <- c("Expt","Group_ID", "Mouse_ID","Female_ID", "Date","Time","Other")
  }
  
  #' Load data 
  #' Set Directories
  #' Set data_dir  to folder in BabyLINC Animal Backup containing excels
  if (is.na(data_dir)) {data_dir <- choose.file(caption = "Select folder containing DeepSqueak exported Excel")} 
  cat("Data are taken from: ", data_dir, fill = T)
  
  #' Create list of relevant file names
  USV_files <- fs::dir_ls(data_dir, regexp = "\\.xlsx$")
  
  fname_df <- data.frame(filename = tools::file_path_sans_ext(basename(USV_files))) %>% 
    separate(filename, into = fname_org, sep = "_| ", extra = "merge") %>% 
    mutate(Group_ID = gsub("ch\\d", "", Group_ID)) %>% 
    select(-Other)
  
  raw_df <-
    USV_files %>%
    map_dfr(read_excel,.id = "filename") %>% 
    mutate(filename = tools::file_path_sans_ext(basename(filename))) %>% 
    separate(filename, into = fname_org, sep = "_| ", extra = "merge") %>% 
    mutate(Group_ID = gsub("ch\\d", "", Group_ID)) %>% 
    # filter(Accepted == T) %>% 
    filter(Accepted == T, !Label == "USV") %>% # Manually excluded missed calls (TEMPORARY SOlUTION)
    # select(Expt, Group_ID, Mouse_ID, Female_ID, Label, ID, `Begin Time (s)`, `End Time (s)`,
    #        `Call Length (s)`, `Principal Frequency (kHz)`, `Mean Power (dB/Hz)`, Date, Time) %>% 
    subset(Label != "Noise") # Remove noise category
  
  output_df <- raw_df %>% distinct(Expt, Group_ID, Mouse_ID) # Check for files with no calls
  
  files_w_no_calls <- setdiff(select(fname_df, c("Expt", "Group_ID", "Mouse_ID")), output_df)
  if (nrow(files_w_no_calls)>1) {
    cat("These files are excluded as they have no calls:", fill=T)
    print(files_w_no_calls)
  }
  
  if (save==T) {
    main_dir <- getwd()
    out_dir <- file.path(main_dir, "USV_summary")
    dir.create(out_dir, showWarnings = FALSE)
    
    output_name <- file.path("USV_summary", paste0(experiment,"_summary.csv"))
    write.csv(raw_df, output_name, row.names = FALSE)
    cat("New dataframes saved in", out_dir, fill = T)
  }
  return (list(fname_df = fname_df, raw_df = raw_df))
}


check_missed_USVs <- function(data_dir, call_order, verbose = T) {
  #' Returns dataframe with basename, number of missed calls with 'USV' label & Call Index for each 'USV' label
  #' 
  #' @param data_drive str: Select directory to choose folder with data extracted excel logs from DeepSqueak.
  #' @param call_order list: List of Syllables to organise call order
  #' @param verbose boolean, T: Prints file name and number of calls in each category.
  #' 
  #' @return dataframe with basename, USV_n, Call_IDs
  
  library(readxl)
  # List all xlsx files in the data directory
  file_list <- list.files(path = data_dir, pattern = "*.xlsx", full.names = TRUE)
  
  # List to store all data frames
  data_frames <- list()
  for (file in file_list) {
    # Read the data from each xlsx file into a data frame
    data <- read_excel(file)
    
    # Store the data frame in the list
    data_frames[[file]] <- data %>% 
      filter(Accepted == T)
    # filter(Accepted == T, Label != "USV")
  }
  
  missed_calls <- data.frame()
  for (i in seq_along(names(data_frames))) {
    df <- data_frames[[i]]
    if (verbose) {
      call_n <- df %>% group_by(Label) %>% summarize(n= n())
      cat(names(data_frames)[i])
      print(call_n)
    }
    
    # Check for extra call categories
    extracalls <- df %>% 
      mutate(extra = ifelse(!Label %in% call_order, T, F)) %>% 
      filter(extra == T) %>% 
      group_by(Label) %>% summarize(n= n())
    
    # Check for duplicate calls
    duplicatecalls <- df %>%
      group_by(across(-ID)) %>%         # Group by everything except 'Value2'
      filter(n_distinct(ID) > 1) %>%    # Keep only rows where Value2 varies
      ungroup() %>% select(ID)
    
    if (dplyr::count(extracalls)>0) {
      missed <- dplyr::filter(df, !Label  %in% call_order) %>% select(ID)
      missed_calls <- extracalls %>%
        mutate(basename = basename(names(data_frames)[i]), .before = Label,
               'missed_Call_IDs' = toString(missed), 
               'duplicate_Call_IDs' = toString(duplicatecalls)) %>%
        rbind(missed_calls,.) 
    }
    
    
    if (dplyr::count(extracalls)==0 & length(duplicatecalls>0)) {
      missed <- dplyr::filter(df, !Label  %in% call_order) %>% select(ID)
      missed_calls <- data.frame(basename = basename(names(data_frames)[i]),
                                 'missed_Call_IDs' = NA, 
                                 'duplicate_Call_IDs' = toString(duplicatecalls),
                                 Label = NA, n = 0) %>%
        rbind(missed_calls,.) 
    }
    
  }
  
  if (any(!is.na(missed_calls$missed_Call_IDs)) | length(duplicatecalls>0)) {
    cat('Unlabelled calls found in', fill=T)
    if (any(!is.na(missed_calls$missed_Call_IDs))){
      missed_calls <- missed_calls %>% pivot_wider(names_from = Label, values_from = n) %>% 
        relocate(missed_Call_IDs, .after = last_col()) %>% 
        mutate(duplicate_Call_IDs = ifelse(duplicate_Call_IDs == "numeric(0)", NA, duplicate_Call_IDs))
    }
    print(missed_calls)
  } else {
    cat("No missed calls or duplicate calls.")
  }
  
  # Save missed calls
  if (length(missed_calls)>0) {
    expt <- missed_calls$basename[1] %>% str_extract("^[A-Z0-9]+")
    output_name <- file.path(plot_dir, paste0(expt,"missed_calls.csv"))
    write.csv(missed_calls, output_name, row.names = FALSE)
    cat("missed calls saved in", plot_dir, fill = T)
  }
  
  return (missed_calls)
} 

# getUSVdata("STFP", data_dir = "H:/Raw_Data/4_STFP/3_USV_Analysis/4_ExportedUSVs/DyPV6-16")

getVidCoding <- function(VCode_file = NA, save = T) {
  #' Returns dataframe with video coding of specific experiments.
  #' 
  #' @param save bool, True: Saves final data frame as .csv file.
  #' @param VCode_file str, NA: Returns directory to choose folder with data. 
  #' 
  #' @description
  #' Create compiled dataframe with filename as new column
  #' Create new columns - Group_ID, Mouse_ID, PD, Phase
  #' Select variables of interest (Group, Mouse, Age, Phase, Call Type, 
  #' start and End time, Call duration, Frequency (kHz) and Mean power)
  #' Remove Calls classified as Noise in DeepSqueak but accidentally accepted 
  #'      during manual pressing of keys
  #' @return list of two items
  #'            behav_df: broken down list of file names
  #'            VCode_fpath: file path of Video coding excel
  #'            dropped_ani: final data frame with selected variables (excludes expt with no calls)
  {
    #### LOAD PACKAGES -------------------------------------------------------------
    library(dplyr)
    library(readxl)
    library(tidyr)
    library(tidyverse)
    library(fs)
    library(hms)
    source("utils.R")
  }
  
  #' Load data 
  #' Set Directories
  #' Set data_dir  to folder in BabyLINC Animal Backup containing excels
  if (is.na(VCode_file)) {VCode_file <- choose.files(caption = "Select file corresponding to Video Coding")} 
  cat("Video Coding data taken from: ", VCode_file, fill = T)
  
  raw_df <-
    VCode_file %>% read_excel() %>% 
    mutate(Treatment = ifelse(Treatment == "Ctrl", "", Treatment)) %>% 
    mutate(Group_ID = paste0(Group_ID, Treatment)) %>% 
    mutate(Aud_Beep = as.numeric(parse_hms(paste0("00:",`Audio Time Beep Start`))), .keep = "unused") %>% 
    mutate(Vid_Beep = as.numeric(parse_hms(paste0("00:",`Video Time Light (Start)`))), .keep = "unused") %>% 
    mutate(Vid_FEntry = as.numeric(parse_hms(paste0("00:",`Video Time Female Added (when all 4 paws touch floor)`))), .keep = "unused") %>% 
    mutate(Vid_Mating = as.numeric(parse_hms(paste0("00:",`Successful Mating?`)))) %>% 
    mutate(`Successful Mating?` = ifelse(is.na(Vid_Mating), `Successful Mating?`, "Yes")) %>% 
    mutate(FEntry = Vid_FEntry - Vid_Beep) %>% 
    mutate(Female_Entry = FEntry + Aud_Beep)
  
  output_df <- raw_df %>% drop_na(-Notes, -Vid_Mating)
  ani_w_no_Alignment <- setdiff(raw_df, output_df) %>% select(-Vid_Mating) # Check for files with no calls
  
  
  if (nrow(ani_w_no_Alignment)>1) {
    cat("Dropped", nrow(ani_w_no_Alignment), "animals due to lack of data.")
    print(ani_w_no_Alignment)
  }
  
  # if (save==T) {
  #   main_dir <- getwd()
  #   out_dir <- file.path(main_dir, "USV_summary")
  #   dir.create(out_dir, showWarnings = FALSE)
  #   
  #   output_name <- file.path("USV_summary", paste0(experiment,"_summary.csv"))
  #   write.csv(raw_df, output_name, row.names = FALSE)
  #   cat("New dataframes saved in", out_dir, fill = T)
  # }
  return (list(behav_df = output_df, VCode_fpath = VCode_file, dropped_ani = ani_w_no_Alignment))
}


getDSSyntaxdata <- function(experiment = "none", data_dir = NA) {
  #' Returns compiled nested list of dataframes extracted from DeepSqueak Syntax excel
  #' 
  #' @param experiment str, "none": Returns reminder to input type of experiment (STFP or FI).
  #' @param data_drive str, NA: Returns directory to choose folder with data. 
  #' 
  #' @description
  #' Extract dataframes from DeepSqueal Syntax excel
  #' Sort and fill Syllables based on call_order
  #' Create nested list of dataframes, file_data_list[[filename_base]]:
  #' - File_Breakdown, 
  #' - Conditional_Probability
  #' - Transition_Count
  #' - Total_Count
  
  #' @return final data frame with selected variables
  {
    #### LOAD PACKAGES -------------------------------------------------------------
    library(dplyr)
    library(readxl)
    library(tidyr)
    library(tidyverse)
    library(fs)
    source("utils.R")
  }
  
  if (experiment == "none") return (cat("Set experiment as FI or STFP."))
  if (experiment == "Hensch") {
    fname_org <- c("Expt","Group_ID", "Mouse_ID","Female_ID", "DateTime")
    call_order <- c(
      'Flat', 
      'Downward', 
      'Upward', 
      'Composite',
      'Chevron', 
      'Reverse Chevron',
      '1 Frequency Step', 
      '2 Frequency Step', 
      '3+ Frequency Step',
      'Complex', 
      'Harmonics', 
      'Short', 
      'Unstructured',
      'Murmur',
      'Noisy'
    )
  }
  
  if (is.na(data_dir)) {data_dir <- choose.dir(caption = "Select folder containing DeepSqueak exported Excel")} 
  cat("Data are taken from: ", data_dir, fill = T)
  
  #' Create list of relevant file names
  USV_files <- fs::dir_ls(data_dir, regexp = "\\.xlsx$")
  
  # Initialize a list to store combined dataframes
  file_data_list <- list()
  
  # Function to extract section indices
  get_section_indices <- function(data, section_title) {
    return(grep(section_title, data[[1]], ignore.case = TRUE))
  }
  
  # Function to report mismatches
  report_mismatch <- function(existing, call_order, filename, section) {
    existing_categories <- unique(existing$Category)
    missing_categories <- setdiff(call_order, existing_categories)
    extra_categories <- setdiff(existing_categories, call_order)
    
    if (length(missing_categories) > 0 || length(extra_categories) > 0) {
      if (length(extra_categories) > 0) {
        cat("Mismatch in", section, "of", filename, ":\n")
        cat("  Extra:", paste(extra_categories, collapse=", "), "\n")
        if (length(missing_categories) > 0) {
          cat("  Missing:", paste(missing_categories, collapse=", "), "\n")
        }
      }
    }
  }
  
  
  
  # Loop over each file
  for (file_path in USV_files) {
    
    # Breakdown file name
    fname_break <-
      data.frame(filename = tools::file_path_sans_ext(basename(file_path))) %>% 
      separate(filename, into = fname_org, sep = "_| ", extra = "merge")
    
    # Load the Excel file
    data <- suppressMessages(read_excel(file_path, col_names = FALSE, trim_ws = TRUE))
    filename_base <- tools::file_path_sans_ext(basename(file_path))
    
    # Determine the start indices of each section
    cond_prob_index <- get_section_indices(data, "Conditional Probability")
    trans_count_index <- get_section_indices(data, "Transition Count")
    total_count_index <- get_section_indices(data, "Total Count")
    
    # Extract the header
    header_index <- cond_prob_index + 1
    headers <- as.character(unlist(data[header_index,]))
    
    # Sort headers to match the call_order, use intersection to ignore extras
    sorted_headers <- intersect(call_order, headers)
    
    # Extract dataframes based on section indices
    conditional_probability <- data[(header_index + 1):(trans_count_index - 4), ]
    transition_count <- data[(trans_count_index + 2):(total_count_index - 4), ]
    total_count <- data[(total_count_index + 1):nrow(data), 1:2]
    
    # Assign column names
    colnames(conditional_probability) <- headers
    colnames(transition_count) <- headers
    colnames(total_count) <- c("Category", "Total Count")
    
    report_mismatch(conditional_probability, call_order, filename_base, "Conditional Probability")
    
    # Apply the function to both dataframes
    conditional_probability <- adjust_dataframe(conditional_probability, call_order)
    transition_count <- adjust_dataframe(transition_count, call_order)
    
    # Total_Count only needs rows filled, not columns
    total_count <- total_count %>%
      mutate(`Total Count` = as.numeric(`Total Count`)) %>% 
      complete(Category = factor(call_order, levels = call_order), fill = list(`Total Count` = 0)) %>%
      arrange(match(Category, call_order))
    
    # Store dataframes as a list within the list
    file_data_list[[filename_base]] <- list(
      File_Breakdown = fname_break,
      Conditional_Probability = conditional_probability,
      Transition_Count = transition_count,
      Total_Count = total_count
    )
  }
  return (file_data_list)
}


getSyntaxdata <- function(USV_excels, fname_org, call_order) {
  #' Returns compiled nested list of Syntax related dataframes extracted from DeepSqueak excel logs
  #' 
  #' @param file_list list: List of file locations to analyse (DeepSqueak USV export logs)
  #' @param fname_org list: List of categories to parition file name
  #' @param call_order list: List of Syllables to organise call order
  #' 
  #' @description
  #' Extract dataframes from DeepSqueal Syntax excel
  #' Sort and fill Syllables based on call_order
  #' Create nested list of dataframes, file_data_list[[filename_base]]:
  #' - File_Breakdown, 
  #' - Conditional_Probability
  #' - Transition_Count
  #' - Total_Count
  
  #' @return final data frame with selected variables
  
  # Initialize a list to store combined dataframes
  file_data_list <- list()
  
  # Loop over each file
  for (i in 1:length(USV_excels)) {
    filename_base <- tools::file_path_sans_ext(basename(names(USV_excels)[i]))
    # Breakdown file name
    fname_break <-
      data.frame(filename = filename_base) %>% 
      separate(filename, into = fname_org, sep = "_| ", extra = "merge")
    
    # Check if accepted calls exists
    if (nrow(USV_excels[[i]])> 0) {
      # Generate Transition Counts df
      transition_counts <- USV_excels[[i]] %>%
        count(Category, next_call_Category) %>%
        pivot_wider(names_from = next_call_Category, values_from = n, values_fill = 0) %>% 
        adjust_dataframe(call_order)
      
      # Generate Conditional Probability df
      transition_matrix <- USV_excels[[i]] %>%
        count(Category, next_call_Category) %>%
        group_by(Category) %>% 
        mutate(transition_probability = n / sum(n)) %>%  # Normalize counts to probabilities
        select(Category, next_call_Category, transition_probability) %>%
        pivot_wider(names_from = next_call_Category, values_from = transition_probability, values_fill = 0) %>% 
        ungroup() %>% 
        adjust_dataframe(call_order)
    } else {
      # If no calls, create dummy dataframe
      transition_counts <- data.frame(Category = call_order, next_call_Category=call_order, n=0) %>% 
        pivot_wider(names_from = next_call_Category, values_from = n, values_fill = 0) %>% 
        adjust_dataframe(call_order)
      transition_matrix <- transition_counts
    }
    
    # Generate Total Counts table
    {
      total_count <- USV_excels[[i]] %>% 
        group_by(Category) %>% 
        summarize('Total Count' = n()) %>% 
        complete(Category = factor(call_order, levels = call_order), fill = list(`Total Count` = 0))
    }
    
    # Store dataframes as a list within the list
    file_data_list[[filename_base]] <- list(
      File_Breakdown = fname_break,
      Conditional_Probability = transition_matrix,
      Transition_Count = transition_counts,
      Total_Count = total_count
    )
  }
  return (file_data_list)
} 


adjust_dataframe <- function(df, call_order) {
  #' Returns df with column and rows filled and sorted with labels from call_order. Missing labels are filled with 0.
  #' 
  #' @param df, dataframe : Input dataframe with Category
  #' @param call_order list: Label order in desired sequence
  #' 
  #' @return filled and sorted data frame of calls
  
  # Adjust headers and fill missing columns with zero
  df_long <- df %>%
    pivot_longer(cols = -Category, names_to = "Call_Type", values_to = "Value") %>%
    mutate(Value = as.numeric(Value)) %>% 
    complete(Category = factor(call_order, levels = call_order), Call_Type = factor(call_order, levels = call_order), fill = list(Value = 0)) %>%
    pivot_wider(names_from = Call_Type, values_from = Value)
  
  return(df_long)
}


getTransitionTables <- function(main_data, seq_order, min_calls_in_seq, min_IBI) {
  #' Returns df transition tables generated for individual call sequences. Sequences < 3 calls are filtered out.
  #' 
  #' @param main_data, dataframe : Input dataframe with Group_ID, Mouse_ID, Label, Begin Time (s), End Time (s), Call Length (s), ISI
  #' @param seq_order, list : list with c(IBI, call_order)
  #' @param min_calls_in_seq, int : value indicating minimum calls in bout to be considered a sequence
  #' 
  #' @return transition_tables : dataframe with Group_ID, Mouse_ID, Count_Matrix, Prob_Matrix
  
  # Generate dataframe of IBI rows
  IBI_df <- main_data %>% 
    dplyr::mutate(prev_endtime = stats::lag(`End Time (s)`)) %>% 
    filter(!is.na(ISI),ISI>min_IBI) %>% 
    mutate(
      `Begin Time (s)` = prev_endtime,
      `End Time (s)` = `Begin Time (s)` + ISI,
      `Call Length (s)` = ISI,
      Label = "IBI"
    ) %>% 
    select(-prev_endtime)
  cat("Transition tables generated with ISI cutoff:", min_IBI, "s", fill = T)
  
  # Add Inter bout interval rows to main_data
  seq_df <- main_data %>% 
    bind_rows(IBI_df) %>% 
    group_by(Group_ID, Mouse_ID) %>% 
    arrange(`Begin Time (s)`) %>% 
    filter(!calls_in_bout<min_calls_in_seq) %>%
    mutate(
      seq_ID = row_number(),
      Next.Label = lead(Label),.before = ID
    ) 
  b <- seq_df %>% filter(Mouse_ID == "C57M8035", Group_ID == "BSL")
  #### Fix Filtering for number of calls in seq ####
  
  
  # Add IBI to start and end of each sequence
  label_seqs <- seq_df %>%
    group_by(Group_ID, Mouse_ID, Female_ID) %>%
    summarise(
      Label = list(c("IBI", as.character(Label), "IBI")),
      .groups = "drop"
    ) %>% unnest(cols = c(Label)) %>% 
    group_by(Group_ID, Mouse_ID, Female_ID) %>%
    mutate(Next.Label = lead(Label)) %>% 
    filter(!(Label == Next.Label & Label == "IBI")) %>% # Filter repeating IBIs due to filtering
    mutate(seq_ID = row_number())
  
  transitions <- label_seqs %>% 
    filter(!is.na(Next.Label)) %>% 
    dplyr::count(Label, Next.Label, name = "Count") %>%  
    mutate(
      Label = factor(Label, levels = seq_order),
      Next.Label = factor(Next.Label, levels = seq_order)
    ) %>%
    complete(Label, Next.Label, fill = list(Count = 0))
  
  
  compute_transitions <- function(transitions, keys) {
    count_matrix <- transitions %>%
      pivot_wider(
        names_from = Next.Label,
        values_from = Count,
        values_fill = 0
      ) %>%
      column_to_rownames("Label") %>%
      as.matrix()
    
    prob_matrix <- prop.table(count_matrix, margin = 1)
    prob_matrix[is.nan(prob_matrix)] <- 0
    
    tibble(
      Group_ID = keys$Group_ID,
      Mouse_ID = keys$Mouse_ID,
      Female_ID = keys$Female_ID,
      Count_Matrix = list(count_matrix),
      Prob_Matrix = list(prob_matrix)
    )
  }
  transition_tables <- transitions %>%
    group_map(compute_transitions) %>%
    bind_rows()
  
  return (transition_tables)
}


calMarkovEnt <- function(transition_tables) {
  #' Returns df with column and rows filled and sorted with labels from call_order. Missing labels are filled with 0.
  #' 
  #' @param transition_tables, dataframe : Input dataframe with Group_ID, Mouse_ID, Count_Matrix, Prob_Matrix
  #' 
  #' @return list(sm_all, entrate_all): list with df of Stationary Matrix and df Markov Chain Entropy values
  
  
  ## Markov Chain Analysis ---------------------------------------------------
  #' Reference: 
  #' Adapted from https://github.com/bvegetabile/ccber/tree/master.
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library('ccber') # For markov chain entropy cal.
  
  # Initialize result data frames
  sm_all <- data.frame()
  entrate_all <- data.frame()
  
  # Loop over each group/subject
  for (i in seq_len(nrow(transition_tables))) {
    cat("Calculating Markov Entropy for:",
        paste(transition_tables$Group_ID[i], transition_tables$Mouse_ID[i]),
        fill = TRUE)
    
    # Get transition probability matrix
    prob_matrix <- transition_tables$Prob_Matrix[[i]]
    
    # Category labels
    categories <- colnames(prob_matrix)
    
    # -----------------------------
    # Calculate Stationary Probabilities
    # -----------------------------
    # Preferred: Empirical (e.g., based on steady-state sequence), but assuming eigen fallback:
    sm <- CalcEigenStationary(prob_matrix)
    
    # Convert to data.frame
    sm_df <- data.frame(Category = categories, sm) %>%
      pivot_wider(names_from = Category, values_from = sm)
    
    # Add identifying info
    sm_df <- sm_df %>%
      mutate(
        Group_ID = transition_tables$Group_ID[i],
        Mouse_ID = transition_tables$Mouse_ID[i],
        Female_ID = transition_tables$Female_ID[i]
      )
    
    sm_all <- bind_rows(sm_all, sm_df)
    
    # -----------------------------
    # Calculate Markov Entropy Rate
    # -----------------------------
    entrate <- CalcMarkovEntropyRate(prob_matrix, sm)
    
    entrate_df <- tibble(
      Group_ID = transition_tables$Group_ID[i],
      Mouse_ID = transition_tables$Mouse_ID[i],
      Female_ID = transition_tables$Female_ID[i],
      entR = entrate
    )
    
    entrate_all <- bind_rows(entrate_all, entrate_df)
    
  }
  
  return(list(sm_all = sm_all, entrate_all = entrate_all))
}
