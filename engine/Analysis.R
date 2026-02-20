#' USV bout Properties Analysis script (Fragmented Care) --------------------------------------------------------------------
#' June 2025 -------------------------------------------------------------------
#' by Gao
#' modified by Jason Lin January 2026
#' 
#' #' -----------------------------------------------------------------------------
#' Pre-requisites:
#'   1.
#' 
#'Overview
#'1. 
#'


## 1.1 Clear Environment -------------------------------------------------------
# can cause errors, comment out if needed
# rm(list = ls())

## 1.2 Load Required Packages --------------------------------------------------
{
  library(dplyr)
  library(stringr)
  library(forcats)
  library("colorspace")
  library(tidyr)
  source("Utils.R")
  source("extractUSV.R")

  
  library(ggplot2)
  library(scales)
  library(ggpubr)
  library(ggdist) # distribution plots - halfeye
  library(lmerTest)
  library(broom.mixed) # 
  library(emmeans) # For multiple pairwise comparisons
  
}

## 1.3 Set Output Directories --------------------------------------------------
{
  # should be passed in by Shiny
  if (!exists("main_dir")) main_dir <- getwd()
  
  if (!exists("data_dir")) {
    data_dir <- dirname(file.choose())
  }
  
  if (!exists("plot_dir")) {
    plot_dir <- file.path(dirname(data_dir), paste0("Plots_", Sys.Date()))
  }
  
  if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  cat("Current main_dir: ", main_dir, fill = TRUE)
  cat("Current data_dir: ", data_dir, fill = TRUE)
  cat("Plots are saved in: ", plot_dir, fill = TRUE)
}

## 1.4 Set User defined variables ----------------------------------------------
{
  fname_org = c("Expt","Group_ID", "Mouse_ID","Female_ID", "Date","Time","Other")
  min_calls <- 200
  min_IGI <- 0.125 # minimum inter group interval in seconds
  min_IBI <- 0.225 # minimum inter bout interval in seconds
  call_order <- c(
    'Short', 
    'Upward', 
    'Downward', 
    '1 Frequency Step', 
    '2 Frequency Step', 
    '3+ Frequency Step',
    'Flat', 
    'Harmonics', 
    'Complex', 
    'Chevron', 
    'Reverse Chevron',
    'Composite',
    'Unstructured',
    'Murmur',
    'Noisy'
  )
  call_pal <- assign_palette_colors(call_order, hue_pal()(length(call_order)))
  call_pal_dark <- assign_palette_colors(call_order, hue_pal(l = 50, c = 100)(length(call_order)))
  show_col(call_pal)
  
  # list of orders for each expt condition.
  cond_order_list <- list(
    Cerevel = c("BSL", "Veh-30Min", "0.32-30Min", "Veh-4Hr", "0.32-4Hr"),
    # Cerevel = c("BSL", "Veh-30Min", "Veh-4Hr", "0.32-30Min",  "0.32-4Hr"),
    NR2A = c("WT", "Het", "KO"),
    C57 = c("CC", "FC", "FCNAC")
    # C57 = c("CC", "FC")
  )
  
  group_pal <- list(
    Cerevel = c("BSL" = "#1B9E77", 
                "Veh-30Min" = "#D95F22", "0.32-30Min" = "#7570C3", 
                "Veh-4Hr" = "#D95F02", "0.32-4Hr" = "#7587D3"),
    NR2A = c("WT" = "#B3DE69", "Het" = "#FCCDE5", "KO" = "#D9D9D9"),
    C57 = c("CC" = "#FB8072", "FC" = "#80B1D3", "FCNAC" = "#FDB462")
    # C57 = c("CC" = "#FB8072", "FC" = "#80B1D3")
  )
  
  group_design <- list(
    Cerevel = c("ABC\n#DE"),
    NR2A = c("ABC"),
    C57 = c("A\nB\nC")
    # C57 = c("A\nB")
  )
  
  # Define label mapping (like a labeller, but for x-axis)
  group_labels_list <- list(
    Cerevel = c(
      "BSL" = "Baseline",
      "Veh-30Min" = "Veh \n(30 min)",
      "0.32-30Min" = "0.32 \n(30 min)",
      "Veh-4Hr" = "Veh \n(4 hr)",
      "0.32-4Hr" = "0.32 \n(4 hr)"
    ),
    NR2A = c("WT" = "WT", "Het" = "Het", "KO" = "KO"),
    C57 = c("CC" = "Control\nCare (CC)", "FC" = "Frag.\nCare (FC)", "FCNAC" = "FCNAC Treated\nFC")
    # C57 = c("CC" = "Control\nCare (CC)", "FC" = "Frag.\nCare (FC)")
  )
  
  ## 1.5 Statistical Comparisons -------------------------------------------------
  # planned_contrasts <- list() # Set to empty list if intending to do all contrast levels
  group_contrasts <- list(
    Cerevel = list(
      "Veh-30Min - 0.32-30Min" = c("BSL" = 0, "Veh-30Min" = -1, "0.32-30Min" = 1, "Veh-4Hr" = 0, "0.32-4Hr" = 0),
      "Veh-4Hr - 0.32-4Hr"     = c("BSL" = 0, "Veh-30Min" = 0, "0.32-30Min" =0, "Veh-4Hr" = -1, "0.32-4Hr" = 1)
    ),
    NR2A = list(),
    C57 = list()
  )
}

# Check for missed calls (Labelled as USV) 
missed_calls <- check_missed_USVs(data_dir, call_order, F)


## 1.6 Load and Prepare Dataset ------------------------------------------------
# Loads all excel logs in drive into a single dataframe
{
  a <- getUSVdata(data_dir = data_dir, save = F, fname_org = fname_org) 
  fname_df <- a$fname_df
  raw_df <- a$raw_df %>% 
    group_by(Group_ID, Mouse_ID) %>% 
    mutate(Category = factor(Label, levels = call_order)) %>% 
    arrange(`Begin Time (s)`) %>% 
    mutate(next_call_time = c(diff(`Begin Time (s)`),NA),
           next_call_Category = lead(Category),
           call_id = row_number(),
           same_bout = next_call_time <= min_IBI,) %>% 
    mutate(bout_id = cumsum(!lag(same_bout, default=TRUE))+1) %>% 
    mutate(ISI = `Begin Time (s)` - lag(`End Time (s)`) # Not the same as next call time. This is the period of silence between calls
    ) %>% 
    add_count(bout_id, name = "calls_in_bout") %>% 
    group_by(Group_ID, Mouse_ID, bout_id) %>% 
    mutate(call_seq = row_number()) %>% 
    filter(!is.na(next_call_time)) %>% 
    ungroup()
}

# Add factor levels to dataframe
{
  raw_df$Group_ID %>% unique() # Display conditions list
  expt <- raw_df$Expt[1]
  cond_order <- cond_order_list[[expt]]
  group_labels <- group_labels_list[[expt]]
  planned_contrasts <- group_contrasts[[expt]]
  sel_pal <- group_pal[[expt]]
  n_cond <- length(cond_order)
  
  # Sort by most calls
  t_data <- raw_df %>% 
    mutate(Group_ID = factor(Group_ID, levels = cond_order)) %>% 
    group_by(Group_ID, Mouse_ID, Female_ID, Date) %>% 
    summarize(tot_calls = n()) %>% 
    left_join(fname_df, .) %>% 
    mutate(tot_calls = replace_na(tot_calls, 0)) %>% 
    mutate(Group_ID = factor(Group_ID, levels = cond_order)) %>% 
    arrange(Group_ID, desc(tot_calls))
  
  # # Comment out to ignore (add to data check object)
  # {
  #   low_BSL_ani <- t_data %>% filter(Group_ID == cond_order[1], tot_calls < min_calls)
  #   if (nrow(low_BSL_ani)) {
  #     print(low_BSL_ani)
  #     stop(paste("^ Consider excluding above from dataset. <", min_calls,"calls in:", cond_order[1]))
  #   }
  # }
  
  mouse_order <- t_data %>% pull(Mouse_ID) %>% unique()
  
  # dplyr::count number of conditions with same animal
  ani_rep <- t_data %>% 
    group_by(Mouse_ID) %>% 
    summarize(n=n())
  
  sort_df <- raw_df %>% 
    mutate(Group_ID = factor(Group_ID, levels = cond_order), 
           Label = factor(Label, levels = call_order),
           Mouse_ID = factor(Mouse_ID, levels = mouse_order)
    ) %>% 
    arrange(Group_ID, Mouse_ID, `Begin Time (s)`)
  
  sort_df %>% distinct(Group_ID, Mouse_ID) %>% group_by(Group_ID) %>% summarize(n()) %>% print()
}
{ 
  # filter animals w < min_calls for sequence analysis
  low_call_ani <- t_data %>% filter(tot_calls < min_calls) %>% pull(Mouse_ID)
  filt_mouse_order <- t_data %>% filter(!Mouse_ID %in% low_call_ani) %>% pull(Mouse_ID) %>% unique()
  filt_df <- sort_df %>% filter(!Mouse_ID %in% low_call_ani) %>% 
    mutate(Mouse_ID = factor(Mouse_ID, levels = filt_mouse_order))
  filt_t_data <- t_data %>% filter(!Mouse_ID %in% low_call_ani) %>% 
    mutate(Mouse_ID = factor(Mouse_ID, levels = filt_mouse_order))
  
  cat(paste("(filt_df & filt_mouse_order), excludes animals <", min_calls,"calls:"), low_call_ani, fill=T)
  cat(paste("Remaining animals:"), fill=T)
  filt_df %>% distinct(Group_ID, Mouse_ID) %>% group_by(Group_ID) %>% summarize(n()) %>% print()
}

# 1.7 Inter Silence Interval Distributions ======================================
# Load data & calculate ISIs
{
  library(mclust)
  library(tidyr)
  
  # Check for negative ISIs (Unable to fit GMM with negative values) ----
  neg_ISIs <- raw_df %>% select(Expt:ID, ISI) %>% filter(ISI<10)
  if (nrow(neg_ISIs)>0) {
    cat("\n#### df contains ISIs <0 ####\n", fill = T)
    print(neg_ISIs)
    
    ISI_df <- filt_df %>% select(Expt:ID, ISI) %>% 
      filter(!is.na(ISI), ISI>0)
    
    cat("(ISI_df) Removed -ve ISIs", fill = T)
  } else {
    ISI_df <- filt_df %>% select(Expt:ID, ISI) 
  }
  # uncomment to save neg_ISI.csvs
  {
    output_name <- file.path(plot_dir, paste0(expt,"_","neg_ISIs.csv"))
    write.csv(neg_ISIs, output_name, row.names = FALSE)
  }
  
  # Plot raw ISI distribution ----
  {
    mean_lines <- filt_df %>%
      group_by(Group_ID) %>%
      summarise(mean_time = mean(ISI, na.rm = TRUE),
                median_time = median(ISI, na.rm = TRUE),
                sd_time = sd(ISI, na.rm = TRUE),
                ref_ISI = min_IGI,
                ref_IBI = min_IBI,
                alt_IBI = .200
      )
    
    title <- paste0(sort_df$Expt[1], "_ISIdist")
    # design <- create_facet_design(filt_t_data)
    
    filt_df %>% 
      group_by(Group_ID, Mouse_ID) %>% 
      # mutate(ISI = ISI*1000) %>% 
      ggplot(aes(x = ISI, after_stat(density))) +
      # geom_histogram(aes(fill = Mouse_ID),position = "dodge", alpha = 0.8, binwidth = 0.1) +
      geom_density(aes(color = Mouse_ID, fill = Mouse_ID), alpha = 0.2) +  # Add density plot
      # Vertical mean line per Group
      geom_vline(data = mean_lines, aes(xintercept = mean_time),  # mapped from summary table
                 linetype = "dashed", color = "black", linewidth = 0.6) +
      # Add mean value as text next to the line
      geom_text(data = mean_lines, aes(x = mean_time, y = 10, label = paste0("mean = ",scales::number(mean_time, accuracy = 0.1))), 
                angle = 90, vjust = 1.5, hjust = 0, size = 3) +
      # Vertical median line per Group
      geom_vline(data = mean_lines, aes(xintercept = median_time),  # mapped from summary table
                 linetype = "dashed", color = "red", linewidth = 0.6) +
      # Add median value as text next to the line
      geom_text(data = mean_lines, 
                aes(x = median_time, y = 10, 
                    label = paste0("median = ",scales::number(median_time, accuracy = 0.01))), 
                angle = 90, vjust = -1.5, hjust = 0, size = 3, color = "red") +
      # reference values (from Castellucci, 2018)
      geom_vline(data = mean_lines, aes(xintercept = ref_ISI),  # mapped from summary table
                 linetype = "dashed", color = "blue", linewidth = 0.6) +
      geom_text(data = mean_lines, aes(x = ref_ISI, y = 10, label = paste0("ISI = ", mean_lines$ref_ISI[1])), 
                angle = 90, vjust = 1.5, hjust = -0.3, size = 3, color = "blue") +
      geom_vline(data = mean_lines, aes(xintercept = ref_IBI),  # mapped from summary table
                 linetype = "dashed", color = "blue", linewidth = 0.6) +
      geom_text(data = mean_lines, aes(x = ref_IBI, y = 10, label = paste0("IBI = ", mean_lines$ref_IBI[1])), 
                angle = 90, vjust = 1.5, hjust = -0.3, size = 3, color = "blue") +
      # geom_vline(data = summary, aes(xintercept = alt_IBI),  # mapped from summary table
      #            linetype = "dashed", color = "blue", linewidth = 0.6) +
      # geom_text(data = summary, aes(x = alt_IBI, y = 10, label = paste0("IBI = ", summary$alt_IBI[1])), 
      #           angle = 90, vjust = 1.5, hjust = -0.3, size = 3, color = "blue") +
      # facet_manual(~Group_ID, design = design) +
      facet_wrap(Group_ID~., nrow = n_cond, labeller = label_wrap_gen(width = 20)) +
      # scale_x_log10(labels = label_number()) +  # Adding logarithmic scale to y-axis
      xlim(0,1)+
      theme_classic() + 
      labs(y = "Density",
           x = "Inter Syllable Interval (s)",
           fill = "Animal",
           color = "Animal",
           title = title
      ) 
    # theme(legend.position = c(0.125,0.2),)
    
    p_path <- file.path(plot_dir, paste0("0_",title, ".png"))
    cat("Saving", p_path, fill=T)
    # ggsave(p_path, width = 11, height = 6*n_cond ,
    #        dpi = 300, bg = "white", units = "cm") # for single column plot
    # ggsave(p_path, height = 8, width = 6*n_cond ,
    #        dpi = 300, bg = "white", units = "cm") # for single column plot
    ggsave(p_path, width = 8*3, height = 8*2 ,
           dpi = 300, bg = "white", units = "cm") # for staggered facet
  }
  
  # Plot raw ISI distribution (Individual) -----
  {
    mean_lines <- filt_df %>%
      group_by(Mouse_ID) %>%
      summarise(mean_time = mean(ISI, na.rm = TRUE),
                median_time = median(ISI, na.rm = TRUE),
                sd_time = sd(ISI, na.rm = TRUE),
                ref_ISI = min_IGI,
                ref_IBI = min_IBI,
                alt_IBI = .200
      )
    
    filt_t_data %>% group_by(Group_ID) %>% summarize(n())
    title <- paste0(sort_df$Expt[1], "_ISIdist_Indiv")
    library(ggh4x)
    
    design <- create_facet_design(filt_t_data)
    
    filt_df %>% 
      group_by(Group_ID, Mouse_ID) %>% 
      # mutate(ISI = ISI*1000) %>% 
      ggplot(aes(x = ISI, after_stat(density))) +
      geom_histogram(aes(fill = Mouse_ID), color = "grey50",position = "dodge", size = .25, alpha = 0.5, binwidth = 0.01) +
      geom_density(aes(color = Mouse_ID), alpha = 1) +  # Add density plot
      geom_text(data = filt_t_data, aes(x = 0.55, y = 50, label=  paste("Total Calls:", tot_calls)),
                size = 3, color = "black") +
      # # Vertical mean line & text per Group
      # geom_vline(data = mean_lines, aes(xintercept = mean_time),  # mapped from summary table
      #   linetype = "dashed", color = "black", linewidth = 0.3) +
      # geom_text(data = mean_lines, 
      #           aes(x = mean_time, y = 15, 
      #               label = paste0("mean = ",scales::number(mean_time, accuracy = 0.1))), 
      #           angle = 90,               # vertical text beside the line
      #           vjust = -0.5,             # position above the base line
      #           hjust = -0.3,             # nudge text slightly right
      #           size = 3) +
      # Vertical median line & text per Group
      geom_vline(data = mean_lines, aes(xintercept = median_time),  # mapped from summary table
                 linetype = "dashed", color = "red", linewidth = 0.3) +
      geom_text(data = mean_lines, 
                aes(x = median_time, y = 15, 
                    label = paste0("median = ",scales::number(median_time, accuracy = 0.01))), 
                angle = 90,               # vertical text beside the line
                vjust = -0.5,             # position above the base line
                hjust = -0.3,             # nudge text slightly right
                size = 3, color = "red") +
      # reference values (from Castellucci, 2018)
      geom_vline(data = mean_lines, aes(xintercept = ref_ISI),  # mapped from summary table
                 linetype = "dashed", color = "blue", linewidth = 0.3) +
      geom_text(data = mean_lines, aes(x = ref_ISI, y = 10, label = paste0("ISI = ", mean_lines$ref_ISI[1])), 
                angle = 90, vjust = 1.5, hjust = -0.3, size = 3, color = "blue") +
      geom_vline(data = mean_lines, aes(xintercept = ref_IBI),  # mapped from summary table
                 linetype = "dashed", color = "blue", linewidth = 0.3) +
      geom_text(data = mean_lines, aes(x = ref_IBI, y = 10, label = paste0("IBI = ", mean_lines$ref_IBI[1])), 
                angle = 90, vjust = 1.5, hjust = -0.3, size = 3, color = "blue") +
      # geom_vline(data = summary, aes(xintercept = alt_IBI),  # mapped from summary table
      #            linetype = "dashed", color = "blue", linewidth = 0.3) +
      # geom_text(data = summary, aes(x = alt_IBI, y = 10, label = paste0("IBI = ", summary$alt_IBI[1])), 
      #           angle = 90, vjust = 1.5, hjust = -0.3, size = 3, color = "blue") +
      facet_manual(Group_ID~Mouse_ID, design = design) +
      # facet_wrap(.~Mouse_ID, ncol = 5, labeller = label_wrap_gen(width = 20)) +
      # scale_x_log10(labels = label_number()) +  # Adding logarithmic scale to y-axis
      xlim(0,0.8)+
      theme_classic() + 
      labs(y = "Density",
           x = "Inter Syllable Interval (s)",
           fill = "Condition",
           color = "Condition",
           title = title
      ) +
      theme(legend.position = "none",)
    max_ani <- t_data %>% group_by(Group_ID) %>% summarize(n=n()) %>% pull(n) %>% max()
    
    p_path <- file.path(plot_dir, paste0("0_",title, ".png"))
    cat("Saving", p_path, fill=T)
    ggsave(p_path, width = 4*max_ani, height = 6*n_cond , 
           dpi = 300, bg = "white", units = "cm")
  }
  
}


# 2. Analysis: Total Calls by Condition ======================================
{
  # Check total call number
  p_data <- sort_df %>% 
    group_by(Group_ID, Mouse_ID, Female_ID, Date) %>% 
    summarize(tot_calls = n()) %>% 
    mutate(is_below_min = tot_calls < min_calls)
  
  # Run lmer
  {
    if (any(ani_rep$n >1)) {
      cat("Some animals are in multiple conditions. Random variable Mouse_ID added for lmer", fill=T)
      model <- p_data %>%
        lmer(tot_calls ~ Group_ID + (1 | Mouse_ID)+ (1 | Female_ID), data = .)
    } else{
      model <- p_data %>% lmer(tot_calls ~ Group_ID + (1 | Female_ID), data = .)
    }
    
    # Pairwise comparisons using emmeans, with Holm adjustment
    em_df <- if (length(planned_contrasts) == 0) {
      pairs(emmeans(model, "Group_ID"), adjust = "holm")} else {
        contrast(emmeans(model, "Group_ID"), method = planned_contrasts, adjust = "holm")}
    
    contrast_df <- em_df %>% 
      as.data.frame() %>%
      separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
      mutate(
        group1 = str_trim(str_replace_all(group1, "^\\(|\\)$", "")),
        group2 = str_trim(str_replace_all(group2, "^\\(|\\)$", "")),
        p.label = symnum(p.value, corr = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns")),
        y.position = max(p_data$tot_calls, na.rm = TRUE) * 1.1
      )
    
    print(contrast_df)
    }
  {
    stat_df <- contrast_df %>% filter(p.label != "ns")
    
    # Create plot
    title <- paste0(sort_df$Expt[1], "_TotCalls")
    p_data %>% 
      ggplot(aes(x = Group_ID, y = tot_calls)) +
      # {if (any(ani_rep$n >1)) {geom_line(aes(group = Mouse_ID), size = 0.5, alpha = 0.8, color = "grey")}} +
      stat_summary(aes(fill = Group_ID),fun = mean, geom = "bar", width = 0.6, alpha = 0.8, color = "black") +  # mean bars
      # geom_jitter(width = 0.1, alpha = 0.4, size = 1) +  # show individual points\
      # # uncomment to annotate cutoff line
      geom_hline(aes(yintercept = min_calls), color = "red", linetype = "dashed")+
      # annotate("text", x = Inf, y = min_calls, label = paste0("min calls:", min_calls), 
      #          hjust = 1.1, vjust = -0.7, color = "red",size = 2, fontface = "italic")+
      geom_jitter(aes(color = is_below_min), width = 0.1, alpha = 0.4, size = 1) +  # color by new column
      # (rest of your code unchanged)
      scale_color_manual(
        values = c(`TRUE` = "red", `FALSE` = "black"),
        guide = FALSE  # removes the color legend for jittered points, optional
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +       # SE error bars
      scale_fill_manual(values = sel_pal) +  # 👈 apply your custom palette
      scale_x_discrete(labels = group_labels) + 
      theme_minimal()+
      labs(
        title = title,
        x = "Condition",
        y = "Total Number of Calls",
        fill = "Condition") +
      {if (nrow(stat_df)>0) {
        stat_pvalue_manual(
          data = stat_df, label = "p.label", xmin = "group1", xmax = "group2",
          y.position = "y.position", tip.length = 0.01, step.increase = 0.1
        )
      }} +
      theme(
        # axis.text.x = element_text(angle = 60, hjust = 0.5),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0.2, "cm")
      )
  }
  
  p_path <- file.path(plot_dir, paste0("1_",title,"_fin.png"))
  cat("Saving", p_path, fill=T)
  ggsave(p_path, width = 3+1.75*n_cond, height = 10 , 
         dpi = 300, bg = "white", units = "cm")
}

# 2.1 Analysis: Total Calls by Condition (removed low callers) ======================================
{
  # Check total call number
  p_data <- filt_df %>% 
    group_by(Group_ID, Mouse_ID, Female_ID, Date) %>% 
    summarize(tot_calls = n()) %>% 
    mutate(is_below_min = tot_calls < min_calls)
  
  # Run lmer
  {
    if (any(ani_rep$n >1)) {
      cat("Some animals are in multiple conditions. Random variable Mouse_ID added for lmer", fill=T)
      model <- p_data %>%
        lmer(tot_calls ~ Group_ID + (1 | Mouse_ID)+ (1 | Female_ID), data = .)
    } else{
      model <- p_data %>% lmer(tot_calls ~ Group_ID + (1 | Female_ID), data = .)
    }
    
    # Pairwise comparisons using emmeans, with Holm adjustment
    em_df <- if (length(planned_contrasts) == 0) {
      pairs(emmeans(model, "Group_ID"), adjust = "holm")} else {
        contrast(emmeans(model, "Group_ID"), method = planned_contrasts, adjust = "holm")}
    
    contrast_df <- em_df %>% 
      as.data.frame() %>%
      separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
      mutate(
        group1 = str_trim(str_replace_all(group1, "^\\(|\\)$", "")),
        group2 = str_trim(str_replace_all(group2, "^\\(|\\)$", "")),
        p.label = symnum(p.value, corr = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns")),
        y.position = max(p_data$tot_calls, na.rm = TRUE) * 1.1
      )
    
    print(contrast_df)
    }
  {
    stat_df <- contrast_df %>% filter(p.label != "ns")
    
    # Create plot
    title <- paste0(sort_df$Expt[1], "_TotCalls_min", min_calls,"Calls")
    p_data %>% 
      ggplot(aes(x = Group_ID, y = tot_calls)) +
      # {if (any(ani_rep$n >1)) {geom_line(aes(group = Mouse_ID), size = 0.5, alpha = 0.8, color = "grey")}} +
      stat_summary(aes(fill = Group_ID),fun = mean, geom = "bar", width = 0.6, alpha = 0.8, color = "black") +  # mean bars
      # geom_jitter(width = 0.1, alpha = 0.4, size = 1) +  # show individual points\
      # # uncomment to annotate cutoff line
      # annotate("text", x = Inf, y = min_calls, label = paste0("min calls:", min_calls), 
      #          hjust = 1.1, vjust = -0.7, color = "red",size = 2, fontface = "italic")+
      geom_jitter(aes(color = is_below_min), width = 0.1, alpha = 0.4, size = 1) +  # color by new column
      # (rest of your code unchanged)
      scale_color_manual(
        values = c(`TRUE` = "red", `FALSE` = "black"),
        guide = FALSE  # removes the color legend for jittered points, optional
      ) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +       # SE error bars
      scale_fill_manual(values = sel_pal) +  # 👈 apply your custom palette
      scale_x_discrete(labels = group_labels) + 
      theme_minimal()+
      labs(
        title = title,
        x = "Condition",
        y = "Total Number of Calls",
        fill = "Condition") +
      {if (nrow(stat_df)>0) {
        stat_pvalue_manual(
          data = stat_df, label = "p.label", xmin = "group1", xmax = "group2",
          y.position = "y.position", tip.length = 0.01, step.increase = 0.1
        )
      }} +
      theme(
        # axis.text.x = element_text(angle = 60, hjust = 0.5),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0.2, "cm")
      )
  }
  
  p_path <- file.path(plot_dir, paste0("1-1_",title,"_fin.png"))
  cat("Saving", p_path, fill=T)
  ggsave(p_path, width = 3+1.75*n_cond, height = 10 , 
         dpi = 300, bg = "white", units = "cm")
}

# 2.2 Analysis: Individual Syllable # by Condition ======================================
{
  # Set call plotting order
  p_rows <- 2 # number of rows to plot other calls
  p_call_order <- c(
    'Upward', 
    'Downward', 
    'Chevron', 
    'Reverse Chevron',
    'Complex', 
    'Flat', 
    '1 Frequency Step', 
    '2 Frequency Step', 
    '3+ Frequency Step',
    'Short', 
    'Harmonics', 
    'Composite',
    'Unstructured',
    'Murmur',
    'Noisy'
  )
  
  # # Option 2: uncomment to sort by freq of occurence
  # {
  #   freq_based_order <- sort_df %>%
  #     group_by(Group_ID, Mouse_ID, Female_ID, Label) %>%
  #     summarize(tot_calls = n()) %>%
  #     filter(Group_ID == "WT") %>% # uncomment to sort based on WT freq only
  #     arrange(desc(tot_calls)) %>%
  #     pull(Label) %>% as.character() %>% unique()
  # 
  #   p_call_order <- freq_based_order
  # }
  
  # Check individual syllable number
  p_data <- filt_df %>% 
    group_by(Group_ID, Mouse_ID, Female_ID, Label) %>% 
    summarize(tot_calls = n()) %>% 
    ungroup() %>% 
    mutate(Mouse_ID = factor(Mouse_ID, levels = mouse_order)) %>% 
    group_by(Group_ID, Mouse_ID, Female_ID) %>% 
    complete(Label, fill = list(tot_calls = 0)) %>% 
    ungroup() %>% 
    mutate(Label = factor(Label, levels = p_call_order))
  
  max_ani <- p_data %>% distinct(Group_ID, Mouse_ID) %>% 
    dplyr::count(Group_ID) %>% pull(n) %>% max()
  
  max_call <- p_data %>% pull(tot_calls) %>% max()
  
  # Run lmer comparisons for individual call types ----
  y_pos_df <- p_data %>%
    group_by(Label) %>%
    summarise(y.position = max(tot_calls) + 3)
  
  model_results <- p_data %>% 
    group_by(Label) %>%
    group_nest() %>%
    mutate(
      model = purrr::map(data, function(df) {
        if (any(ani_rep$n > 1)) {
          lmer(tot_calls ~ Group_ID + (1 | Mouse_ID) + (1 | Female_ID), data = df)
        } else {
          lmer(tot_calls ~ Group_ID + (1 | Female_ID), data = df)
        }
      }),
      # Compute estimated marginal means
      emmeans_res = purrr::map(model, ~ emmeans(.x, "Group_ID")),
      
      # Conditionally switch between contrast types
      emmeans_pairwise = purrr::map(emmeans_res, function(em) {
        if (length(planned_contrasts) == 0) {
          pairs(em, adjust = "holm")
        } else {
          contrast(em, method = planned_contrasts, adjust = "holm")
        }
      }),
      pairs_df = purrr::map(emmeans_pairwise, ~ as.data.frame(.x))) %>% 
    select(Label, pairs_df) %>%
    unnest(pairs_df) %>% 
    separate(contrast, into = c("group1", "group2"), sep = " - ", remove= F) %>%
    mutate(
      group1 = str_trim(group1),
      group2 = str_trim(group2),
      p.label = symnum(p.value, corr = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*", "ns")
      )) %>% 
    mutate(
      group1 = str_replace_all(group1, "^\\(|\\)$", ""),
      group2 = str_replace_all(group2, "^\\(|\\)$", "")
    ) %>% 
    left_join(y_pos_df, by = "Label")
  
  # Create plots
  title <- paste0(sort_df$Expt[1], "_IndivSyllables_min", min_calls,"Calls")
  stat_df <- model_results %>% filter(p.label != "ns")
  p_data %>% 
    separate(Group_ID, into = c("Treatment", "Timepoint"), sep = "-", fill = "right", remove = F) %>% 
    mutate(
      Timepoint = if_else(is.na(Timepoint), "BSL", Timepoint),
      Timepoint = factor(Timepoint, levels = c("BSL", "30Min", "4Hr"))
    ) %>% 
    ggplot(aes(x = Group_ID, y = tot_calls)) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 0.5) +  # show individual points
    # {if (any(ani_rep$n >1)) {geom_line(aes(group = Mouse_ID), size = 0.5, alpha = 0.8, color = "grey")}} +
    stat_summary(aes(fill = Group_ID),fun = mean, geom = "bar", width = 0.6, alpha = 0.8, color = "black") +  # mean bars
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +       # SE error bars
    scale_fill_manual(values = sel_pal) +  # 👈 apply your custom palette
    scale_x_discrete(labels = group_labels) + 
    theme_minimal()+
    facet_wrap(.~Label, scales = "free_y", nrow = p_rows)+
    labs(
      title = title,
      x = "Condition",
      y = "Total Number of Calls",
      fill = "Condition") +
    {if (nrow(stat_df)>0) {
      stat_pvalue_manual(
        data = stat_df, label = "p.label", xmin = "group1", xmax = "group2",
        y.position = "y.position", tip.length = 0.01, step.increase = 0, size = 5
      )
    }} +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0.5, size = 7),
      legend.position = c(0.925,0.2), # nrows =2
      # legend.position = c(0.925,0.1),
      # axis.text.x = element_text(size = 6.5),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing = unit(0.2, "cm")
    )
  
  p_path <- file.path(plot_dir, paste0("2_", title,".png"))
  cat("Saving", p_path, fill=T)
  ggsave(p_path, 
         # width = 5+2*length(call_order), height = 15 , # 2 rows
         width = 5+4*ceiling(length(p_call_order)/p_rows), height = 10 *p_rows, 
         dpi = 300, bg = "white", units = "cm")
}

# 2.2.1 Analysis: Individual Syllable % by Condition ======================================
{
  # Set call plotting order
  p_rows <- 2 # number of rows to plot other calls
  p_call_order <- c(
    'Upward', 
    'Downward', 
    'Chevron', 
    'Reverse Chevron',
    'Complex', 
    'Flat', 
    '1 Frequency Step', 
    '2 Frequency Step', 
    '3+ Frequency Step',
    'Short', 
    'Harmonics', 
    'Composite',
    'Unstructured',
    'Murmur',
    'Noisy'
  )
  
  # Get total call number
  tot_df <- filt_df %>% 
    group_by(Group_ID, Mouse_ID, Female_ID, Date) %>% 
    summarize(tot_calls = n()) %>% 
    mutate(is_below_min = tot_calls < min_calls)
  
  # Check individual syllable number
  p_data <- filt_df %>% 
    group_by(Group_ID, Mouse_ID, Female_ID, Label) %>% 
    summarize(indiv_calls = n()) %>% 
    left_join(tot_df) %>% 
    mutate(per_calls = indiv_calls/tot_calls*100) %>% 
    ungroup() %>% 
    mutate(Mouse_ID = factor(Mouse_ID, levels = mouse_order)) %>% 
    group_by(Group_ID, Mouse_ID, Female_ID) %>% 
    complete(Label, fill = list(per_calls = 0)) %>% 
    ungroup() %>% 
    mutate(Label = factor(Label, levels = p_call_order))
  
  max_ani <- p_data %>% distinct(Group_ID, Mouse_ID) %>% 
    dplyr::count(Group_ID) %>% pull(n) %>% max()
  
  max_call <- p_data %>% pull(per_calls) %>% max()
  
  # Run lmer comparisons for individual call types ----
  y_pos_df <- p_data %>%
    group_by(Label) %>%
    summarise(y.position = max(per_calls) + 3)
  
  model_results <- p_data %>% 
    group_by(Label) %>%
    group_nest() %>%
    mutate(
      model = purrr::map(data, function(df) {
        if (any(ani_rep$n > 1)) {
          lmer(per_calls ~ Group_ID + (1 | Mouse_ID) + (1 | Female_ID), data = df)
        } else {
          lmer(per_calls ~ Group_ID + (1 | Female_ID), data = df)
        }
      }),
      # Compute estimated marginal means
      emmeans_res = purrr::map(model, ~ emmeans(.x, "Group_ID")),
      
      # Conditionally switch between contrast types
      emmeans_pairwise = purrr::map(emmeans_res, function(em) {
        if (length(planned_contrasts) == 0) {
          pairs(em, adjust = "holm")
        } else {
          contrast(em, method = planned_contrasts, adjust = "holm")
        }
      }),
      pairs_df = purrr::map(emmeans_pairwise, ~ as.data.frame(.x))) %>% 
    select(Label, pairs_df) %>%
    unnest(pairs_df) %>% 
    separate(contrast, into = c("group1", "group2"), sep = " - ", remove= F) %>%
    mutate(
      group1 = str_trim(group1),
      group2 = str_trim(group2),
      p.label = symnum(p.value, corr = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*", "ns")
      )) %>% 
    mutate(
      group1 = str_replace_all(group1, "^\\(|\\)$", ""),
      group2 = str_replace_all(group2, "^\\(|\\)$", "")
    ) %>% 
    left_join(y_pos_df, by = "Label")
  
  # Create plots
  title <- paste0(sort_df$Expt[1], "_IndivSyllablesPercent_min", min_calls,"Calls")
  stat_df <- model_results %>% filter(p.label != "ns")
  p_data %>% 
    separate(Group_ID, into = c("Treatment", "Timepoint"), sep = "-", fill = "right", remove = F) %>% 
    mutate(
      Timepoint = if_else(is.na(Timepoint), "BSL", Timepoint),
      Timepoint = factor(Timepoint, levels = c("BSL", "30Min", "4Hr"))
    ) %>% 
    ggplot(aes(x = Group_ID, y = per_calls)) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 0.5) +  # show individual points
    # {if (any(ani_rep$n >1)) {geom_line(aes(group = Mouse_ID), size = 0.5, alpha = 0.8, color = "grey")}} +
    stat_summary(aes(fill = Group_ID),fun = mean, geom = "bar", width = 0.6, alpha = 0.8, color = "black") +  # mean bars
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +       # SE error bars
    scale_fill_manual(values = sel_pal) +  # 👈 apply your custom palette
    scale_x_discrete(labels = group_labels) + 
    theme_minimal()+
    facet_wrap(.~Label, scales = "free_y", nrow = p_rows)+
    labs(
      title = title,
      x = "Condition",
      y = "% of Total Calls",
      fill = "Condition") +
    {if (nrow(stat_df)>0) {
      stat_pvalue_manual(
        data = stat_df, label = "p.label", xmin = "group1", xmax = "group2",
        y.position = "y.position", tip.length = 0.01, step.increase = 0, size = 5
      )
    }} +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0.5, size = 7),
      legend.position = c(0.925,0.2), # nrows =2
      # legend.position = c(0.925,0.1),
      # axis.text.x = element_text(size = 6.5),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing = unit(0.2, "cm")
    )
  
  p_path <- file.path(plot_dir, paste0("2-2_", title,".png"))
  cat("Saving", p_path, fill=T)
  ggsave(p_path, 
         # width = 5+2*length(call_order), height = 15 , # 2 rows
         width = 5+4*ceiling(length(p_call_order)/p_rows), height = 10 *p_rows, 
         dpi = 300, bg = "white", units = "cm")
}

# 2.3 Analysis: Individual Syllable Variable by Condition ======================================
{
  # Set call plotting order
  p_call_order <- c(
    'Upward', 
    'Chevron', 
    'Harmonics', 
    '1 Frequency Step', 
    '2 Frequency Step', 
    '3+ Frequency Step',
    'Short', 
    'Downward', 
    'Complex', 
    'Flat', 
    'Composite',
    'Reverse Chevron',
    'Unstructured',
    'Murmur',
    'Noisy'
  )
  
  # # Option 2: uncomment to sort by freq of occurence
  # {
  #   freq_based_order <- sort_df %>%
  #     group_by(Group_ID, Mouse_ID, Female_ID, Label) %>%
  #     summarize(tot_calls = n()) %>%
  #     filter(Group_ID == "WT") %>% # uncomment to sort based on WT freq only
  #     arrange(desc(tot_calls)) %>%
  #     pull(Label) %>% as.character() %>% unique()
  # 
  #   p_call_order <- freq_based_order
  # }
  
  # Check individual syllable number
  p_data <- filt_df %>% 
    group_by(Group_ID, Mouse_ID, Female_ID, Label) %>% 
    summarize(`Principal Freq (kHz)` = mean(`Principal Frequency (kHz)`),
              `Power (dB/Hz)` = mean(`Mean Power (dB/Hz)`),
              `Delta Freq (kHz)` =  mean(`Delta Freq (kHz)`),
              `Call Duration (ms)` =  mean(`Call Length (s)`*1000),
              `Sinuosity` =  mean(`Sinuosity`),
              `Tonality` =  mean(`Tonality`)
    ) %>% 
    ungroup() %>% 
    mutate(Mouse_ID = factor(Mouse_ID, levels = mouse_order)) %>% 
    mutate(Label = factor(Label, levels = p_call_order))
  
  max_ani <- p_data %>% distinct(Group_ID, Mouse_ID) %>% 
    dplyr::count(Group_ID) %>% pull(n) %>% max()
  
  for (var in colnames(p_data)[5:ncol(p_data)]) {
    # Run lmer comparisons for each var ----
    y_pos_df <- p_data %>%
      group_by(Label) %>%
      summarise(y.position = max(!!sym(var))*1.15)
    
    model_results <- p_data %>% 
      group_by(Label) %>%
      group_nest() %>%
      mutate(
        model = purrr::map(data, function(df) {
          if (any(ani_rep$n > 1)) {
            lmer(!!sym(var) ~ Group_ID + (1 | Mouse_ID) + (1 | Female_ID), data = df)
          } else {
            lmer(!!sym(var) ~ Group_ID + (1 | Female_ID), data = df)
          }
        }),
        # Compute estimated marginal means
        emmeans_res = purrr::map(model, ~ emmeans(.x, "Group_ID")),
        
        # Conditionally switch between contrast types
        emmeans_pairwise = purrr::map(emmeans_res, function(em) {
          if (length(planned_contrasts) == 0) {
            pairs(em, adjust = "holm")
          } else {
            contrast(em, method = planned_contrasts, adjust = "holm")
          }
        }),
        pairs_df = purrr::map(emmeans_pairwise, ~ as.data.frame(.x))) %>% 
      select(Label, pairs_df) %>%
      unnest(pairs_df) %>% 
      separate(contrast, into = c("group1", "group2"), sep = " - ", remove= F) %>%
      mutate(
        group1 = str_trim(group1),
        group2 = str_trim(group2),
        p.label = symnum(p.value, corr = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns")
        )) %>% 
      mutate(
        group1 = str_replace_all(group1, "^\\(|\\)$", ""),
        group2 = str_replace_all(group2, "^\\(|\\)$", "")
      ) %>% 
      left_join(y_pos_df, by = "Label")
    
    # Create plots
    title <- paste0(sort_df$Expt[1], "_", var, "_IndivSyllables")
    stat_df <- model_results %>% filter(p.label != "ns") %>% 
      group_by(Label) %>%
      mutate(n_comparisons = n(), 
             y.position = ifelse(n_comparisons > 1, y.position + y.position* 0.15 * (row_number() - 1), y.position)) %>%
      ungroup() %>%
      select(-n_comparisons) 
    
    p_data %>% 
      separate(Group_ID, into = c("Genotype", "Treatment"), sep = "-", fill = "right", remove = F) %>% 
      mutate(
        Treatment = if_else(is.na(Treatment), "BSL", Treatment),
        Treatment = factor(Treatment, levels = c("BSL", "30Min", "4Hr"))
      ) %>% 
      ggplot(aes(x = Group_ID, y = !!sym(var))) +
      geom_boxplot(aes(fill = Group_ID), width = 0.1, alpha = 1) +
      geom_violin(aes(color = Group_ID), width = 0.5, alpha = 0.5) +
      geom_jitter(width = 0.1, alpha = 0.4, size = 0.5) +  # show individual points
      # {if (any(ani_rep$n >1)) {geom_line(aes(group = Mouse_ID), size = 0.5, alpha = 0.8, color = "grey")}} +
      # stat_summary(aes(fill = Group_ID),fun = mean, geom = "bar", width = 0.6, alpha = 0.8, color = "black") +  # mean bars
      # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +       # SE error bars
      scale_fill_manual(values = sel_pal) +
      scale_color_manual(values = sel_pal) +
      scale_x_discrete(labels = group_labels) + 
      theme_minimal()+
      facet_wrap(.~Label, scales = "free_y", nrow = 3)+
      labs(
        title = title,
        x = "Condition",
        y = var,
        fill = "Condition",
        color = "Condition") +
      {if (nrow(stat_df)>0) {
        stat_pvalue_manual(
          data = stat_df, label = "p.label", xmin = "group1", xmax = "group2",
          y.position = "y.position", tip.length = 0.0, step.increase = 0, size = 5
        )
      }} +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 0.5, size = 7),
        # legend.position = c(0.925,0.2), # nrows =2
        legend.position = c(0.925,0.1),
        # axis.text.x = element_text(size = 6.5),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0.2, "cm")
      )
    
    p_path <- file.path(plot_dir, paste0("2-3_", gsub("[^A-Za-z0-9_]", "_", gsub("[ ()/]", "", title)),".png"))
    cat("Saving", p_path, fill=T)
    ggsave(p_path, 
           # width = 5+2*length(call_order), height = 15 , # 2 rows
           width = 5+2*9, height = 20 , 
           dpi = 300, bg = "white", units = "cm")
    
    model_results %>% 
      mutate(Variable = var,.before = 1) %>% 
      write.csv(., sub("\\.png$", ".csv", p_path))
  }
}


# 3.1 Prepare bout dataset ====================================================
# Initialize relevant dataframes
{
  min_calls_in_seq <- 3
  seq_order <- c('IBI', call_order)
  
  bout_call_order <- c(
    'Short', 
    'Upward', 
    'Downward', 
    '1 Frequency Step', 
    '2 Frequency Step', 
    '3+ Frequency Step',
    'Flat', 
    'Harmonics', 
    'Complex', 
    'Chevron', 
    'Reverse Chevron',
    'Composite',
    'Unstructured',
    'Murmur',
    'Noisy'
  )
  
  # Select 3 top call animals from each condition to plot
  sel_anis <- filt_t_data %>% group_by(Group_ID) %>% slice_max(tot_calls, n=3) %>% pull(Mouse_ID) %>% as.character()
  # sel_anis <- c("C57M8234","C57M8147", "C57M8035", "C57M8034")
  
  
  bout_t_data <- filt_t_data %>% 
    filter(tot_calls>min_calls) # Only include conditions when calls > min_calls (excludes by per experiment not by animal)
  
  bout_data <- filt_df %>% 
    semi_join(bout_t_data) %>% # Only include conditions when calls > min_calls
    group_by(Group_ID, Mouse_ID, bout_id) %>% 
    mutate(relative_start = `Begin Time (s)`[1]) %>% # Set relative start time
    mutate(
      Label = factor(Label, levels = bout_call_order), # Change order of first call type
      Start.Time = `Begin Time (s)`- relative_start,
      End.Time = `Begin Time (s)`- `End Time (s)`,
      Start.seq = coalesce(lag(call_seq), 0) # Fills with 0 if NA
    ) %>% 
    filter(!calls_in_bout<min_calls_in_seq) %>%
    arrange(Start.seq)
  
  shift_df <- bout_data %>% 
    filter(Mouse_ID %in% sel_anis) %>% # Select specific animal examples
    group_by(Group_ID, Mouse_ID, bout_id, Label) %>% 
    summarize(pos1 = min(call_seq)) %>% 
    mutate(Label_order = as.integer(Label), .keep = "unused") %>% 
    pivot_wider(names_from = Label_order, names_prefix = "Lab", values_from = pos1) %>% 
    mutate(across(everything(), ~replace(., is.na(.), 100)))
  
  shift_bout_data <- bout_data %>% left_join(shift_df) %>% 
    filter(Mouse_ID %in% sel_anis) %>% # Select specific animal examples
    mutate(shift = ifelse(
      Lab1<Lab2, -Lab1+1, ifelse(
        Lab2<Lab3, -Lab2+1, ifelse(
          Lab3<Lab4, -Lab3+1, ifelse(
            Lab4<Lab5, -Lab4+1, +1))))) %>% 
    mutate(Start.seq.L = Start.seq+shift, call_seq.L = call_seq+shift) %>%  # Shift start seq
    mutate(L1 = Label[1], L2 = Label[2], L3 = Label[3], L4 = Label[4])
}

# 3.2 Analysis: Bout Raster by Condition ======================================
{
  n_bout <- 150 #Set number of bouts to plot
  title <- paste0(shift_bout_data$Expt[1], "_next_seq_raster_minIBI", min_IBI)
  
  # Arrange plot order
  p_order <- shift_bout_data %>% 
    group_by(Group_ID, Mouse_ID, bout_id) %>% 
    slice_head(n=1) %>% 
    group_by(Group_ID) %>% 
    # arrange(L1, desc(calls_in_bout),) %>%
    arrange(desc(calls_in_bout),) %>%
    # arrange(L1, Lab1,  desc(calls_in_bout),) %>%
    # arrange(desc(Lab1), desc(Lab2), desc(Lab3), desc(Lab4), desc(calls_in_bout)) %>%
    mutate(bout_plot_y = row_number()) %>% 
    select(Group_ID, Mouse_ID, Female_ID, bout_id, bout_plot_y)
  
  #TODO: Add number of animals in raster (n in diff conditions are unbalanced)
  
  p_data <- shift_bout_data %>% 
    left_join(p_order) %>% 
    group_by(Group_ID, Mouse_ID, bout_id)
  
  p_data %>% group_by(Group_ID, Mouse_ID, bout_id) %>% slice_sample(n=1) %>% group_by(Group_ID, Mouse_ID) %>% summarize(n())
  p_tot <- p_data %>% ungroup() %>% distinct(Group_ID, Mouse_ID) %>% group_by(Group_ID) %>% summarize(n = n())
  bout_max <- p_data %>% pull(bout_plot_y) %>% max()
  seq_max <- p_data %>% pull(call_seq) %>% max()
  # Segment version of raster
  p_data %>% 
    ggplot() +
    geom_segment(aes(x = Start.seq, xend = call_seq, 
                     y = bout_plot_y, yend = bout_plot_y, 
                     color = Label), size = 1) +
    geom_text(data = p_tot, aes(x = seq_max*.5, y = bout_max*.9, label = paste0("n = ", n)), 
              vjust = 1.5, hjust = -0.3, size = 3, color = "black") +
    theme_minimal() +
    labs(
      title = paste0("Call sequences per bout (ISI>", min_IBI, "s)"),
      x = "Syllable Order",
      y = paste("Longest", n_bout, "Bouts")
    ) +
    scale_y_continuous(limits = c(NA, n_bout)) +
    # scale_y_reverse() +
    scale_color_manual(values = call_pal) +
    facet_wrap(Group_ID~., ncol = n_cond, labeller = label_wrap_gen(width = 20)) +
    # facet_manual(Group_ID~., design = "ABC\n#DE")+
    theme(
      axis.text.y = element_text(size = 6),
      panel.grid.major.y = element_blank(),
      legend.position = c(0.225, 0.8),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.size = unit(.5, "lines")
    ) +
    guides(y = "none", color = guide_legend(title = "Call Type"))
  
  max_ani <- p_data %>% ungroup() %>% distinct(Group_ID, Mouse_ID) %>% 
    dplyr::count(Group_ID) %>% pull(n) %>% max()
  
  p_path <- file.path(plot_dir, paste0("3_", title,".png"))
  cat("Saving", p_path, fill=T)
  ggsave(p_path, width = 3+5*n_cond, height = 12.5 , 
         dpi = 300, bg = "white", units = "cm")
}

# 4.1 Prepare Markov Chain Analysis dataset ====================================
# Initialize relevant dataframes
{  
  source("extractUSV.R")
  
  transition_tables <- getTransitionTables(bout_data, seq_order, min_calls_in_seq, min_IBI)
  
  # transition_tables$dplyr::count_Matrix[1]
  # transition_tables$Prob_Matrix[1]
  
  markov_results <- calMarkovEnt(transition_tables)
  sm_all<- markov_results$sm_all
  entrate_all <- markov_results$entrate_all %>% left_join(bout_t_data)
}

# 4.2 Analysis: EntR by Condition ==============================================
{
  p_data <- entrate_all
  
  # { # Uncomment to exclude inappropriate BSL condition
  #   p_data <- entrate_all %>% filter(!Group_ID == "BSL") 
  #   planned_contrasts <- list(
  #     "Veh-30Min - 0.32-30Min" = c("Veh-30min" = -1, "0.32-30Min" = 1, "Veh-4Hr" = 0, "0.32-4Hr" = 0),  # ✅ matched names
  #     "Veh-4Hr - 0.32-4Hr"     = c("Veh-30min" = 0, "0.32-30Min" =0, "Veh-4Hr" = -1, "0.32-4Hr" = 1)
  #   )
  # }
  
  # Run lmer
  {
    if (any(ani_rep$n >1)) {
      cat("Some animals are in multiple conditions. Random variable Mouse_ID added for lmer", fill=T)
      model <- p_data %>%
        lmer(entR ~ Group_ID + (1 | Mouse_ID)+ (1 | Female_ID), data = .)
    } else{
      model <- p_data %>% lmer(entR ~ Group_ID + (1 | Female_ID), data = .)
    }
    
    # Pairwise comparisons using emmeans, with Holm adjustment
    em_df <- if (length(planned_contrasts) == 0) {
      pairs(emmeans(model, "Group_ID"), adjust = "holm")} else {
        contrast(emmeans(model, "Group_ID"), method = planned_contrasts, adjust = "holm")}
    
    contrast_df <- em_df %>% 
      as.data.frame() %>%
      separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
      mutate(
        group1 = str_trim(str_replace_all(group1, "^\\(|\\)$", "")),
        group2 = str_trim(str_replace_all(group2, "^\\(|\\)$", "")),
        p.label = symnum(p.value, corr = FALSE,
                         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns")),
        y.position = max(p_data$tot_calls, na.rm = TRUE) * 1.1
      )
    
    print(contrast_df)
  }
  {
    stat_df <- contrast_df %>% filter(p.label != "ns")
    
    # Create plot
    title <- paste0(sort_df$Expt[1], "_EntR_", min_IBI, "s_minCalls", min_calls_in_seq)
    p_data %>% 
      ggplot(aes(x = Group_ID, y = entR)) +
      # stat_halfeye(adjust = .5, width = 0.75, .width = 0, justification = -.2, alpha = .8) +
      # {if (any(ani_rep$n >1)) {geom_line(aes(group = Mouse_ID), size = 0.5, alpha = 0.8, color = "grey")}} +
      geom_boxplot(aes(fill = Group_ID), width = 0.3, alpha = 0.5) +
      geom_jitter(width = 0.1, alpha = 0.4, size = 1) +  # show individual points
      # stat_summary(aes(fill = Group_ID),fun = mean, geom = "bar", width = 0.6, alpha = 0.8, color = "black") +  # mean bars
      # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +       # SE error bars
      scale_fill_manual(values = sel_pal) +  # 👈 apply your custom palette
      scale_x_discrete(labels = group_labels) + 
      theme_minimal()+
      labs(
        title = title,
        x = "Condition",
        y = "Markov Chain Entropy Rate",
        fill = "Condition") +
      {if (nrow(stat_df)>0) {
        stat_pvalue_manual(
          data = stat_df, label = "p.label", xmin = "group1", xmax = "group2",
          y.position = "y.position", tip.length = 0.01, step.increase = 0.1
        )
      }} +
      theme(
        # axis.text.x = element_text(angle = 60, hjust = 0.5),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing = unit(0.2, "cm")
      )
  }
  
  p_path <- file.path(plot_dir, paste0("4_", title,"_fin.png"))
  cat("Saving", p_path, fill=T)
  ggsave(p_path, width = 3+1.75*n_cond, height = 10 , 
         dpi = 300, bg = "white", units = "cm")
}


# 4.2 Prepare mean Syntax data ==============================================
{
  
  mean_mat <- function(mats) {
    Reduce("+", mats) / length(mats)
  }
  
  facet_ids <- LETTERS[1:n_cond] # e.g. A, B, C, etc.
  
  avg_transitions <- transition_tables %>%
    group_by(Group_ID) %>%
    summarize(
      matrices = list(Prob_Matrix)
    ) %>%  mutate(
      Avg_Prob_Matrix = lapply(matrices, mean_mat)
    ) %>% 
    mutate(facet_id = facet_ids, 
           facet_label = as.character(Group_ID))
}


# 4.2 Analysis: Syntax by Condition ==============================================
{
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggh4x)
  
  min_prob <- 0.1
  design <- group_design[[expt]]
  # design <- c("ABC")
  
  # Define label mappings
  node_labeller <- c(
    "IBI" = "Start/End",
    "Flat" = "Flat",
    "Downward" = "Down",
    "Upward" = "Up",
    "Composite" = "Cmpst",
    "Chevron" = "Chev",
    "Reverse Chevron" = "Rev \nChev",
    "1 Frequency Step" = "1 Freq\nStep",
    "2 Frequency Step" = "2 Freq\nStep",
    "3+ Frequency Step" = "3+ Freq\nStep",
    "Complex" = "Cmplx",
    "Harmonics" = "Hrm",
    "Short" = "Short",
    "Unstructured" = "Unstr",
    "Murmur" = "Murm",
    "Noisy" = "Nois"
  )
  
  # ---- Function to build graph data per animal ----
  build_graph_data <- function(mat, label_order, facet_id = NA, edge_weights = "prob", min_prob) {
    n_labels <- length(label_order)
    theta <- seq(from = pi, by = -2 * pi / n_labels, length.out = n_labels)
    
    # Circular node position
    nodes <- tibble(
      name = label_order,
      x = round(cos(theta), 10),
      y = round(sin(theta), 10),
      direction = (atan2(y, x) * 180 / pi) %% 360,
      facet_id = facet_id
    )
    
    edges <- as.data.frame(as.table(mat)) %>%
      filter(Freq > min_prob) %>%
      rename(from = Var1, to = Var2, !!edge_weights := Freq) %>%
      mutate(is_loop = from == to) %>% 
      mutate(
        from_idx = match(from, label_order),
        span = 90
      ) %>% 
      left_join(nodes %>% select(name, direction, facet_id),
                by = c("from" = "name")) %>%
      mutate(
        is_loop = from == to,
        strength = ifelse(as.character(from) < as.character(to), -0.4, 0.2)
        # strength = 0.2
      )
    
    list(nodes = nodes, edges = edges)
  }
  
  # Build graph data for each condition
  graph_data <- pmap(
    list(
      mat = avg_transitions$Avg_Prob_Matrix,
      facet_id = avg_transitions$facet_id,
      edge_weights = "prob"
    ),
    build_graph_data,
    label_order = seq_order,
    min_prob = min_prob
  )
  
  all_nodes <- bind_rows(purrr::map(graph_data, "nodes"))
  all_edges <- bind_rows(purrr::map(graph_data, "edges"))
  max_prob <- max(all_edges$prob)
  tg <- tbl_graph(nodes = all_nodes, edges = all_edges, directed = TRUE)
  
  #' Plotting code ----
  title_map <- setNames(avg_transitions$facet_label, avg_transitions$facet_id)
  cap_size <- circle(6, 'mm')
  title <- paste0(entrate_all$Expt[1], "_Syntax_ByCondition_minIBI", min_IBI)
  
  ggraph(tg, layout = "manual", x = x, y = y) +
    geom_edge_arc(
      aes(edge_width = prob, edge_colour = prob),
      arrow = arrow(length = unit(1, "mm"), type = "closed"),
      end_cap = cap_size,
      start_cap = cap_size,
      lineend = "round",
      strength = tg %>% activate(edges) %>% filter(!is_loop) %>% pull(strength)
    ) +
    geom_edge_loop(
      aes(edge_width = prob, edge_colour = prob, direction = direction, span = span),
      arrow = arrow(length = unit(1, "mm"), type = "closed"),
      end_cap = cap_size,
      start_cap = cap_size,
    ) +
    geom_node_point(aes(fill = name), size = 12, shape = 21, alpha = 0.8, show.legend = FALSE) +
    geom_node_text(aes(label = node_labeller[name]), size = 2, color = "black") +
    scale_edge_width(range = c(0.1, 2), guide = "none") +
    scale_edge_colour_gradientn(
      colours = rev(c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154")),
      limits = c(min_prob, max_prob),
      name = paste0("Transition\nProbability\n(>",min_prob, ")")
    ) +
    facet_manual(~facet_id, design = design, labeller = as_labeller(title_map), strip.position = "bottom") +
    coord_fixed() +
    theme_void() + 
    # guides(color="Syllables")+
    theme(
      strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 15),
      # legend.position = c(0.2,0.25),
      legend.title = element_text(size = 11, margin = margin(b = 15)),
      legend.text = element_text(size = 8),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.text.align = 0,
      plot.title = element_text(hjust = 0.5)
    )
  
  p_path <- file.path(plot_dir, paste0("4-1_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- strsplit(design, "\n")[[1]] %>% {data.frame(nrow = length(.), ncol = max(nchar(.)))}
  ggsave(p_path, width = 5+10*layout_df$ncol, height = 2+10*layout_df$nrow, 
         dpi = 300, bg = "transparent", units = "cm")
}

# 4.2 Analysis: Alternative Syntax by Condition (better) ==============================================
{
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggh4x)
  
  min_prob <- 0.1
  
  syn_order <- seq_order
  # # Uncomment & Reorder this to change arrangement of syllables on syntax circle
  # syn_order <- c(
  #   'IBI',
  #   'Short',
  #   'Upward',
  #   'Downward',
  #   'Flat',
  #   'Harmonics',
  #   'Complex',
  #   'Chevron',
  #   'Reverse Chevron',
  #   'Composite',
  #   'Unstructured',
  #   '1 Frequency Step',
  #   '2 Frequency Step',
  #   '3+ Frequency Step',
  #   'Murmur',
  #   'Noisy'
  # )
  
  # Define label mappings
  node_labeller <- c(
    "IBI" = "Start/\nEnd",
    "Flat" = "Flat",
    "Downward" = "Down",
    "Upward" = "Up",
    "Composite" = "Cmpst",
    "Chevron" = "Chev",
    "Reverse Chevron" = "Rev \nChev",
    "1 Frequency Step" = "1 Freq\nStep",
    "2 Frequency Step" = "2 Freq\nStep",
    "3+ Frequency Step" = "3+ Freq\nStep",
    "Complex" = "Cmplx",
    "Harmonics" = "Hrm",
    "Short" = "Short",
    "Unstructured" = "Unstr",
    "Murmur" = "Murm",
    "Noisy" = "Nois"
  )
  
  # ---- Function to build graph data per animal ----
  build_graph_data <- function(mat, label_order, facet_id = NA, edge_weights = "prob", min_prob) {
    n_labels <- length(syn_order)
    theta <- seq(from = pi, by = -2 * pi / n_labels, length.out = n_labels)
    offset <- 1.1 # how far out to place the label
    
    # Circular node position
    nodes <- tibble(
      name = syn_order,
      x = round(cos(theta), 10),
      y = round(sin(theta), 10),
      label_x = x * offset,
      label_y = y * offset,
      direction = (atan2(y, x) * 180 / pi) %% 360,
      hjust = ifelse(direction > 90 & direction < 270, 1, 0),
      angle = ifelse(direction > 90 & direction < 270, direction + 180, direction),
      facet_id = facet_id
    )
    
    edges <- as.data.frame(as.table(mat)) %>%
      filter(Freq > min_prob) %>%
      rename(from = Var1, to = Var2, !!edge_weights := Freq) %>%
      mutate(is_loop = from == to) %>% 
      mutate(
        from_idx = match(from, label_order)
      ) %>% 
      left_join(nodes %>% select(name, direction, facet_id),
                by = c("from" = "name")) %>%
      mutate(
        is_loop = from == to,
        strength = ifelse(as.character(from) < as.character(to), -0.4, 0.2)
        # strength = 0.2
      )
    
    list(nodes = nodes, edges = edges)
  }
  
  # Build graph data for each condition
  graph_data <- pmap(
    list(
      mat = avg_transitions$Avg_Prob_Matrix,
      facet_id = avg_transitions$facet_id,
      edge_weights = "prob"
    ),
    build_graph_data,
    label_order = syn_order,
    min_prob = min_prob
  )
  
  all_nodes <- bind_rows(purrr::map(graph_data, "nodes"))
  all_edges <- bind_rows(purrr::map(graph_data, "edges"))
  max_prob <- max(all_edges$prob)
  tg <- tbl_graph(nodes = all_nodes, edges = all_edges, directed = TRUE)
  
  #' Plotting code ----
  title_map <- setNames(avg_transitions$facet_label, avg_transitions$facet_id)
  cap_size <- circle(2, 'mm')
  title <- paste0(entrate_all$Expt[1], "_Syntax_ByCondition_minIBI", min_IBI)
  
  ggraph(tg, layout = "manual", x = x, y = y) +
    geom_edge_arc(
      aes(edge_width = prob, edge_colour = prob),
      arrow = arrow(length = unit(1, "mm"), type = "closed"),
      end_cap = cap_size,
      start_cap = cap_size,
      lineend = "round",
      strength = tg %>% activate(edges) %>% filter(!is_loop) %>% pull(strength),
      alpha = 0.7
    ) +
    geom_edge_loop(
      aes(edge_width = prob, edge_colour = prob, direction = direction, span = 85, strength =1.1),
      arrow = arrow(length = unit(1, "mm"), type = "closed"),
      end_cap = cap_size,
      start_cap = cap_size,
      alpha = 0.7
    ) +
    geom_node_point(aes(fill = name), size = 2, shape = 21, alpha = 0.8, show.legend = F) +
    geom_node_text(aes(x = label_x, y = label_y, label = node_labeller[name], 
                       color = name, angle = angle, hjust = hjust), size = 4.6, show.legend = F) +
    scale_edge_width(range = c(0.1, 2.5), guide = "none",) +
    scale_edge_colour_gradientn(
      colours = rev(c("#fde725", "#5ec962", "#21918c", "#3b528b", "#440154")),
      limits = c(min_prob, max_prob),
      name = paste0("Transition\nProbability\n(>",min_prob, ")")
    ) +
    scale_color_manual(values = c(call_pal_dark,c("IBI" = "grey30")))+
    scale_fill_manual(values = c(call_pal_dark,c("IBI" = "grey30")))+
    facet_manual(~facet_id, design = design, labeller = as_labeller(title_map), strip.position = "bottom") +
    coord_fixed(clip = "off", xlim = c(-1.5, 1.5), ylim = c(-1.4, 1.5)) +
    theme_void() +
    labs(fill = "Syllables")+
    theme(
      strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 15),
      # legend.position = c(0.2,0.25),
      legend.title = element_text(size = 11, margin = margin(b = 15)),
      legend.text = element_text(size = 8),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.text.align = 0,
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(20, 30, 20, 30),   # top, right, bottom, left (in points)
      panel.spacing = unit(2, "lines")
    )
  
  p_path <- file.path(plot_dir, paste0("4-2_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- strsplit(design, "\n")[[1]] %>% {data.frame(nrow = length(.), ncol = max(nchar(.)))}
  ggsave(p_path, width = 5+10*layout_df$ncol, height = 2+10*layout_df$nrow, 
         dpi = 300, bg = "white", units = "cm")
}

# 4.3 Analysis: Transitions Heatmap by Condition ============================
#' Compare transition probabilities by condition
{
  # Define label mappings
  heatmap_labeller <- c(
    "IBI" = "Start/End",
    "Flat" = "Flat",
    "Downward" = "Down",
    "Upward" = "Up",
    "Composite" = "Cmpst",
    "Chevron" = "Chev",
    "Reverse Chevron" = "Rev Chev",
    "1 Frequency Step" = "1 FreqStep",
    "2 Frequency Step" = "2 FreqStep",
    "3+ Frequency Step" = "3+ FreqStep",
    "Complex" = "Cmplx",
    "Harmonics" = "Hrm",
    "Short" = "Short",
    "Unstructured" = "Unstr",
    "Murmur" = "Murm",
    "Noisy" = "Nois"
  )
  
  freq_order <- filt_df %>% filter(Group_ID == "CC") %>% group_by(Mouse_ID, Label) %>% summarize(tot_calls = n()) %>% 
    group_by(Label) %>% summarize(mean = mean(tot_calls)) %>% pull(Label) %>% as.character() %>%  c("IBI",.)
  heatmap_order <- heatmap_labeller[freq_order]
  
  # Combine matrices to long format
  heatmap_df <- map2_dfr(
    avg_transitions$Avg_Prob_Matrix,
    avg_transitions$facet_id,
    ~ as.data.frame(as.table(.x)) %>%
      setNames(c("from","to","prob")) %>%
      mutate(facet_id = .y)
  )
  
  # Factor ordering
  heatmap_df <- heatmap_df %>%
    mutate(
      from = recode_factor(from, !!!heatmap_order),
      to   = recode_factor(to,   !!!heatmap_order)
    )
  
  title <- paste0(entrate_all$Expt[1], "_SyntaxTransitions_ByCondition_minIBI", min_IBI)
  
  ggplot(heatmap_df, aes(x = to, y = from, fill = prob)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colours = c("white", "#88CCEE", "#332288", "#CC6677", "#882255"),
      limits = c(0, 1),
      name = "Transition\nProbability"
    ) +
    facet_manual(~ facet_id, design = design, labeller = as_labeller(title_map)) + # Use your custom title_map!
    theme_minimal(base_size=14) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    labs(x = "To Syllable", y = "From Syllable",
         title = title)
  
  p_path <- file.path(plot_dir, paste0("4-3_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- strsplit(design, "\n")[[1]] %>% {data.frame(nrow = length(.), ncol = max(nchar(.)))}
  ggsave(p_path, width = 5+10*layout_df$ncol, height = 2+10*layout_df$nrow, 
         dpi = 300, bg = "white", units = "cm")
}

# 5.1 Analysis: Difference in Individual Transitions by Condition ============================
# Run LMEM stats for Transition Prob by Group ---
{
  # Define label mappings
  heatmap_labeller <- c(
    "IBI" = "Start/End",
    "Flat" = "Flat",
    "Downward" = "Down",
    "Upward" = "Up",
    "Composite" = "Cmpst",
    "Chevron" = "Chev",
    "Reverse Chevron" = "Rev Chev",
    "1 Frequency Step" = "1 FreqStep",
    "2 Frequency Step" = "2 FreqStep",
    "3+ Frequency Step" = "3+ FreqStep",
    "Complex" = "Cmplx",
    "Harmonics" = "Hrm",
    "Short" = "Short",
    "Unstructured" = "Unstr",
    "Murmur" = "Murm",
    "Noisy" = "Nois"
  )
  heatmap_order <- heatmap_labeller[freq_order]
  
  long_stats <- transition_tables %>%
    mutate(
      transitions = purrr::map(Prob_Matrix, ~ as.data.frame(as.table(.x)) %>%
                                 setNames(c("from", "to", "prob")))
    ) %>%
    select(Group_ID, Mouse_ID, Female_ID, transitions) %>%
    unnest(transitions) %>%
    mutate(
      from = recode_factor(from, !!!heatmap_order),
      to   = recode_factor(to,   !!!heatmap_order)
    )
  
  # All possible from-to pairs
  transition_grid <- expand.grid(
    from = syn_order,
    to = syn_order,
    stringsAsFactors = FALSE
  )
  
  # LMEM stats
  stat_df <- long_stats %>%
    group_by(from, to) %>%
    nest() %>%
    mutate(
      lmer_model = purrr::map(data, ~ tryCatch(
        lmer(prob ~ Group_ID + (1 | Female_ID), data = .x),
        error = function(e) NULL
      ))
    ) %>%
    filter(!map_lgl(lmer_model, is.null)) %>%
    mutate(
      # All pairwise contrasts, no adjustment yet
      em_contrasts = purrr::map(lmer_model, ~ {
        em <- try(emmeans(.x, "Group_ID"), silent = TRUE)
        if(inherits(em, "try-error")) return(tibble())
        as.data.frame(pairs(em, adjust = "none"))
      })
    ) %>%
    unnest(em_contrasts, keep_empty = TRUE) %>%
    separate(col = contrast, into = c("group1", "group2"), sep = " - ") %>%
    mutate(
      group1 = str_trim(group1),
      group2 = str_trim(group2)
    ) %>%
    arrange(p.value) %>% # optional
    group_by(from, to) %>%
    mutate(
      p_adj = p.adjust(p.value, method = "holm"),
      p_label = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      y.position = 1.1 + 0.02 * row_number()
    ) %>%
    ungroup()
  
  plot_df <- stat_df %>%
    mutate(signif_flag = p_label != "ns") %>% 
    mutate(
      group1 = factor(group1, levels = cond_order),
      group2 = factor(group2, levels = cond_order),
      pairs = paste(group1, group2, sep = "-")
    )
  
  # Significant transitions highlighted
  ggplot(plot_df, aes(x=to, y=from, fill=signif_flag)) +
    geom_tile(colour="grey80") +
    facet_grid(group2~group1) +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "white"), 
                      name = "Group difference") +
    theme_minimal(base_size=11) +
    theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1)) +
    labs(title="Significant Group-to-Group Differences in Transition Probabilities",
         x="To", y="From")
  
  p_path <- file.path(plot_dir, paste0("5-0_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- data.frame(nrow = 1, ncol = 1)
  ggsave(p_path, width = 4+16*layout_df$ncol, height = 2+14*layout_df$nrow, 
         dpi = 300, bg = "white", units = "cm")
}

# 5.1: Heatmap with significant transitions 'stared' ============================
{
  
  title <- paste0(entrate_all$Expt[1], "_TransHM_ByCondition_minIBI", min_IBI)
  
  sig_plot_df <- plot_df %>% 
    filter(p_label != "ns") %>%
    mutate(transition = paste(from, to, sep = " -> ")) %>% 
    rowwise() %>%
    mutate(maxProb = max(data$prob)) %>%
    group_by(transition) %>% 
    mutate(y.position = maxProb+ 0.03*row_number()) %>%
    select(-c(data,lmer_model)) %>% 
    ungroup()
  
  ggplot(plot_df, aes(x=to, y=from, fill=-estimate)) +
    geom_tile(colour="grey80") +
    facet_grid(group2~group1) +
    scale_fill_gradient2(name="Transition \nDifference\n(Top - Right)",
                         low="blue", mid="white", high="red", midpoint=0, na.value="grey90") +
    geom_text(data=sig_plot_df, aes(label=p_label), 
              color="red", size=5, fontface="bold", vjust=0.75) +
    theme_minimal(base_size=11) +
    theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1)) +
    labs(title="Significant Group-Pair Transitions",
         x="To", y="From")
  
  p_path <- file.path(plot_dir, paste0("5-1_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- data.frame(nrow = 1, ncol = 1)
  ggsave(p_path, width = 4+16*layout_df$ncol, height = 2+14*layout_df$nrow, 
         dpi = 300, bg = "white", units = "cm")
}

# 5.2 Analysis: Significant Difference in Individual Transitions by Condition ============================
# Individual Significant Transitions
{
  title <- paste0(entrate_all$Expt[1], "_SigTrans_ByCondition_minIBI", min_IBI)
  plot_bar_df <- long_stats %>%
    semi_join(sig_plot_df, by = c("from", "to")) %>% 
    mutate(transition = paste(from, to, sep = " -> "))
  
  ggplot(plot_bar_df) +
    geom_boxplot(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)),width = 0.25, alpha = 0.5) +
    # geom_violin(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)), width = 0.5, alpha = 0.3) +
    geom_jitter(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)), alpha = 0.5, size = 1) +
    # facet_grid(from ~ to) +
    facet_grid(. ~ transition) +
    stat_pvalue_manual(
      sig_plot_df,
      label = "p_label",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01, bracket.size = 0.5, size = 4
    ) +
    theme_minimal(base_size=13) +
    theme(
      axis.text.x = element_text(angle=40, hjust=1, vjust=1),
      axis.text.y = element_text(size=9)
    ) +
    labs(x = "Condition", y = "Transition probability",,
         fill = "Condition",
         title = "Significant Changes in Syllable Transitions") +
    scale_fill_manual(values = sel_pal)
  
  p_path <- file.path(plot_dir, paste0("5-2_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- data.frame(nrow = 1, ncol = plot_bar_df$transition %>% unique() %>% length())
  ggsave(p_path, width = 5+5*layout_df$ncol, height = 2+16*layout_df$nrow, 
         dpi = 300, bg = "white", units = "cm")
}

# 5.3 Analysis: All Difference in Individual Transitions by Condition ============================
# Individual Transitions differences
{
  title <- paste0(entrate_all$Expt[1], "_AllTrans_ByCondition_minIBI", min_IBI)
  all_plot_df <- plot_df %>% 
    mutate(transition = paste(from, to, sep = " -> ")) %>% 
    rowwise() %>%
    mutate(maxProb = max(data$prob)) %>%
    group_by(transition) %>% 
    mutate(y.position = maxProb+ 0.04*row_number()) %>%
    select(-c(data,lmer_model)) %>% 
    ungroup()
  
  plot_bar_df <- long_stats %>%
    semi_join(all_plot_df, by = c("from", "to")) %>% 
    mutate(transition = paste(from, to, sep = " -> "))
  
  ggplot(plot_bar_df) +
    geom_boxplot(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)),width = 0.25, alpha = 0.5) +
    # geom_violin(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)), width = 0.5, alpha = 0.3) +
    geom_jitter(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)), alpha = 0.5, size = 1) +
    facet_grid(from ~ to) +
    # facet_grid(. ~ transition) +
    stat_pvalue_manual(
      sig_plot_df,
      label = "p_label",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.03, bracket.size = 0.5, size = 4
    ) +
    theme_minimal(base_size=13) +
    theme(
      axis.text.x = element_text(angle=40, hjust=1, vjust=1),
      axis.text.y = element_text(size=9)
    ) +
    labs(x = "Condition", 
         y = "Transition probability",
         subtitle = "(Current -> Next)",
         fill = "Condition",
         title = "All Syllable Transitions") +
    scale_fill_manual(values = sel_pal)
  
  p_path <- file.path(plot_dir, paste0("5-3_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- data.frame(nrow = length(heatmap_labeller), ncol = length(heatmap_labeller))
  ggsave(p_path, width = 5+4*layout_df$ncol, height = 2+4*layout_df$nrow, 
         dpi = 300, bg = "white", units = "cm")
}


# 5.4 Analysis: Selected Difference in Individual Transitions by Condition ============================

# Top 10 mean (estimate) differences
selected_transitions <- all_plot_df %>% mutate(absEst = abs(estimate)) %>% slice_max(absEst,n=10) %>% arrange(estimate)

# print(names(heatmap_labeller)) # Printed to help selection
# selected_transitions <- tibble::tibble(
#   from = c("Unstructured", "Harmonics", "Harmonics"),
#   to   = c("Flat", "Flat", "Composite")
# )

# Individual selected Transitions
{
  title <- paste0(entrate_all$Expt[1], "_SelTrans_ByCondition_minIBI", min_IBI)
  plot_bar_df <- long_stats %>%
    semi_join(selected_transitions, by = c("from", "to")) %>% 
    mutate(transition = paste(from, to, sep = " -> "))
  
  ggplot(plot_bar_df) +
    geom_boxplot(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)),width = 0.25, alpha = 0.5) +
    # geom_violin(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)), width = 0.5, alpha = 0.3) +
    geom_jitter(aes(x = Group_ID, y = prob, fill = fct_rev(Group_ID)), alpha = 0.5, size = 1) +
    # facet_grid(from ~ to) +
    facet_grid(. ~ transition) +
    stat_pvalue_manual(
      selected_transitions,
      label = "p_label",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01, bracket.size = 0.5, size = 4
    ) +
    theme_minimal(base_size=13) +
    theme(
      axis.text.x = element_text(angle=40, hjust=1, vjust=1),
      axis.text.y = element_text(size=9)
    ) +
    labs(x = "Condition", y = "Transition probability",,
         fill = "Condition",
         title = "Largest Mean Changes in Syllable Transitions") +
    scale_fill_manual(values = sel_pal)
  
  p_path <- file.path(plot_dir, paste0("5-4_", title,".png"))
  cat("Saving", p_path, fill=TRUE)
  layout_df <- data.frame(nrow = 1, ncol = plot_bar_df$transition %>% unique() %>% length())
  ggsave(p_path, width = 5+5*layout_df$ncol, height = 2+18*layout_df$nrow, 
         dpi = 300, bg = "white", units = "cm")
}

