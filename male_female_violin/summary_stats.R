##############################################################################
### Generates median and mean values found in Supplementary Material       ###
###                                                                        ###
### Ciara A. Higham Oct 2024 - Indoor Environments paper                   ###
##############################################################################

#install.packages("dplyr")
library(dplyr)

# comparison times
t1 <- "0"
t2 <- "60"
t3 <- "240"

# find current working directory
current_wd <-
  dirname(rstudioapi::getActiveDocumentContext()$path) # find current wd

# Initialize empty data frames for summary statistics
summary_stats_sars_cov_2 <- data.frame()
summary_stats_norovirus <- data.frame()

# Loop over each pathogen
for (pathogen in c("SARS-CoV-2", "Norovirus")) {
  for (counter_loc in c("counter_a", "counter_b")) {
    # wd for dose response data from dr_mc_male_female.R
    # and wd for location to save file
    if (pathogen == "SARS-CoV-2") {
      data_wd <- paste0(current_wd, "/dr_output/sars_cov_2/")
      save_wd <- paste0(current_wd, "/plots/sars_cov_2/")
      file_name <- "sars_cov_2"
      setwd(data_wd) # set to relevant dose-response data wd
      if (counter_loc == "counter_a") {
        male_S1 <- read.csv("male_s1_counter_a.csv", header = TRUE)
        male_S2 <- read.csv("male_s2_counter_a.csv", header = TRUE)
        female_S1 <-
          read.csv("female_s1_counter_a.csv", header = TRUE)
        female_S2 <-
          read.csv("female_s2_counter_a.csv", header = TRUE)
        file_name <- paste0(file_name, "_counter_a.png")
        combined <- rbind(female_S1, female_S2, male_S1, male_S2)
        filtered_combined <-
          subset(combined,
                 t_enter %in% c(as.numeric(t1), as.numeric(t2), as.numeric(t3)))
      } else if (counter_loc == "counter_b") {
        male_S1 <- read.csv("male_s1_counter_b.csv", header = TRUE)
        male_S2 <- read.csv("male_s2_counter_b.csv", header = TRUE)
        female_S1 <-
          read.csv("female_s1_counter_b.csv", header = TRUE)
        female_S2 <-
          read.csv("female_s2_counter_b.csv", header = TRUE)
        file_name <- paste0(file_name, "_counter_b.png")
        combined <- rbind(female_S1, female_S2, male_S1, male_S2)
        filtered_combined <-
          subset(combined,
                 t_enter %in% c(as.numeric(t1), as.numeric(t2), as.numeric(t3)))
      }
    } else if (pathogen == "Norovirus") {
      data_wd <- paste0(current_wd, "/dr_output/norovirus/")
      save_wd <- paste0(current_wd, "/plots/norovirus/")
      file_name <- "norovirus"
      setwd(data_wd) # set to relevant dose-response data wd
      if (counter_loc == "counter_a") {
        male_S1 <- read.csv("male_s1_counter_a.csv", header = TRUE)
        male_S2 <- read.csv("male_s2_counter_a.csv", header = TRUE)
        female_S1 <-
          read.csv("female_s1_counter_a.csv", header = TRUE)
        female_S2 <-
          read.csv("female_s2_counter_a.csv", header = TRUE)
        file_name <- paste0(file_name, "_counter_a.png")
        combined <- rbind(female_S1, female_S2, male_S1, male_S2)
        filtered_combined <-
          subset(combined,
                 t_enter %in% c(as.numeric(t1), as.numeric(t2), as.numeric(t3)))
      } else if (counter_loc == "counter_b") {
        male_S1 <- read.csv("male_s1_counter_b.csv", header = TRUE)
        male_S2 <- read.csv("male_s2_counter_b.csv", header = TRUE)
        female_S1 <-
          read.csv("female_s1_counter_b.csv", header = TRUE)
        female_S2 <-
          read.csv("female_s2_counter_b.csv", header = TRUE)
        file_name <- paste0(file_name, "_counter_b.png")
        combined <- rbind(female_S1, female_S2, male_S1, male_S2)
        filtered_combined <-
          subset(combined,
                 t_enter %in% c(as.numeric(t1), as.numeric(t2), as.numeric(t3)))
      }
    }
    
    # Compute summary statistics
    summary_stats <- filtered_combined %>%
      group_by(ACH, t_enter, Arrangement) %>%
      summarize(
        median_response = formatC(signif(median(
          Response, na.rm = TRUE
        ) * 100, 2), format = "e"),
        mean_response = formatC(signif(mean(
          Response, na.rm = TRUE
        ) * 100, 2), format = "e")
      )
    
    # Add summary statistics to the appropriate data frame
    if (pathogen == "SARS-CoV-2") {
      summary_stats_sars_cov_2 <-
        rbind(summary_stats_sars_cov_2, summary_stats)
    } else if (pathogen == "Norovirus") {
      summary_stats_norovirus <-
        rbind(summary_stats_norovirus, summary_stats)
    }
    
    # Print summary statistics for each pathogen and counter location
    print(paste("Summary statistics for", pathogen, "at", counter_loc))
    print(summary_stats)
  }
}

# Display final summary statistics data frames
print("Summary statistics for SARS-CoV-2:")
print(summary_stats_sars_cov_2)

print("Summary statistics for Norovirus:")
print(summary_stats_norovirus)
