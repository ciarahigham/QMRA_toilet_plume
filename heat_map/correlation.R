##############################################################################
### Calculates correlation coefficients between variables                  ###
### Need to run dr_mc_uniform_tdur.R first                                 ###
### Ciara A. Higham Oct 2024 - Indoor Environments paper                   ###
##############################################################################

###############################
###       Main code         ###
###############################

# find current wd
current_wd <- dirname(rstudioapi::getActiveDocumentContext()$path)

# vector we want to find the correlation coefficients for
variables <-
  c("RH",
    "v",
    "t_dur",
    "t_enter",
    "Breath_Rate",
    "Vol_Faeces",
    "rho",
    "f")

# scenarios
arrangements <- c("S1", "S2")

# ventilation rates
ach_levels <- c(1.5, 3, 6)

# groupings for t_enter
t_enter_conditions <- list(
  "<600" = combined$t_enter < 600,
  "<60" = combined$t_enter < 60,
  "60-600" = combined$t_enter >= 60 & combined$t_enter < 600
)

# data frame for the results
results <- data.frame(
  Variable = character(),
  Arrangement = character(),
  ACH = numeric(),
  t_enter_condition = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (pathogen in c("SARS-CoV-2")) {
  for (counter_loc in c("counter_a")) {
    # wd for dose response data from dose_response_uniform_tdur.R
    # and wd for location to save file
    if (pathogen == "SARS-CoV-2") {
      data_wd <- paste0(current_wd, "/dr_output/sars_cov_2/")
      save_wd <- paste0(current_wd, "/plots/sars_cov_2/")
      file_name <- "sars_cov_2"
    } else if (pathogen == "Norovirus") {
      data_wd <- paste0(current_wd, "/dr_output/norovirus/")
      save_wd <- paste0(current_wd, "/plots/norovirus/")
      file_name <- "norovirus"
    }
    
    setwd(data_wd) # set to relevant dose-response data wd
    
    # read in the dose-response data
    if (counter_loc == "counter_a") {
      data_S1 <- read.csv("all_data_s1_counter_a.csv", header = TRUE)
      data_S2 <-
        read.csv("all_data_s2_counter_a.csv", header = TRUE)
      file_name <- paste0(file_name, "_counter_a.png")
    } else if (counter_loc == "counter_b") {
      data_S1 <- read.csv("all_data_s1_counter_b.csv", header = TRUE)
      data_S2 <-
        read.csv("all_data_s2_counter_b.csv", header = TRUE)
      file_name <- paste0(file_name, "_counter_b.png")
    }
    
    # combine all 4 sets of data
    combined <- rbind(data_S1, data_S2)
    
    # response as a percentage rather than a factor
    combined$Response <- combined$Response * 100
    
    # for each scenario
    for (arrangement in arrangements) {
      # and each ventilation rate
      for (ach in ach_levels) {
        for (t_enter_cond_name in names(t_enter_conditions)) {
          t_enter_cond <- t_enter_conditions[[t_enter_cond_name]]
          
          # find the subset of data
          subset_data <-
            combined[combined$Arrangement == arrangement &
                       combined$ACH == ach &
                       t_enter_cond,]
          
          # then find the correlation coefficient between infection risk for each variable
          # using spearman correlation ranking
          for (variable in variables) {
            correlation <-
              cor(subset_data[[variable]], subset_data$Response, method = "spearman")
            
            # add the data to results data frame
            results <- rbind(
              results,
              data.frame(
                Variable = variable,
                Arrangement = arrangement,
                ACH = ach,
                t_enter_condition = t_enter_cond_name,
                Correlation = correlation
              )
            )
          }
        }
      }
    }
  }
}

# print results
print(results)