##############################################################################
### Plots for data of dr_mc_uniform_tdur.R                                 ###
### Plots t_dur as a function on t_enter and normalised risk               ###
### Particle counter location A                                            ###
### Normalised risk is the risk divided by median risk across an ACH and   ###
### scenario for a given virus type                                        ###
### Ciara A. Higham Oct 2024 - Indoor Environments paper                   ###
##############################################################################

###############################
###        Pre-lims         ###
###############################

# # install packages
#
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("scales")
# install.packages("extrafont")

library(ggplot2)
library(dplyr)
library(scales)
library(extrafont)

loadfonts()

###############################
###    Input parameters     ###
###############################

# add ACH to strip labels
strip_labels_new <- c("1.5 ACH" , "3 ACH" , "6 ACH")
strip_labels_current <- c("1.5", "3", "6")
strip_labels <- setNames(strip_labels_new, strip_labels_current)

# counter location A
counter_loc <- "counter_a"

# find current wd
current_wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
colours = c("#0571b0",
            "#92c5de",
            "#f7f7f7",
            "#f4a582",
            "#ca0020")
for (pathogen in c("SARS-CoV-2", "Norovirus")) {
  # each pathogen
  
  for (arrangement in c("S1", "S2")) {
    # wd for dose response data from all_times_dose_response_monte_carlo.R
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
    
    # set to relevant dose-response data wd
    setwd(data_wd)
    
    # read in the dose-response data
    # and set the save file name
    if (arrangement == "S1") {
      dose_response_df <-
        read.csv("all_data_s1_counter_a.csv", header = TRUE)
      file_name <- paste0(file_name, "_S1_counter_a.png")
      median_response = median(dose_response_df$Response)
      max_response = max(dose_response_df$Response) / median_response
      min_response = min(dose_response_df$Response[dose_response_df$Response > 0]) / median_response
      log_max_response = log10(max_response)
      log_min_response = log10(min_response)
      legend_name <- bquote("Normalised risk")
    } else if (arrangement == "S2") {
      dose_response_df <-
        read.csv("all_data_s2_counter_a.csv", header = TRUE)
      file_name <- paste0(file_name, "_S2_counter_a.png")
    }
    
    # calculate normalised reponse and create a column for when normalised is 0
    dose_response_df <- dose_response_df %>%
      group_by(ACH) %>%
      mutate(Normalised = Response / median_response) %>%
      ungroup() %>%
      mutate(IsZero = Normalised == 0)
    
    # order so the maximum normalised value is first
    dose_response_df <-
      dose_response_df[order(dose_response_df$Normalised), ]
    
    # define the legend values ensuring 1 (median) is in the middle
    legend_values <-
      c(min_response, 10 ^ ((log_min_response + log(1)) / 2), 1, 10 ^ ((log10(1) +
                                                                          log_max_response) / 2), max_response)
    
    # transform the values to the log scale
    log_legend_values <- log10(legend_values)
    
    # rescale the log-transformed legend values to the range [0, 1]
    rescaled_log_values <-
      rescale(log_legend_values,
              to = c(0, 1),
              from = range(log_legend_values))
    
    # create the plot
    plot <-
      ggplot(dose_response_df,
             aes(x = t_enter, y = t_dur, fill = Normalised)) +
      
      # yellow points for Normalised = 0, as they don't go on the log scale
      geom_point(
        data = subset(dose_response_df, IsZero == TRUE),
        aes(x = t_enter, y = t_dur, color = "0"),
        size = 0.01,
        shape = '.'
      ) +
      # generate the heat map
      geom_raster(
        data = subset(dose_response_df, IsZero == FALSE),
        aes(x = t_enter, y = t_dur, fill = Normalised),
        interpolate = FALSE
      ) +
      
      # colour gradient for the normalised rsk
      scale_fill_gradientn(
        colors = colours,
        values = rescaled_log_values,
        limits = c(min_response, max_response),
        breaks = c(
          1e-7,
          1e-6,
          1e-5,
          1e-4,
          1e-3,
          1e-2,
          1e-1,
          1e0,
          1e1,
          1e2,
          1e3,
          1e4,
          1e5
        ),
        labels = c(
          1e-7,
          1e-6,
          1e-5,
          1e-4,
          1e-3,
          1e-2,
          1e-1,
          1e0,
          1e1,
          1e2,
          1e3,
          1e4,
          1e5
        ),
        name = legend_name,
        trans = "log10"
      ) +
      
      # axis labels
      labs(x = expression("t"[enter] ~ "(min)"), y = expression("t"[dur] ~
                                                                  "(min)")) +  # Adding (minutes) to axis labels
      
      # plot each ACH as an individual facet
      facet_wrap(~ ACH, labeller = as_labeller(strip_labels)) +
      
      # change x and y labels to minutes
      scale_x_continuous(
        limits = c(0, 600),
        breaks = seq(0, 1200, by = 60),
        labels = seq(0, 1200, by = 60) / 60,
        minor_breaks = seq(0, 1200, by = 30)
      ) +
      scale_y_continuous(
        limits = c(0, 900),
        breaks = seq(0, 900, by = 120),
        labels = seq(0, 900, by = 120) / 60,
        minor_breaks = seq(0, 900, by = 60)
      ) +
      
      # aesthetics e.g. text size
      theme_bw() +
      theme(
        axis.text = element_text(size = 8, family = "Arial"),
        strip.text = element_text(
          size = 8,
          margin = margin(2.65, 0.5, 2.65, 0.5),
          family = "Arial"
        ),
        axis.title.x = element_text(size = 8, family = "Arial"),
        # x-axis title
        axis.title.y = element_text(size = 8, family = "Arial"),
        # y-axis title
        strip.background = element_rect(color = "black", linewidth = 0.2),
        legend.position = "top",
        legend.text = element_text(size = 8, family = "Arial"),
        legend.title = element_text(size = 8, family = "Arial"),
        panel.border = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(2, "cm"),
        legend.margin = margin(0.01, 0.1, 0.01, 0.1)
      ) +
      
      
      guides(
        color = guide_legend(
          title = " ",
          # no title
          title.position = "top",
          label.position = "bottom",
          title.hjust = 3,
          override.aes = list(
            size = 4.5,
            shape = 22,
            fill = "yellow" ,
            colour = "black"
          ),
          # Reduce point size in legend
          keywidth = unit(0.5, "cm"),
          keyheight = unit(0.3, "cm"),
          order = 1
        ),
        fill = guide_colourbar(
          title.position = "top",
          title.hjust = 0.5,
          order = 2
        )
      )  +
      scale_color_manual(values = c("0" = "yellow"))
    # set save wd
    setwd(save_wd)
    
    # save the plot
    ggsave(
      paste0(
        tolower(pathogen),
        "_",
        arrangement,
        "_",
        counter_loc,
        ".png"
      ),
      plot = plot,
      width = 19,
      height = 9,
      units = "cm",
      dpi = 500
    )
  }
}
