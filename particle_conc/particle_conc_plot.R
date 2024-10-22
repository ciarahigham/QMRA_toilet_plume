##############################################################################
### Plots for the particle concentration for counter A and B               ###
### With S1 and S2 on same plot                                            ###
### Ciara A. Higham Oct 2024 - Indoor Environments                         ###
##############################################################################

###############################
###        Pre-lims         ###
###############################

#install.packages("ggplot2")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("reshape2")
#install.packages("ggh4x")

library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggh4x)

###############################
###     Input parameters    ###
###############################

# colors for the S2 vs S1
colors <- c("#af58ba", "#00cd6c")

# labels for the particle sizes
strip_labels <-
  c(
    "0.3" = "0.3 micron",
    "0.5" = "0.5 micron",
    "1" = "1 micron",
    "3" = "3 micron",
    "5" = "5 micron",
    "10" = "10 micron"
  )

#######################################
### Read in the experimental data #####
#######################################

# define file names and ACH values
file_names <- c(
  "1.5ach-a.txt",
  "1.5ach-b.txt",
  "1.5ach-c.txt",
  "3ach-a.txt",
  "3ach-b.txt",
  "3ach-c.txt",
  "6ach-a.txt",
  "6ach-b.txt",
  "6ach-c.txt"
)

# ach values
ach_values <- c(1.5, 3, 6)

#current wd
current_wd <- dirname(rstudioapi::getActiveDocumentContext()$path)

inc = 5e3 # the amount by which the data is shifted - allows two scales
limits = c(inc, 2e+7) # limits of plot
breaks = c(inc, 1e+5 + inc, 1e6 + inc, 1e7 + inc) # divisions on y axis
labels = c("0", expression(10 ^ 5), expression(10 ^ 6), expression(10 ^
                                                                     7))

conc_dist_df <- data.frame() # data frame for c_{i,j,t}

for (counter_loc in c("counter_a", "counter_b")) {
  # do for each counter
  
  # plot both S2 and S1 scenario on same plot
  for (arr in c("S2", "S1")) {
    # set the wd to where the data is located
    if (arr == "S2") {
      wd <- paste0(current_wd, "/../raw_data/S2")
    } else if (arr == "S1") {
      wd <- paste0(current_wd, "/../raw_data/S1")
    }
    
    if (counter_loc == "counter_a") {
      wd <- paste0(wd, "/counter_a")
      save_wd <- paste0(current_wd, "/plots/counter_a")
    } else if (counter_loc == "counter_b") {
      wd <- paste0(wd, "/counter_b")
      save_wd <- paste0(current_wd, "/plots/counter_b")
    }
    
    setwd(wd)
    
    # read in file
    file_data <- lapply(file_names, read.table, header = TRUE)
    
    for (ach in ach_values) {
      # iterate over 1.5, 3 and 6 ACH
      
      temp_post_flush <-
        data.frame() # temp list to store the post-flush data frames
      
      # read in data file
      ach_files <-
        file_names[grep(paste0(ach, "ach-"), file_names)] # filter file to the corresponding ach in loop
      
      for (i in 1:3) {
        # loop through the replicates for each ach
        
        data <- file_data[[which(file_names == ach_files[i])]]
        
        pre_flush_rows <- data[892:921,] # data 5 minutes pre-flush
        post_flush_rows <- data[922:nrow(data),] # data post-flush
        
        colnames(post_flush_rows) <-
          c("0.3", "0.5", "1", "3", "5", "10") # column names to particle sizes
        
        bckg_conc_df <-
          apply(pre_flush_rows, 2, median) # assign median to data frame
        
        for (j in 1:6) {
          # remove background for each size
          post_flush_rows[, j] <-
            post_flush_rows[, j] - bckg_conc_df[j]
        }
        
        post_flush_rows[post_flush_rows <= 0] <- 0 # if < 0 set to 0
        post_flush_rows$Time <-
          seq(1, nrow(post_flush_rows)) * 10 # time post-flush column
        
        max_time <-
          max(post_flush_rows$Time) # max time of post-flush data i.e. when counter stops running
        
        if (max_time > 1500) {
          # if max time is greater than 25 mins
          valid_rows <-
            post_flush_rows$Time <= 1500 # just keep data up to 25 mins
          post_flush_rows <- post_flush_rows[valid_rows,]
        } else if (max_time < 1500) {
          # if the counter was stopped before 25 mins
          new_times <- seq(max_time + 10, 1500, by = 10) #
          new_rows <- data.frame(Time = new_times)
          new_rows[, c("0.3", "0.5", "1", "3", "5", "10")] <-
            0 # set the values after the max time but before 25 min to zero
          post_flush_rows <- rbind(post_flush_rows, new_rows)
        }
        
        post_flush_rows$ACH <-
          rep(ach, nrow(post_flush_rows)) # add ACH column
        post_flush_rows$Replicate <-
          rep(i, nrow(post_flush_rows)) # add a column indicating replicate number
        
        # melt so that there is a column for size
        post_flush_rows <-
          melt(post_flush_rows, id.vars = c("Replicate", "Time", "ACH"))
        
        colnames(post_flush_rows)[4] <- "Size" # name size col
        colnames(post_flush_rows)[5] <- "Conc" # name conc col
        post_flush_rows$Arrangment <- arr # column for S2/S1 scenario
        
        temp_post_flush <-
          bind_rows(temp_post_flush, post_flush_rows) # combine data for each replicate
      }
      conc_dist_df <-
        bind_rows(conc_dist_df, temp_post_flush) # combine each ACH frame
    }
  }
  
  # calculate mean and standard error across replicates
  plot_df <- conc_dist_df %>%
    group_by(Time, ACH, Size, Arrangment) %>%
    summarise(Mean_Conc = mean(Conc) + inc,
              # add inc to readjust the scale
              SE_Mean_Conc = sd(Conc) / sqrt(n()))
  
  setwd(save_wd)
  
  plots <-
    lapply(unique(plot_df$ACH), function(ach_value) {
      # plot for each ACH
      
      subset_data <-
        subset(plot_df, ACH == ach_value) # select data for relevant ACH
      
      p <-
        ggplot(data = subset(plot_df, ACH == ach_value),
               aes(x = Time, y = Mean_Conc, group = Arrangment)) +
        
        # plot standard error as a ribbon
        geom_ribbon(
          aes(
            ymin = Mean_Conc - SE_Mean_Conc,
            ymax = Mean_Conc + SE_Mean_Conc,
            fill = factor(Arrangment)
          ),
          alpha = 0.4,
          show.legend = FALSE
        ) +
        
        # plot mean_conc as a line for both arrangements
        geom_line(aes(color = factor(Arrangment)), size = 0.25)  +
        
        # colours for the plot
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        
        # change x scale between 0-10 mins and change to mins
        scale_x_continuous(
          limits = c(0, 600),
          breaks = seq(0, 1500, by = 60),
          labels = seq(0, 1500, by = 60) / 60,
          minor_breaks = seq(0, 1500, by = 30)
        ) +
        
        # add labels
        labs(x = "Time post flush (min)",
             y = "Concentration (#/m³)",
             fill = "Particle Size") +
        theme_bw() +
        
        # aesthetics eg text size etc
        theme(
          text = element_text(family = "Arial"),
          axis.text = element_text(size = 8),
          strip.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          legend.position = "top",
          legend.text = element_text(size = 8),
          plot.title = element_blank(),
          legend.title = element_blank(),
          panel.border = element_blank(),
          legend.margin = margin(c(0.001, 0.001, 0.001, 0.001)),
          legend.key.size = unit(0.4, 'cm')
        ) +
        
        # facets for each particle size
        facet_wrap( ~ factor(Size),
                    scales = "free",
                    labeller = as_labeller(strip_labels)) +
        
        # all on the same scale and remove the scales for 0.5, 1, 5, 10 micron sizes
        facetted_pos_scales(
          y = list(
            scale_y_log10(
              limits = limits,
              breaks = breaks,
              labels = labels
            ),
            scale_y_log10(
              limits = limits,
              breaks = breaks,
              guide = "none"
            ),
            scale_y_log10(
              limits = limits,
              breaks = breaks,
              guide = "none"
            ),
            scale_y_log10(
              limits = limits,
              breaks = breaks,
              labels = labels
            ),
            scale_y_log10(
              limits = limits,
              breaks = breaks,
              guide = "none"
            ),
            scale_y_log10(
              limits = limits,
              breaks = breaks,
              guide = "none"
            )
          )
        ) +
        
        guides(colour = guide_legend(override.aes = list(
          linewidth = 3, linetype = 1
        )))
      
      #save plot
      ggsave(
        paste0(ach_value, "_", counter_loc, ".png"),
        plot = p,
        width = 19,
        height = 6.9,
        units = "cm",
        dpi = 500
      )
      print(p)
    })
}