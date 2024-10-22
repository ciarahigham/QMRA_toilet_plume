##############################################################################
### Plots for the infection risk for counter A and B                       ###
### Plots infection risk for SARS-CoV-2 and norovirus                      ###
### i.e. generates 4 plots with violins for male/female                    ###
### Need to run dr_mc_male_female.R first                                  ###
###                                                                        ###
### Ciara A. Higham Oct 2024 - Indoor Environments paper                   ###
##############################################################################

###############################
###        Pre-lims         ###
###############################

#install.packages("ggh4x")
#install.packages("ggplot2")
#install.packages("extrafont")

library(ggh4x)
library(ggplot2)
library(extrafont)

loadfonts()

##############################################################################################################
##############################################################################################################

# The following blocks of code extend ggplot2 GeomViolin and make adjustments to add quantiles to the plot
# see stackoverflow (https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2) (@jan-glx, @Axeman)
# and stackoverflow (https://stackoverflow.com/questions/47651868/split-violin-plot-with-ggplot2-with-quantiles) (@Axeman)

GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <-
      transform(
        data,
        xminv = x - violinwidth * (x - xmin),
        xmaxv = x + violinwidth * (xmax - x)
      )
    grp <- data[1, "group"]
    newdata <-
      plyr::arrange(transform(data, x = if (grp %% 2 == 1)
        xminv
        else
          xmaxv), if (grp %% 2 == 1)
            y
        else
          - y)
    newdata <-
      rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
      round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 &
        !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <-
        create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
      aesthetics <-
        data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <-
        rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <-
        GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

create_quantile_segment_frame <-
  function(data,
           draw_quantiles,
           split = FALSE,
           grp = NULL) {
    dens <- cumsum(data$density) / sum(data$density)
    ecdf <- stats::approxfun(dens, data$y)
    ys <- ecdf(draw_quantiles)
    violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
    violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
    violin.xs <- (stats::approxfun(data$y, data$x))(ys)
    if (grp %% 2 == 0) {
      data.frame(
        x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
        y = rep(ys, each = 2),
        group = rep(ys, each = 2)
      )
    } else {
      data.frame(
        x = ggplot2:::interleave(violin.xminvs, violin.xs),
        y = rep(ys, each = 2),
        group = rep(ys, each = 2)
      )
    }
  }

geom_split_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "identity",
           ...,
           draw_quantiles = NULL,
           trim = TRUE,
           scale = "area",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomSplitViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        draw_quantiles = draw_quantiles,
        na.rm = na.rm,
        ...
      )
    )
  }

##############################################################################################################
##############################################################################################################

# comparison times
t1 <- "0"
t2 <- "60"
t3 <- "240"

# find current working directory
current_wd <-
  dirname(rstudioapi::getActiveDocumentContext()$path) # find current wd


# list of the new strip labels
strip_labels_new <-
  c("1.5 ACH" ,
    "3 ACH" ,
    "6 ACH" ,
    paste(t1, "s"),
    paste(t2, "s"),
    paste(t3, "s"),
    "S2",
    "S1")

# list of the current labels
strip_labels_current <-
  c("1.5", "3", "6", t1, t2, t3, "S2", "S1")

# set the strip_labels to the new strip labels
strip_labels <- setNames(strip_labels_new, strip_labels_current)

###############################
###       Main code         ###
###############################

for (pathogen in c("Norovirus", "SARS-CoV-2")) {
  for (counter_loc in c("counter_a")) {
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
        max_sars <- max(filtered_combined$Response * 100)
        min_sars <- min(filtered_combined$Response * 100)
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
        # the maximum dose response value for counter a for norovirus
        max_counter_a <- max(filtered_combined$Response * 100)
        max_counter_a_formatted <- signif(max_counter_a, 3)
        max_norovirus <- max(filtered_combined$Response * 100)
        min_norovirus <- min(filtered_combined$Response * 100)
        y_label <- bquote("Normalised risk")
        
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
    
    
    filtered_combined$Response <- filtered_combined$Response * 100
    
    # calculate a normalised response value by dividing by max of counter A norovirus
    filtered_combined$Normalised_Response <-
      filtered_combined$Response / max_counter_a
    
    # set the colour for female and male
    filtered_combined$fill_color <-
      ifelse(filtered_combined$Gender == "Female", "#009ADE", "#FF1F5B")
    
    # create the plot
    p <- ggplot() +
      
      # create the split violin plot
      geom_split_violin(
        data = filtered_combined,
        
        # select the data to plot
        aes(x = Arrangement, y = Normalised_Response, fill = fill_color),
        
        # 25th, 50th and 75th quantile
        draw_quantiles = c(0.25, 0.5, 0.75),
        
        
        alpha = 0.8,
        size = 0.15
      ) +
      
      # nest the facets it is split by ACH, then split by t_enter, then split by S1/S2
      # and add the labels
      facet_nested(
        . ~ ACH + t_enter + Arrangement,
        scales = "free",
        labeller = as_labeller(strip_labels)
      ) +
      
      # plot on a log scale
      scale_y_log10(
        limits = c(1e-09, 1),
        breaks = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1),
        labels = c(
          expression(10 ^ {
            -10
          }) ,
          expression(10 ^ {
            -8
          }),
          expression(10 ^ {
            -6
          }),
          expression(10 ^ {
            -4
          }),
          expression(10 ^ {
            -2
          }),
          "1"
        )
      ) +
      
      theme_bw() +
      
      # colours for the violin
      scale_fill_manual(
        values = c("#009ADE", "#FF1F5B"),
        labels = c("Female", "Male"),
        name = NULL,
        guide = "legend"
      ) +
      scale_x_discrete(breaks = NULL) + # remove x breaks
      
      # aesthetics for the plot text size etc
      theme(
        text = element_text(family = "Arial"),
        # Arial font
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 8, margin = margin(2.65, 0.5, 2.65, 0.5)),
        strip.background = element_rect(color = "black", size = 0.2),
        legend.position = "top",
        legend.text = element_text(size = 8),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, family = "Arial"),
        axis.text.x = element_blank(),
        legend.margin = margin(0.01, 0.1, 0.01, 0.1),
        legend.key.size = unit(0.3, 'cm'),
      ) +
      
      # y label
      ylab(y_label)
    
    # print plot
    print(p)
    
    setwd(save_wd) # set save wd
    
    # save plot
    ggsave(
      file_name,
      plot = p,
      width = 19,
      height = 6,
      units = "cm",
      dpi = 500
    )
  }
}