dt <- 0.1
motiondata <-  generate_data(T = 2, frequency = 0.5,
                             amplitude = 20,
                             sensor_sd = 2.0)
# plot_trajectories(motiondata = motiondata)

motiondata_hidden <- motiondata %>%
    gather(type, value = value, -time) %>%
    filter(type != "observations") %>%
    mutate(type = factor(type,
                         levels = c("position",
                                    "velocity",
                                    "acceleration")))
motiondata_hidden <- motiondata %>%
    gather(type, value = value, -time) %>%
    # filter(type != "observations") %>%
    mutate(type = factor(type,
                         levels = c("position",
                                    "velocity",
                                    "acceleration",
                                    "observations")))
# df_tmp <- data_frame(time = motiondata$time,
#                      observations = motiondata$observations)
# df_tmp$param <- as.factor("Measurements")

p1 <- motiondata_hidden %>%
    filter(type != "observations") %>%
    ggplot(aes(x = time, y = value)) +
    geom_line(aes(linetype = type), size = 1.4) +
    geom_point(data = filter(motiondata_hidden, type == "observations"),
               aes(x = time, y = value),
               shape = 21, colour = "white",
               fill = "white", size = 6, stroke = 2) +
    geom_point(data = filter(motiondata_hidden, type == "observations"),
               aes(x = time, y = value),
               colour = "black",
               fill = "white", size = 4, stroke = 2) +
    # scale_shape_manual(name = element_blank(),
    #                    labels = "Measurements",
    #                    values = 21) +
    guides(linetype = guide_legend(keywidth = 2.5, keyheight = 1)) +
    xlab("Time (sec)") +
    ylab("Angular position (deg)") +
    facet_wrap(~ type, ncol = 1) +
    theme_apa() +
    theme(legend.title=element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour="black",
                                   size = 0.75, lineend="square"))
p1
