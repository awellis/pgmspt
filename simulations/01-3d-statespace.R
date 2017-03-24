# Rbiips model: 3-d state  ----

library(ggplot2)
library(dplyr)
library(tidyr)
library(Rbiips)

# library(viridis)
# color_palette <- viridis::plasma(n = 9)

sd2precision <- function(sd) {
    prec <- 1/(sd^2)
    prec
}

neff <- function(weights) {
    1/sum(weights^2)
}

Mode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

#
# trdata <- vestcog::generate_data(T = 2, amplitude = 20, sensor_sd = 1.7,
#                                  as_df = TRUE, seed = TRUE)

dt <- 0.1
t_max <- 2 # how long is the movement (in seconds)? #length(trdata$observations)
freq <- .5 # what is the frequency (in Hz)?

mean_a_init <- 0
prec_a_init <- sd2precision(1)
mean_w_init <- 0
prec_w_init <- sd2precision(1)
mean_p_init <- 0
prec_p_init <- sd2precision(1)

prec_a <- sd2precision(0.5)
prec_w <- sd2precision(0.5)
prec_p <- sd2precision(0.3)
prec_y <- sd2precision(1.7)

# f_eval <- function(x, k) { .5 * x + 25*x/(1+x^2) + 8*cos(1.2*k) }
# f_dim <- function(x_dim, k_dim) { 1 }
#
# f_eval_2 <- function(A, freq, t) {
#     A/2 * sin(two_pi*(t-1) * freq) + A * sin(two_pi*(t-1) * 4 * freq)
# }

data = list(t_max = t_max/dt,
            dt = dt,
            freq = dt * freq,
            A = 20,
            # y = trdata$observations,
            #a = trajectory_data$acceleration/20,
            mean_a_init = mean_a_init,
            prec_a_init = prec_a_init,
            mean_w_init = mean_w_init,
            prec_w_init = prec_w_init,
            mean_p_init = mean_p_init,
            prec_p_init = prec_p_init,
            prec_a = prec_a,
            prec_w = prec_w,
            prec_p = prec_p,
            prec_y = prec_y,
            two_pi = 2*pi)

model_string <- "
data {
  p_true[1] ~ dnorm(mean_p_init, prec_p_init)
  w_true[1] ~ dnorm(mean_w_init, prec_w_init)

  y[1] ~ dnorm(w_true[1], prec_y)

  for (t in 2:t_max) {
    # p: position
    # w: velocity
    a_true[t-1] <- A * sin(two_pi*(t-1) * freq) + A * sin(two_pi*(t-1) * 4 * freq)
    p_true[t] ~ dnorm(p_true[t-1] + dt * w_true[t-1] + dt^2/2 * a_true[t-1], prec_p)
    w_true[t] ~ dnorm(w_true[t-1] + dt * a_true[t-1], prec_w)

    y[t] ~ dnorm(w_true[t], prec_y)
  }
}

model {
  p[1] ~ dnorm(mean_p_init, prec_p_init)
  w[1] ~ dnorm(mean_w_init, prec_w_init)
  a[1] ~ dnorm(mean_a_init, prec_a_init)

  y[1] ~ dnorm(w[1], prec_y)

  for (t in 2:t_max) {
    # p: position
    # w: velocity
    # a[t-1] <- A * sin(two_pi*(t-1)/t_max)
    p[t] ~ dnorm(p[t-1] + dt * w[t-1] + dt^2/2 * 1 * a[t-1], prec_p)
    w[t] ~ dnorm(w[t-1] + dt * A * a[t-1], prec_w)
    a[t] ~ dnorm(a[t-1], prec_a)

    y[t] ~ dnorm(w[t], prec_y)
  }
}"


model = biips_model(file = textConnection(model_string), data = data,
                    sample_data = TRUE)

n_particles <- 3000 # Number of particles

variables <- c('a', 'w', 'p') # Variables to be monitored
mn_type <- 'f'
rs_type <- c('systematic', 'multinomial', 'stratified')
rs_thres <- 0.5 # n_particles # Optional parameters

# Run SMC ----
out_smc <- biips_smc_samples(model, variables, n_particles,
                             type = mn_type, rs_type = rs_type[2],
                             rs_thres = rs_thres)
summ_smc <- biips_summary(out_smc, probs=c(.025, .975))

## save model fit ----
m1_3d_statespace <- list(model = model,
                         out_smc = out_smc,
                         summ_smc = summ_smc)
save(m1_3d_statespace, file = "simulations/model-fits/m1-3d-statespace.Rda")


# vel_f_mean <- summ_smc$w$f$mean
# vel_f_lower <- summ_smc$w$f$quant$`0.025`
# vel_f_upper <- summ_smc$w$f$quant$`0.975`
#
# pos_f_mean <-  summ_smc$p$f$mean
# pos_f_lower <- summ_smc$p$f$quant$`0.025`
# pos_f_upper <- summ_smc$p$f$quant$`0.975`

# df <- with(summ_smc, {
#     data_frame(t = 1:length(w$f$mean),
#                vel_mean = w$f$mean,
#                vel_lower = w$f$quant$`0.025`,
#                vel_upper = w$f$quant$`0.975`,
#                pos_mean =  p$f$mean,
#                pos_lower = p$f$quant$`0.025`,
#                pos_upper = p$f$quant$`0.975`)
# })

df_w <- with(summ_smc, {
    data_frame(t = 1:length(w$f$mean),
               mean = w$f$mean,
               lower = w$f$quant$`0.025`,
               upper = w$f$quant$`0.975`)
})
df_w$param <- "Velocity"
df_w$true <- model$data()$w_true
df_w$observations <- model$data()$y

df_p <- with(summ_smc, {
    data_frame(t = 1:length(p$f$mean),
               mean = p$f$mean,
               lower = p$f$quant$`0.025`,
               upper = p$f$quant$`0.975`)
})
df_p$param <- "Position"
df_p$true <- model$data()$p_true
df_p$observations <- model$data()$y

df <- rbind(df_p, df_w)

# sensor <- data_frame(velocity = model$data()$w,
#                      observations = model$data()$y)

# df <- data.frame(t = 1:t_max,
#                  vel_mean = vel_f_mean,
#                  vel_lower = vel_f_lower,
#                  vel_upper = vel_f_upper,
#                  pos_mean = pos_f_mean,
#                  pos_lower = pos_f_lower,
#                  pos_upper = pos_f_upper,
#                  velocity = model$data()$w, #trdata$velocity,
#                  observations = model$data()$y)
#                  # observations = trdata$observations)


# df_tidy <- df %>%
#     gather(param, value, -t)
#
# plot_filtering_estimates <- function(df) {
#     p <- ggplot(data = df, aes(x = t)) +
#
#         geom_ribbon(aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.3,
#                     fill = color_palette[2]) +
#         geom_line(aes(y = velocity), colour = color_palette[7], alpha = 0.9,
#                   linetype = 2, size = 1.2) +
#
#         geom_line(aes(y = vel_mean), colour = color_palette[2], size = 1.6) +
#
#         geom_point(aes(y = observations), colour = color_palette[1],
#                    size = 2, shape = 15, alpha = 0.6) +
#         theme_bw() +
#         theme(legend.title = element_blank()) +
#         ylab("") + xlab("")
#     # ggtitle("Velocity storage")
#
#     # geom_vline(xintercept=as.numeric(as.Date("1959-12-01")), linetype=2) +
#     # ggtitle(paste0("ARIMA -- Holdout MAPE = ", round(100*MAPE,2), "%")) +
#     # theme(axis.text.x=element_text(angle = -90, hjust = 0))
#
#     print(p)
# }
#
#
# plot_filtering_estimates <- function(df, pos = TRUE) {
#     p <- ggplot(data = df, aes(x = t)) +
#         geom_ribbon(aes(ymin = vel_lower, ymax = vel_upper), alpha = 0.3,
#                     fill = color_palette[2]) +
#         geom_line(aes(y = velocity), colour = color_palette[7], alpha = 0.9,
#                   linetype = 2, size = 1.2) +
#
#         geom_line(aes(y = vel_mean), colour = color_palette[2], size = 1.6) +
#
#         geom_point(aes(y = observations), colour = color_palette[1],
#                    size = 2, shape = 15, alpha = 0.6) +
#         theme_bw() +
#         theme(legend.title = element_blank()) +
#         ylab("") + xlab("")
#     # ggtitle("Velocity storage")
#
#     # geom_vline(xintercept=as.numeric(as.Date("1959-12-01")), linetype=2) +
#     # ggtitle(paste0("ARIMA -- Holdout MAPE = ", round(100*MAPE,2), "%")) +
#     # theme(axis.text.x=element_text(angle = -90, hjust = 0))
#
#     if (pos) {
#         p <- p + geom_ribbon(aes(ymin = pos_lower, ymax = pos_upper), alpha = 0.3,
#                              fill = color_palette[6]) +
#             geom_line(aes(y = pos_mean), colour = color_palette[6], size = 1.6)
#     }
#
#     print(p)
# }


# plot_filtering_estimates(df, pos = TRUE)

p1 <- df %>%
    ggplot(data = ., aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey90") +

    geom_line(aes(y = true), alpha = 0.9,
              linetype = 4, size = 1.2) +

    geom_line(aes(y = mean), size = 1.6) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    facet_wrap(~ param, nrow = 1) +
    scale_x_continuous(limits = c(0, 20),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "0.5", "1", "1.5", "2")) +
    xlab("Time (sec)") + ylab("Angular position") +
    papaja::theme_apa(base_size = 18)

p2 <- p1 + geom_point(data = filter(df, param == "velocity"),
                      aes(y = observations),
                      alpha = 1, fill = "white",
                      colour = "white",
                      shape = 21, size = 6) +
    geom_point(data = filter(df, param == "velocity"),
               aes(y = observations), alpha = 1,
               fill = "black",
               colour = "grey40", shape = 21, size = 4)

p2





p3 <- df %>%
    filter(param == "Position") %>%
    ggplot(data = ., aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey90") +

    geom_line(aes(y = true), alpha = 0.9,
              linetype = 4, size = 1.2) +

    geom_line(aes(y = mean), size = 1.6) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    scale_x_continuous(limits = c(0, 20),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "0.5", "1", "1.5", "2")) +
    xlab("Time (sec)") + ylab("Angular position (deg)") +
    papaja::theme_apa(base_size = 18)
p3

p4 <- df %>%
    filter(param == "Velocity") %>%
    ggplot(data = ., aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey90") +

    geom_line(aes(y = true), alpha = 0.9,
              linetype = 4, size = 1.2) +

    geom_line(aes(y = mean), size = 1.6) +
    geom_point(aes(y = observations),
               alpha = 1, fill = "white",
               colour = "white",
               shape = 21, size = 6) +
    geom_point(aes(y = observations), alpha = 1,
               fill = "black",
               colour = "grey40", shape = 21, size = 4) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    scale_x_continuous(limits = c(0, 20),
                       breaks = c(0, 5, 10, 15, 20),
                       labels = c("0", "0.5", "1", "1.5", "2")) +
    xlab("Time (sec)") + ylab("Angular Velocity (deg/s)") +
    papaja::theme_apa(base_size = 18)
p4

p5 <- cowplot::plot_grid(p3, p4)
p5
# cowplot::plot_grid(p5, p5, ncol = 1)
