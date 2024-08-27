
# methods (1) status quo, (2) status quo using only the WGOA, (3) AR1, (4) linear interpolation
library(tidyverse)

sable <- read_csv('data/sablearearpns.csv')
source('R/helper_fxns.R')
theme_set(theme_minimal(base_size = 17))


# method 1 ----
# status quo goa ratio

m1 <- sable %>% 
  tidyr::expand(species, year, strata) %>% 
  left_join(sable) %>% 
  left_join(
    sable %>% 
      filter(!(strata %in% c("Aleutians", "Bering Sea"))) %>% 
      group_by(year) %>% 
      # GOA wide CPUE/RPN/RPWs
      summarize(goa_rpn = sum(rpn, na.rm = TRUE),
                goa_rpw = sum(rpw, na.rm = TRUE)) %>% 
      ungroup() %>%
      # Lag GOA CPUE/RPN/RPWs, i.e. GOA_{j} and GOA_{j-1}
      mutate(goa_rpn_lyr = lag(goa_rpn, default = goa_rpn[1]),
             goa_rpw_lyr = lag(goa_rpw, default = goa_rpw[1]))
  ) %>%
  # Lag RPN/RPW values by Council area (only AI and BS end up getting
  # used); also get modern mean CV for interpolation
  group_by(strata) %>% 
  mutate(rpn_lyr = lag(rpn, default = rpn[1]),
         rpw_lyr = lag(rpw, default = rpw[1]),
         mean_cv = ifelse(year >= 1996, mean(sqrt(rpn_var) / rpn, na.rm = TRUE), 0.05),
         mean_cv = ifelse(is.nan(mean_cv), 0, mean_cv)) %>% 
  ungroup() %>% 
  # Extrapolate missing BSAI years using GOA_{j} / GOA_{j-1}
  mutate(rpn = ifelse(strata %in% c("Aleutians", "Bering Sea") & is.na(rpn),
                      rpn_lyr * (goa_rpn / goa_rpn_lyr), rpn),
         rpw = ifelse(strata %in% c("Aleutians", "Bering Sea") & is.na(rpw),
                      rpw_lyr * (goa_rpw / goa_rpw_lyr), rpw)) %>% 
  # Assumed mean CV for these interpolations
  mutate(rpn_var = ifelse(strata %in% c("Aleutians", "Bering Sea") & 
                            between(year, 1996, max(sable$year)) & 
                            !is.na(rpn) & is.na(rpn_var), 
                          (mean_cv * rpn)^2, rpn_var),
         rpw_var = ifelse(strata %in% c("Aleutians", "Bering Sea") & 
                            between(year, 1996, max(sable$year)) & 
                            !is.na(rpw) & is.na(rpw_var), 
                          (mean_cv * rpw)^2, rpw_var)) %>% 
  select(-c(contains(c("goa", "lyr")), mean_cv)) %>% 
  filter(year >= 1996) %>% 
  mutate(method = "statquo_goaratio") %>% 
  arrange(strata, year)

# method 2 ----

# wgoa ratio

m2 <- sable %>% 
  tidyr::expand(species, year, strata) %>% 
  left_join(sable) %>% 
  left_join(
    sable %>% 
      filter(strata == 'Western Gulf of Alaska') %>% 
      # Lag WGOA CPUE/RPN/RPWs, i.e. GOA_{j} and GOA_{j-1}
      mutate(wgoa_rpn = rpn,
             wgoa_rpw = rpw,
             wgoa_rpn_lyr = lag(rpn, default = rpn[1]),
             wgoa_rpw_lyr = lag(rpw, default = rpw[1])) %>% 
      select(species, year, contains(c("goa", "lyr")))
  ) %>%
  # Lag RPN/RPW values by Council area (only AI and BS end up getting
  # used); also get modern mean CV for interpolation
  group_by(strata) %>% 
  mutate(rpn_lyr = lag(rpn, default = rpn[1]),
         rpw_lyr = lag(rpw, default = rpw[1]),
         mean_cv = ifelse(year >= 1996, mean(sqrt(rpn_var) / rpn, na.rm = TRUE), 0.05),
         mean_cv = ifelse(is.nan(mean_cv), 0, mean_cv)) %>% 
  ungroup() %>% 
  # Extrapolate missing BSAI years using GOA_{j} / GOA_{j-1}
  mutate(rpn = ifelse(strata %in% c("Aleutians", "Bering Sea") & is.na(rpn),
                      rpn_lyr * (wgoa_rpn / wgoa_rpn_lyr), rpn),
         rpw = ifelse(strata %in% c("Aleutians", "Bering Sea") & is.na(rpw),
                      rpw_lyr * (wgoa_rpw / wgoa_rpw_lyr), rpw)) %>% 
  # Assumed mean CV for these interpolations
  mutate(rpn_var = ifelse(strata %in% c("Aleutians", "Bering Sea") & 
                            between(year, 1996, max(sable$year)) & 
                            !is.na(rpn) & is.na(rpn_var), 
                          (mean_cv * rpn)^2, rpn_var),
         rpw_var = ifelse(strata %in% c("Aleutians", "Bering Sea") & 
                            between(year, 1996, max(sable$year)) & 
                            !is.na(rpw) & is.na(rpw_var), 
                          (mean_cv * rpw)^2, rpw_var)) %>% 
  select(-c(contains(c("wgoa", "lyr")), mean_cv)) %>% 
  filter(year >= 1996) %>% 
  mutate(method = "wgoaratio") %>% 
  arrange(strata, year)

# method 3 ----
# linear interpolation
m3 <- f_linapprox(df = sable, method_lab = "linear_approx") %>% 
  select(-contains('cv'))

# method 4 -----
library(TMB)
compile("src/ar1.cpp")
dyn.load(dynlib("src/ar1"))

arsable <- sable %>% 
  filter(strata %in% c('Aleutians','Bering Sea')) %>% 
  mutate(rpn_cv = sqrt(rpn_var) / rpn,
         rpw_cv = sqrt(rpw_var) / rpw)

model_yrs = sort(unique(sable$year))

arsable <- arsable %>% 
  tidyr::expand(strata, year = model_yrs) %>% 
  left_join(arsable) 


# rpns first
data <- list(model_yrs = model_yrs,
             obs = arsable %>% 
               pivot_wider(id_cols = year, names_from = strata, values_from = rpn) %>% 
               select(-year) %>% 
               as.matrix(),
             obs_cv = arsable %>% 
               pivot_wider(id_cols = year, names_from = strata, values_from = rpn_cv) %>% 
               select(-year) %>% 
               as.matrix())

par <- list(log_PE = rep(log(1), ncol(data$obs)), 
            logit_rho = rep(genlogit(0.5, -1, 1), ncol(data$obs)),
            log_pred = log(apply(X = data$obs, MARGIN = 2, FUN = zoo::na.approx, maxgap = 100, rule = 2)))

map <- par
map$log_PE <- as.factor(1:length(map$log_PE))
map$log_PE <- as.factor(rep(1, length(map$log_PE)))
map$logit_rho <- as.factor(1:length(map$logit_rho))
map$logit_rho <- as.factor(rep(1, length(map$logit_rho)))
map$log_pred <- as.factor(1:length(map$log_pred))

mod <- MakeADFun(data, par, map, random = "log_pred", DLL = "ar1")
# input$data = mod$simulate(complete=TRUE)
opt = nlminb(mod$par, mod$fn, mod$gr)

mod$rep = mod$report()
mod$rep$pred
mod$sdrep = sdreport(mod)
summary(mod$sdrep)
mod$env$parList()$log_pred

pars = as.list(mod$sdrep, "Est")
sd = as.list(mod$sdrep, "Std")

ci(pars$log_PE[1], sd$log_PE[1], type = "exp") 
ci(pars$logit_rho[1], sd$logit_rho[1], lo = -1, hi = 1, type = "expit", k = 2) 
# rho is ~1, likely this model is overparameterized

# AI
pred <- matrix(data = NA, nrow = nrow(pars$log_pred), ncol = 4)
for(i in 1:nrow(pars$log_pred)) pred[i,] <- ci(pars$log_pred[i,1], sd$log_pred[i,1], type = "exp")
pred <- as.data.frame(pred); names(pred) <- c('pred', 'se', 'lci', 'uci')
preddf <- pred %>% mutate(year = model_yrs, strata = 'Aleutians')
# preddf <- pred %>% mutate(year = model_yrs, strata = 'Bering Sea')
# Bering
pred <- matrix(data = NA, nrow = nrow(pars$log_pred), ncol = 4)
for(i in 1:nrow(pars$log_pred)) pred[i,] <- ci(pars$log_pred[i,2], sd$log_pred[i,2], type = "exp")
pred <- as.data.frame(pred); names(pred) <- c('pred', 'se', 'lci', 'uci')
preddf <- preddf %>% bind_rows(pred %>% mutate(year = model_yrs, strata = 'Bering Sea'))

outrpn <- arsable %>% 
  # filter(strata %in% c('Aleutians', 'Bering Sea')) %>% 
  mutate(log_rpn = ifelse(rpn > 0, log(rpn), NA),
         sd_log_rpn = ifelse(rpn > 0, sqrt(log(rpn_cv^2 + 1)), NA),
         rpn_lci = exp(log_rpn - qnorm(0.975) * sd_log_rpn),
         rpn_uci = exp(log_rpn + qnorm(0.975) * sd_log_rpn)) %>% 
  left_join(preddf %>% 
              rename(rpn_pred = pred, rpn_se = se, rpn_pred_lci = lci, rpn_pred_uci = uci)) 

outrpn %>% 
  ggplot(aes(x = year, y = rpn)) +
  geom_point() + 
  geom_errorbar(aes(ymin = rpn_lci, ymax = rpn_uci)) +
  geom_line(aes(y = rpn_pred), col = 'red') +
  geom_ribbon(aes(ymin = rpn_pred_lci, ymax = rpn_pred_uci), col = NA, fill = 'red', alpha = 0.2) +
  facet_wrap(~strata)

# rpws second
data <- list(model_yrs = model_yrs,
             obs = arsable %>% 
               pivot_wider(id_cols = year, names_from = strata, values_from = rpw) %>% 
               select(-year) %>% 
               as.matrix(),
             obs_cv = arsable %>% 
               pivot_wider(id_cols = year, names_from = strata, values_from = rpw_cv) %>% 
               select(-year) %>% 
               as.matrix())

par <- list(log_PE = rep(log(1), ncol(data$obs)), 
            logit_rho = rep(genlogit(0.5, -1, 1), ncol(data$obs)),
            log_pred = log(apply(X = data$obs, MARGIN = 2, FUN = zoo::na.approx, maxgap = 100, rule = 2)))

map <- par
map$log_PE <- as.factor(1:length(map$log_PE))
map$log_PE <- as.factor(rep(1, length(map$log_PE)))
map$logit_rho <- as.factor(1:length(map$logit_rho))
map$logit_rho <- as.factor(rep(1, length(map$logit_rho)))
map$log_pred <- as.factor(1:length(map$log_pred))

mod <- MakeADFun(data, par, map, random = "log_pred", DLL = "ar1")
# input$data = mod$simulate(complete=TRUE)
opt = nlminb(mod$par, mod$fn, mod$gr)

mod$rep = mod$report()
mod$rep$pred
mod$sdrep = sdreport(mod)
summary(mod$sdrep)
mod$env$parList()$log_pred

pars = as.list(mod$sdrep, "Est")
sd = as.list(mod$sdrep, "Std")

ci(pars$log_PE[1], sd$log_PE[1], type = "exp") 
ci(pars$logit_rho[1], sd$logit_rho[1], lo = -1, hi = 1, type = "expit", k = 2) 
# rho is ~1, likely this model is overparameterized

# AI
pred <- matrix(data = NA, nrow = nrow(pars$log_pred), ncol = 4)
for(i in 1:nrow(pars$log_pred)) pred[i,] <- ci(pars$log_pred[i,1], sd$log_pred[i,1], type = "exp")
pred <- as.data.frame(pred); names(pred) <- c('pred', 'se', 'lci', 'uci')
preddf <- pred %>% mutate(year = model_yrs, strata = 'Aleutians')
# preddf <- pred %>% mutate(year = model_yrs, strata = 'Bering Sea')
# Bering
pred <- matrix(data = NA, nrow = nrow(pars$log_pred), ncol = 4)
for(i in 1:nrow(pars$log_pred)) pred[i,] <- ci(pars$log_pred[i,2], sd$log_pred[i,2], type = "exp")
pred <- as.data.frame(pred); names(pred) <- c('pred', 'se', 'lci', 'uci')
preddf <- preddf %>% bind_rows(pred %>% mutate(year = model_yrs, strata = 'Bering Sea'))

outrpw <- arsable %>% 
  # filter(strata %in% c('Aleutians', 'Bering Sea')) %>% 
  mutate(log_rpw = ifelse(rpw > 0, log(rpw), NA),
         sd_log_rpw = ifelse(rpw > 0, sqrt(log(rpw_cv^2 + 1)), NA),
         rpw_lci = exp(log_rpw - qnorm(0.975) * sd_log_rpw),
         rpw_uci = exp(log_rpw + qnorm(0.975) * sd_log_rpw)) %>% 
  left_join(preddf %>% 
              rename(rpw_pred = pred, rpw_se = se, rpw_pred_lci = lci, rpw_pred_uci = uci)) 

outrpw %>% 
  ggplot(aes(x = year, y = rpw)) +
  geom_point() + 
  geom_errorbar(aes(ymin = rpw_lci, ymax = rpw_uci)) +
  geom_line(aes(y = rpw_pred), col = 'red') +
  geom_ribbon(aes(ymin = rpw_pred_lci, ymax = rpw_pred_uci), col = NA, fill = 'red', alpha = 0.2) +
  facet_wrap(~strata)

m4 <- outrpn %>% 
  mutate(species = 'Sablefish',
         rpn = ifelse(is.na(rpn), rpn_pred, rpn)) %>% 
  select(species, year, strata, rpn, rpn_var) %>% 
  left_join(outrpw %>% mutate(rpw = ifelse(is.na(rpw), rpw_pred, rpw)) %>% 
              select(year, strata, rpw, rpw_var)) %>% 
  bind_rows(sable %>% filter(!strata %in% c('Aleutians','Bering Sea'))) %>% 
  mutate(method = 'ar1') %>% 
  filter(year >= 1996)

# compare ----

allm <- m1 %>% bind_rows(m2, m3, m4) 

allm <- allm %>% 
  left_join(allm %>% 
              filter((strata == "Aleutians" & year %% 2 == 0) |
                       (strata == "Bering Sea" & year %% 2 != 0) |
                       grepl("Gulf of Alaska", strata)) %>% 
              select(species, year, strata, obs_rpn = rpn, obs_rpw = rpw, method)) %>% 
  mutate(strata = factor(strata, ordered = TRUE, levels = c('Aleutians', 'Bering Sea', 
                                                            'Western Gulf of Alaska', 'Central Gulf of Alaska', 
                                                            'Eastern Gulf of Alaska')),
         method = factor(method, ordered = TRUE, levels = c('statquo_goaratio', 'wgoaratio', 'linear_approx', 'ar1')))

perf <- allm %>% 
  filter(strata %in% c("Aleutians", "Bering Sea")) %>% 
  group_by(species, method, strata) %>% 
  mutate(next_obs_rpn = lead(obs_rpn, 1),
         next_obs_rpw = lead(obs_rpw, 1),
         rpn_pdiff = 100 * abs(next_obs_rpn - rpn)/((next_obs_rpn + rpn)/2),
         rpw_pdiff =  100 * abs(next_obs_rpw - rpw)/((next_obs_rpw + rpw)/2)) %>% 
  group_by(species, strata, method) %>% 
  summarise(bsai_rpn_stat = mean(rpn_pdiff, na.rm = TRUE),
         bsai_rpw_stat = mean(rpw_pdiff, na.rm = TRUE)) %>%
  arrange(species, strata, bsai_rpn_stat)
perf

allm %>% 
  ggplot(aes(x = year, y = obs_rpn)) +
  geom_line(aes(y = rpn, col = method, lty = method), size = 1) +
  geom_point(aes(y = rpn, col = method, shape = method), size = 2) +
  geom_point() +
  facet_wrap(~strata, ncol = 1) +
  ggthemes::scale_color_colorblind()

allm %>% 
  ggplot(aes(x = year, y = obs_rpw)) +
  geom_line(aes(y = rpw, col = method, lty = method), size = 1) +
  geom_point(aes(y = rpw, col = method, shape = method), size = 2) +
  geom_point() +
  facet_wrap(~strata, ncol = 1) +
  ggthemes::scale_color_colorblind()

# akwide ----

akwide <- allm %>% 
  group_by(species, method, year) %>% 
  summarise(rpn = sum(rpn),
            rpw = sum(rpw))# %>% 
  # bind_rows(read_csv('data/sableakwide.csv') %>% 
  #             select(species, year, rpn, rpw) %>% 
  #             mutate(method = 'db_akwide')) 

akwide %>% 
  ggplot(aes(x = year, y = rpn)) +
  geom_line(aes(y = rpn, col = method, lty = method), size = 1) +
  geom_point(aes(y = rpn, col = method, shape = method), size = 2) +
  ggthemes::scale_color_colorblind()

akwide %>% 
  ggplot(aes(x = year, y = rpw)) +
  geom_line(aes(y = rpw, col = method, lty = method), size = 1) +
  geom_point(aes(y = rpw, col = method, shape = method), size = 2) +
  ggthemes::scale_color_colorblind()

  

