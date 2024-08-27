# different methods for calculating off-year surveys

# use long-term mean to get average proportion by area to approximate BS or AI in missing survey years
f_ltmean <- function(df, # cols = fmp ('Aleutians', 'Bering Sea'), year, rpn
                     method_lab = "ltmean_olddat", # label for method
                     mean_syr = 1996, # starting year for long-term mean calcs
                     index_syr = 1996 # starting year for index
) {
  
  df <- df %>% 
    right_join(expand_grid(fmp = unique(df$fmp),
                           year = min(df$year):max(df$year)))
  
  df <- df %>% 
    pivot_wider(id_cols = year, values_from = rpn, names_from = fmp) %>% 
    left_join(df %>% 
                group_by(fmp) %>% 
                mutate(mean_rpn = mean(rpn[year >= mean_syr], na.rm = TRUE)) %>% 
                ungroup() %>% 
                mutate(prop = mean_rpn / sum(unique(mean_rpn), na.rm = TRUE)) %>% 
                filter(!is.na(rpn)) %>% 
                select(year, prop)) %>% 
    mutate(bsairpn = ifelse(!is.na(Aleutians), Aleutians / prop, `Bering Sea` / prop),
           Aleutians = ifelse(is.na(Aleutians), bsairpn - `Bering Sea`, Aleutians),
           `Bering Sea` = ifelse(is.na(`Bering Sea`), bsairpn - Aleutians, `Bering Sea`)) %>% 
    pivot_longer(cols = c(Aleutians, `Bering Sea`), names_to = "fmp", values_to = "rpn") %>% 
    select(year, fmp, rpn) %>% 
    mutate(method = method_lab) %>% 
    arrange(fmp, year) %>% 
    filter(year >= index_syr)
  
  return(df)
}

# use linear interpolation to get rpns/cvs in off-survey years; end years set
# equal to nearest year
f_linapprox <- function(df, # cols = fmp ('Aleutians', 'Bering Sea'), year, rpn, rpn_var
                        method_lab = "lininterp_newdat", # label for method
                        index_syr = 1996 # starting year for index
) {
  
  df <- df %>% 
    mutate(rpn_cv = sqrt(rpn_var) / rpn,
           rpw_cv = sqrt(rpw_var) / rpw) %>% 
    right_join(df %>% tidyr::expand(species, strata, year)) %>% 
    # expand_grid(strata = unique(df$strata),
    #                        year = unique(df$year))) %>% 
    filter(year >= index_syr) %>% 
    arrange(strata, year) %>% 
    group_by(strata) %>% 
    mutate(rpn = zoo::na.approx(rpn, maxgap = 20, rule = 2),
           rpn_cv = zoo::na.approx(rpn_cv, maxgap = 20, rule = 2),
           rpn_var = ifelse(is.na(rpn_var), (rpn_cv * rpn)^2, rpn_var),
           rpw = zoo::na.approx(rpw, maxgap = 20, rule = 2),
           rpw_cv = zoo::na.approx(rpw_cv, maxgap = 20, rule = 2),
           rpw_var = ifelse(is.na(rpw_var), (rpw_cv * rpw)^2, rpw_var)
           ) %>%
    mutate(method = method_lab) %>% 
    ungroup()
  
  return(df)
}

genlogit <- function(x, low=0, upp=1){
  p = (x-low)/(upp-low)
  return(log(p/(1-p)))
}

geninvlogit <- function(x, low=0, upp=1){
  pprime = exp(x)/(1+exp(x))
  return(pprime*(upp-low)+low)
}

ci = function(par, se, p=0.975, lo = 0, hi = 1, type = "I", k = 1){
  ci = par + c(-1,1) * qnorm(0.975) * se
  if(type == "I") {
    return(c(par, se, ci))
  }
  if(type == "exp") {
    return(c(exp(par), exp(par)*se, exp(ci)))
  }
  if(type == "expit") { #Delta-method: V(lo + (hi-lo)/(1 + exp(-x))) ~ ((hi-lo) * p * (1-p))^2 * V(x)
    p = 1/(1 + exp(- k * par))
    dm.se = k * abs(hi-lo)*p*(1-p)*se
    return(c(geninvlogit(par, lo, hi), dm.se, lo + (hi-lo)/(1 + exp(-ci))))
  }
}

