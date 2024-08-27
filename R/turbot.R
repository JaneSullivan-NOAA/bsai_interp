
# compare BS and AI ----

compare_rpns <- f_ltmean(rpn_old, method_lab = 'statusquo_meanratio_olddata', 
                         mean_syr = 1996, index_syr = 1996) %>% 
  bind_rows(f_ltmean(rpn, method_lab = 'statusquo_meanratio_newdat', 
                     mean_syr = 1996, index_syr = 1996)) %>% 
  bind_rows(f_linapprox(rpn, method_lab = 'new_linearapprox_newdat', index_syr = 1996))

ggplot(compare_rpns, aes(x = year, y = rpn / 1e3, 
                         col = method, lty = method, shape = method)) +
  geom_point() +
  geom_line() +
  facet_wrap(~fmp) +
  labs(x = "Year", y = "RPN",
       title ="Greenland turbot Relative Population Numbers",
       subtitle = "Comparing methods for interpolating missing survey years (odd = BS, even = AI)")

ggsave("results/turbot_rpnmethods_fmp.png", units = "in", 
       width = 8, height = 5)  

# compare full BSAI index ----

index <- compare_rpns %>% 
  group_by(year, method) %>% 
  dplyr::summarise(rpn = sum(rpn, na.rm = TRUE),
                   rpn_var = sum(rpn_var, na.rm = TRUE)) %>% 
  # assume fixed cv = 0.2 for old 'longterm mean' method
  mutate(rpn_var = ifelse(is.na(rpn_var) | rpn_var == 0, (rpn * 0.2)^2, rpn_var),
         lci = rpn - 1.96 * sqrt(rpn_var),
         uci = rpn + 1.96 * sqrt(rpn_var)) #%>% View()

ggplot(index, aes(x = year, y = rpn / 1e3,
                  ymin = lci / 1e3, ymax = uci / 1e3, 
                  fill = method, group = method,
                  col = method, lty = method, shape = method)) +
  geom_ribbon(col = 'white', alpha = 0.2) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "RPN",
       title ="Greenland turbot Relative Population Numbers with 95% CI",
       subtitle = "Comparing methods for interpolating missing survey years (odd = BS, even = AI)")

ggsave("results/turbot_rpnmethods.png", units = "in", 
       width = 8, height = 5)  
