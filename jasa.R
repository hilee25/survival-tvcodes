### 데이터 적용 ###
data <- jasa

cp <- make_cp_from_dates(
  data = jasa,
  accept_col = "accept.dt",
  tx_col     = "tx.date",
  fu_col     = "fu.date",
  status_col = "fustat",
  eps = 1e-5,
  fix_zero_followup = "bump"
)

head(cp)

res <- plot_all_landmarks(
  cp_df = cp,
  landmarks = c(30, 60, 90, 120),
  show_risktable = TRUE,         
  risk_table_height = 0.35,      
  risk_table_base_size = 9,      
  theme_fn = ggplot2::theme_bw,
  legend_position = "top",
  ncol = 2,
  conf_int = FALSE
)
print(res$plot)



### 동일성 확인 ###
## Cox(time-dependent) & Score test(β=0)
fit <- coxph(Surv(tstart, tstop, event) ~ response,
             data = cp, ties = "efron")
s <- summary(fit)                         # s$sctest: Score (logrank) test
# summary.coxph의 sctest는 보통 c(Chisq, df, p) 순서
chisq_score <- as.numeric(s$sctest[1]); chisq_score
df_score    <- as.integer(s$sctest[2]); df_score
p_score     <- as.numeric(s$sctest[3]); p_score

## Mantel–Byar (log-rank on counting-process)
mb_res <- mantel_byar_cp(cp)

# 값 비교하여 표시
cat(sprintf("Mantel–Byar (custom): chisq=%.4f, p=%.6f\n", mb_res$chisq, mb_res$p))
cat(sprintf("Cox Score test: chisq=%.4f, p=%.6f\n", chisq_score, p_score))
