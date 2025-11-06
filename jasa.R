# =========================================================
# jasa.R: 데이터 분석 전체 재현 코드
# =========================================================

## 필요한 패키지
library(survival); library(dplyr); library(ggplot2); library(purrr); library(tibble); library(patchwork); library(survminer)
options(scipen = 999) # 지수표기 억제


## (3.2) 데이터 적용
cp <- make_cp_from_dates(
  data = jasa,
  accept_col = "accept.dt",
  tx_col     = "tx.date",
  fu_col     = "fu.date",
  status_col = "fustat",
  eps = 1e-5,
  fix_zero_followup = "bump"
)


## (3.3.1) jasa 데이터 적용: L = 60일
result_60 <- fit_for_L(cp, L = 60, robust = TRUE)

cat(sprintf("60일 랜드마크: HR=%.2f (95%% CI: %.2f-%.2f), p=%.3f\n",
            result_60$hr, result_60$ci_low, result_60$ci_high, result_60$p))


## (3.3.2) jasa 데이터 적용: L = 30, 60, 90, 120일
# 다중 랜드마크 분석 및 시각화
res <- plot_all_landmarks(
  cp_df = cp,
  landmarks = c(30, 60, 90, 120),
  robust = TRUE,
  conf_int = FALSE,
  show_risktable = TRUE,
  risk_table_height = 0.35,
  ncol = 2,
  annotate_cox_p = FALSE
)

# 패널 그래프 출력
print(res$plot)

# 요약 테이블
print(res$table)


## (3.4.1)
mb_result <- mantel_byar_cp(cp, ref_val = 0)

cat(sprintf("Mantel-Byar 검정:\n"))
cat(sprintf("  관찰 사망자 (비이식): %d\n", mb_result$D_ref))
cat(sprintf("  기대 사망자: %.2f\n", mb_result$E))
cat(sprintf("  χ² = %.4f\n", mb_result$chisq))
cat(sprintf("  p-value = %.4f\n", mb_result$p))


## (3.4.2) 
# 시간종속 Cox 모형
fit_tdcox <- coxph(Surv(tstart, tstop, event) ~ response,
                   data = cp, ties = "efron")

# 요약
summary(fit_tdcox)

# Score 검정
s <- summary(fit_tdcox)
cat(sprintf("\nCox Score 검정:\n"))
cat(sprintf("  χ² = %.4f\n", s$sctest[1]))
cat(sprintf("  p-value = %.4f\n", s$sctest[3]))


## (3.4.3)
# Schoenfeld 잔차 검정
test_ph <- cox.zph(fit_tdcox)
print(test_ph)

# 시각화
ggcoxzph(test_ph)
