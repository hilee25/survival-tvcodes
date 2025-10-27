# ==============================================
# Step 2) fit_for_L: 단일 랜드마크 L에서 추정/검정
# ==============================================

#' fit_for_L
#' @param cp_df data.frame. columns(id, tstart, tstop, event, response)
#' @param L numeric. landmark 시간(원 단위와 동일)
#' @param robust logical. cluster(id) 기반 분산 사용 (기본 TRUE)
#'
#' @return list(L, n_total, n_events, n_group, hr, ci_low, ci_high, p, fit)

fit_for_L <- function(cp_df, L, robust = TRUE) {
  stopifnot(all(c("id","tstart","tstop","event","response") %in% names(cp_df)))
  
  # 1) L 시점 at-risk 정의
  had_event_before_L <- cp_df |>
    group_by(id) |>
    summarize(ev_before = any(event == 1 & tstop <= L), .groups = "drop")
  
  at_risk_ids <- cp_df |>
    group_by(id) |>
    summarize(has_followup = any(tstop > L), .groups = "drop")
  
  keep_ids <- had_event_before_L |>
    inner_join(at_risk_ids, by = "id") |>
    filter(!ev_before & has_followup) |>
    pull(id)
  
  df_risk <- cp_df |> filter(id %in% keep_ids)
  
  # 2) L을 덮는 인터벌의 response 값으로 그룹 고정
  zL <- df_risk |>
    filter(tstart <= L, tstop > L) |>
    select(id, response) |>
    distinct() |>
    rename(Z_L = response)
  
  # 3) L 이후 구간만 남기기
  postL <- df_risk |>
    filter(tstop > L) |>
    mutate(tstart = pmax(tstart, L)) |>
    filter(tstop > tstart) |>
    left_join(zL, by = "id") |>
    mutate(Z_L = factor(Z_L, levels = c(0,1), labels = c("Z=0","Z=1")))
  
  # 4) Cox 적합 (Mantel–Byar = Score test)
  f <- as.formula(Surv(tstart, tstop, event) ~ Z_L)
  if (robust) {
    fit <- coxph(f, data = postL, robust = TRUE, cluster = id)
  } else {
    fit <- coxph(f, data = postL)
  }
  
  s  <- summary(fit)
  hr <- unname(exp(coef(fit))[1])
  ci <- unname(exp(confint(fit))[1, ])
  p  <- as.numeric(s$sctest[3])
  
  list(
    L        = L,
    n_total  = dplyr::n_distinct(postL$id),
    n_events = sum(postL$event),
    n_group  = table(postL$Z_L) |> unclass(),
    hr       = as.numeric(hr),
    ci_low   = as.numeric(ci[1]),
    ci_high  = as.numeric(ci[2]),
    p        = p,
    fit      = fit
  )
}
