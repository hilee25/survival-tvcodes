# ==============================================
# fit_for_L(): 단일 랜드마크 분석 함수
# 단일 랜드마크 시점 L에서 그룹 고정 후 Cox 추정/검정
# - 입력: counting-process long 데이터(id, tstart, tstop, event, response)
# - 처리: L 직후 at-risk 집단만 남기고, L을 덮는 구간의 response 값을 Z_L로 고정
# - 출력: HR, 95% CI, score test p-value 등 요약 + coxph 적합 객체
# ==============================================
#' @param cp_df data.frame. columns(id, tstart, tstop, event, response).
#' @param L numeric. 랜드마크 시점.
#' @param robust logical. cluster(id) 기반 강건 분산 사용 여부 (기본 TRUE).
#'
#' @return list(L, n_total, n_events, n_group, hr, ci_low, ci_high, p, fit)

fit_for_L <- function(cp_df, L, robust = TRUE) {
  # --- (0) 입력 유효성 검사 ---
  stopifnot(all(c("id","tstart","tstop","event","response") %in% names(cp_df)))
  
  # --- (1) L 이전 시점 발생자 제외 ---
  # L 이전에 이미 사건이 있었는지 확인
  had_event_before_L <- cp_df |>
    group_by(id) |>
    summarize(ev_before = any(event == 1 & tstop <= L), .groups = "drop")

  # L 이후로 관찰이 이어지는지 확인
  at_risk_ids <- cp_df |>
    group_by(id) |>
    summarize(has_followup = any(tstop > L), .groups = "drop")
  
  # --- (2) L 이후 추적 존재 확인 ---  
  # L 이전 사건 없음 & L 이후 추적 있음 → 랜드마크 at-risk 집단
  keep_ids <- had_event_before_L |>
    inner_join(at_risk_ids, by = "id") |>
    filter(!ev_before & has_followup) |>
    pull(id)
  
  df_risk <- cp_df |> filter(id %in% keep_ids)

  # --- (3) L 시점의 공변량 상태 고정 (Z_L) ---  
  # tstart <= L < tstop 인 구간의 response를 Z_L로 가져옴 (개체별 1개 값)
  zL <- df_risk |>
    filter(tstart <= L, tstop > L) |>
    select(id, response) |>
    distinct() |>
    rename(Z_L = response)
  
  # --- (4) 추적 구간을 [L, tstop)로 좌절단 ---  
  # 시작점을 pmax(tstart, L)로 잘라서 L 직후부터의 위험구간만 유지
  postL <- df_risk |> 
    filter(tstop > L) |>
    mutate(tstart = pmax(tstart, L)) |>
    filter(tstop > tstart) |>
    left_join(zL, by = "id") |>
    mutate(Z_L = factor(Z_L, levels = c(0,1), labels = c("Z=0","Z=1")))
  
  # --- (5) Cox 모형 적합 Surv(tstart, tstop, event) ~ Z_L ---
  f <- as.formula(Surv(tstart, tstop, event) ~ Z_L)
  if (robust) {
    fit <- coxph(f, data = postL, robust = TRUE, cluster = id)
  } else {
    fit <- coxph(f, data = postL)
  }
  
  # --- (6) HR, 신뢰구간, p값 추출 및 반환 ---
  s  <- summary(fit)
  hr <- unname(exp(coef(fit))[1])  # 위험비(HR)
  ci <- unname(exp(confint(fit))[1, ])  # 95% CI
  p  <- as.numeric(s$sctest[3])  # Score test p-value
  
  list(
    L        = L,
    n_total  = n_distinct(postL$id),  # L 이후 분석에 실제로 들어간 개체 수
    n_events = sum(postL$event),  # L 이후 관찰된 사건 수
    n_group  = table(distinct(postL, id, Z_L)$Z_L) |> unclass(),  # L 시점 그룹 분포
    hr       = as.numeric(hr),
    ci_low   = as.numeric(ci[1]),
    ci_high  = as.numeric(ci[2]),
    p        = p,
    fit      = fit
  )
}
