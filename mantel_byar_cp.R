# ==============================================
# mantel_byar_cp(): Mantel-Byar 검정 함수
# 각 사건 시점에서 위험집단을 재구성하여 Mantel-Byar 검정을 수행
# - 입력: counting-process 데이터(id, tstart, tstop, event, response), ref_val(기준 상태 값)
# - 처리: 사건 발생 시점들의 위험집단에서 계산
# - 출력: list(chisq, p, D_ref, E, V, contrib(시점별 표), ref_val)
# ==============================================
#' @param cp data.frame. columns(id, tstart, tstop, event, response).
#' @param ref_val 0 또는 1. 통상 0을 참조군으로 둠.
#'
#' @return list(chisq, p, D_ref, E, V, contrib, ref_val)

mantel_byar_cp <- function(cp, ref_val = 0L) {
  # --- (0) 입력 점검 ---
  req <- c("tstart","tstop","event", "response")
  stopifnot(all(req %in% names(cp)))
  if (!is.numeric(cp$tstart) || !is.numeric(cp$tstop)) {
    stop("start/stop must be numeric times (e.g., days since baseline).")
  }
  if (any(cp$tstop <= cp$tstart, na.rm = TRUE)) {
    stop("Found non-positive intervals: require start < stop for all rows.")
  }
  
  # --- (1) response를 이진(0/1)으로 강제 ---
  z <- cp[["response"]]
  if (is.factor(z)) z <- as.character(z)
  z <- as.integer(z)
  badz <- which(!is.na(z) & !(z %in% c(0L,1L)))
  if (length(badz)) stop("response must be binary 0/1; offending rows: ", paste(head(badz,5), collapse=", "))
  cp[["response"]] <- z
  
  # --- (2) 사건 시점 추출 ---
  # event==1인 행들의 tstop을 오름차순 고유값으로 정리
  evt_times <- sort(unique(cp$tstop[cp$event == 1]))
  if (length(evt_times) == 0) {
    # 사건이 하나도 없으면 통계량=0, p=1 (정보 없음)
    return(list(chisq = 0, p = 1, D_ref = 0, E = 0, V = 0,
                contrib = data.frame(time = numeric(0), d = integer(0),
                                     N = integer(0), Nref = integer(0),
                                     dref = integer(0), OE = numeric(0), V = numeric(0))))
  }
  
  # --- (3) 누적량 초기화 ---
  D_ref <- 0  # Σ d_ref : 기준(ref) 상태에서 관측된 사건 수 합
  E_sum <- 0  # Σ E_ref : 기준(ref) 상태에서의 기대 사건 수 합
  V_sum <- 0  # Σ V_ref : 분산 합
  
  contrib <- vector("list", length(evt_times))
  
  # --- (4) 사건 시점별 합산 ---
  for (i in seq_along(evt_times)) {
    tt <- evt_times[i]
    
    # 위험집단: tstart < tt <= tstop
    at_risk <- cp[cp$tstart < tt & cp$tstop >= tt, , drop = FALSE]
    if (nrow(at_risk) == 0) {
      contrib[[i]] <- data.frame(time = tt, d = 0L, N = 0L, Nref = 0L,
                                 dref = 0L, OE = 0, V = 0)
      next
    }
    
    N    <- nrow(at_risk)
    Nref <- sum(at_risk[["response"]] == ref_val, na.rm = TRUE)
    
    # 해당 시점의 사건들
    at_evt <- cp[cp$tstop == tt & cp$event == 1L, , drop = FALSE]
    d      <- nrow(at_evt)
    if (d == 0L) {
      contrib[[i]] <- data.frame(time = tt, d = 0L, N = N, Nref = Nref,
                                 dref = 0L, OE = 0, V = 0)
      next
    }
    dref   <- sum(at_evt[["response"]] == ref_val, na.rm = TRUE)
    
    # 한쪽 상태만 존재하거나 N<=1이면 정보 없음
    if (N <= 1L || Nref == 0L || Nref == N) {
      contrib[[i]] <- data.frame(time = tt, d = d, N = N, Nref = Nref,
                                 dref = dref, OE = 0, V = 0)
      next
    }
    
    # 기대/분산 (동시사건 보정)
    E  <- d * (Nref / N)
    V  <- d * (Nref * (N - Nref) / (N^2)) * ((N - d) / (N - 1))
    
    OE <- (dref - E)
    
    # 누적
    D_ref <- D_ref + dref
    E_sum <- E_sum + E
    V_sum <- V_sum + V
    
    contrib[[i]] <- data.frame(time = tt, d = d, N = N, Nref = Nref,
                               dref = dref, OE = OE, V = V)
  }
  
  contrib <- do.call(rbind, contrib)
  
  # --- (5) 최종 검정통계량 ---
  if (!is.finite(V_sum) || V_sum <= 0) {
    chisq <- 0
    pval  <- 1
  } else {
    chisq <- ( (D_ref - E_sum)^2 ) / V_sum
    pval  <- pchisq(chisq, df = 1, lower.tail = FALSE)
  }
  
  # ---- (6) 반환 ---
  list(chisq = chisq, p = pval, D_ref = D_ref, E = E_sum, V = V_sum,
       contrib = contrib, ref_val = ref_val)
}
