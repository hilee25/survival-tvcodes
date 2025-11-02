# =========================================================
# make_cp_from_dates(): Counting process 자료 변환 함수
# 날짜 기반 원자료 → counting-process + time-dependent binary covariate 생성
# - 입력: landmark_days, accept, response, follow-up end, status
# - 출력: id, tstart, tstop, event, response 가 있는 data.frame
# - 용도: coxph(Surv(tstart, tstop, event) ~ response) 분석에 바로 사용
# =========================================================
#' @param data            분석할 data.frame
#' @param landmark_days   numeric 벡터. baseline(accept)로부터 며칠 후를 landmark로 지정 (예: c(60, 120)).
#' @param accept_col      character. baseline 날짜 변수명 (예: "accept.dt").
#' @param tx_col          character. 반응 날짜 변수명 (예: "tx.date").
#' @param fu_col          character. 추적 종료(사망·중도절단 발생 시점) 날짜 변수명 (예: "fu.date").
#' @param status_col      character. 사건 지표 변수명 (예: "fustat", 1=사망, 0=생존).
#' @param eps             0-길이 구간을 피하기 위해 아주 작은 시간을 더해 주는 보정값.
#' @param fix_zero_followup futime==0(accept=fu 같은 날) 처리 방식(bump: 아주 살짝 늘림, drop: 해당 사례 제외, error: 중단)
#'
#' @return data.frame. columns(id, tstart, tstop, event, response).

make_cp_from_dates <- function(data,
                               accept_col,
                               tx_col,
                               fu_col,
                               status_col,
                               eps = 1e-8,
                               fix_zero_followup = c("bump","drop","error")) {
  
  # --- (1) 입력 검증 및 id 생성 ---
  fix_zero_followup <- match.arg(fix_zero_followup) 
  
  # data.frame 여부
  stopifnot(is.data.frame(data))
  
  # 필수 열 존재 확인
  req <- c(accept_col, tx_col, fu_col, status_col)
  if (!all(req %in% names(data))) {
    stop("필수 열 누락: ", paste(setdiff(req, names(data)), collapse = ", "))
  }
  
  # id 없으면 생성(1..n)
  if (!("id" %in% names(data))) data$id <- seq_len(nrow(data))
  
  # 날짜형으로 변환
  acc <- as.Date(data[[accept_col]])
  tx  <- as.Date(data[[tx_col]])
  fu  <- as.Date(data[[fu_col]])
  
  # status → 정수(0/1)로 정규화
  st <- data[[status_col]]
  if (is.factor(st)) st <- as.integer(as.character(st))
  st <- as.integer(st)
  
  # --- (2) 날짜를 수치 시간으로 변환 ---
  futime <- as.numeric(fu - acc) # 추적기간
  txtime <- as.numeric(tx - acc)  # 반응까지의 시간 (NA 허용: 반응 없음)
  
  
  # --- (3) 데이터 품질검증 및 0-길이 구간처리 ---
  # 추적 종료가 baseline 이전이면 데이터 오류
  bad <- which(futime < 0)
  if (length(bad)) {
    stop("fu_col < accept_col 인 레코드가 있습니다. id: ",
         paste(head(data$id[bad], 5), collapse = ", "),
         if (length(bad) > 5) ", ...", "")
  }
  
  # 추적기간 0일(futime==0) 처리
  zero_fu <- which(futime == 0 | is.na(futime))
  if (length(zero_fu)) {
    if (fix_zero_followup == "bump") { 
      # 같은 날 종료 → 아주 작은 eps 만큼 늘려 0-길이 구간 방지
      futime[zero_fu] <- futime[zero_fu] + eps
    } else if (fix_zero_followup == "drop") { 
      # 해당 사례 제외
      keep <- setdiff(seq_len(nrow(data)), zero_fu)
      data <- data[keep, , drop = FALSE]
      acc  <- acc[keep]; tx <- tx[keep]; fu <- fu[keep]
      st   <- st[keep];  futime <- futime[keep]; txtime <- txtime[keep]
    } else { 
      # 에러메시지 출력 후 중단
      stop(sprintf("추적 0일 사례 %d건 error 발생", length(zero_fu)))
    }
  }
  
  # --- (4) 반응시점 조정(음수, 0, 또는 추적 종료 이후 처리) ---
  # baseline 이전 tx → 변화 없음으로 간주 (NA)
  txtime[!is.na(txtime) & txtime < 0] <- NA_real_
  
  # baseline 당일 tx → eps 만큼 뒤로(0-길이 분할 방지)
  txtime[!is.na(txtime) & txtime == 0] <- eps
  
  # tx 가 fu 와 같거나 이후(tx >= fu) → 추적구간 내 변화 없음으로 간주(NA)
  idx_ge <- which(!is.na(txtime) & !is.na(futime) & txtime >= futime)
  if (length(idx_ge)) {
    txtime[idx_ge] <- NA_real_
  }
  
  # --- (5) tmerge()를 이용한 구간 분할 ---
  base <- data.frame(id = data$id, tstart = 0, tstop = futime)
  
  # 사건(event): 개인별 tstop 시점에서 status==1 이면 사건
  tm <- tmerge(data1 = base, data2 = base, id = id,
               event = event(tstop, st == 1L))
  
  # 시간종속 이진 공변량(response): txtime 시점에 0→1로 전환
  tm <- tmerge(tm, base, id = id,
               response = tdc(txtime))
  
  # --- (6) 극소 구간 제거 및 반환 ---  
  # 0-길이 구간 제거
  tiny <- 1e-12
  tm <- tm[(tm$tstop - tm$tstart) > tiny, , drop = FALSE]
  
  if (any(tm$tstop <= tm$tstart)) {
    stop("여전히 0-길이 구간이 남아 있습니다. eps를 키우세요.")
  }
  
  # 반환
  tm[, c("id","tstart","tstop","event","response")]
}
