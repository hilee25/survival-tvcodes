# 필요한 패키지
library(survival); library(dplyr); library(ggplot2)
library(purrr); library(tibble); library(patchwork); library(survminer)
options(scipen = 999) # 지수표기 억제

# =========================================================
# Step 1) tmerge: 날짜 열을 counting-process(long) 형식으로 변환
# =========================================================

#' @param data            분석할 data.frame
#' @param landmark_days   numeric 벡터. baseline(accept)로부터 며칠 후를 landmark로 지정 (예: c(60, 120)).
#' @param accept_col      character. baseline 날짜 변수명 (예: "accept.dt").
#' @param tx_col          character. 반응(이식) 날짜 변수명 (예: "tx.date").
#' @param fu_col          character. 추적 종료(사망·중도절단 발생 시점) 날짜 변수명 (예: "fu.date").
#' @param status_col      character. 사건 지표 변수명 (예: "fustat", 1=사망, 0=생존).
#' @param eps             0-길이 구간을 피하기 위해 아주 작은 시간을 더해 주는 보정값
#' @param fix_zero_followup futime==0(accept=fu 같은 날) 처리 방식(bump: 아주 살짝 늘림, drop: 해당 사례 제외, error: 중단)

make_cp_from_dates <- function(data,
                               accept_col,
                               tx_col,
                               fu_col,
                               status_col,
                               eps = 1e-8,
                               fix_zero_followup = c("bump","drop","error")) {
  
  fix_zero_followup <- match.arg(fix_zero_followup) 
  
  # --- 유효성 점검 ---
  stopifnot(is.data.frame(data))
  req <- c(accept_col, tx_col, fu_col, status_col)
  if (!all(req %in% names(data))) {
    stop("필수 열 누락: ", paste(setdiff(req, names(data)), collapse = ", "))
  }
  
  # --- id가 없으면 만들기(식별) ---
  if (!("id" %in% names(data))) data$id <- seq_len(nrow(data))
  
  # --- 날짜형으로 변환 ---
  acc <- as.Date(data[[accept_col]])
  tx  <- as.Date(data[[tx_col]])
  fu  <- as.Date(data[[fu_col]])
  
  st <- data[[status_col]]
  if (is.factor(st)) st <- as.integer(as.character(st))
  st <- as.integer(st)
  
  # 일(day) 단위 시간축 만들기
  futime <- as.numeric(fu - acc)
  txtime <- as.numeric(tx - acc)  # NA 허용=반응없음
  
  # --- 추적종료가 베이스라인보다 과거면 데이터 오류로 간주 ---
  # fu < accept이면 데이터 오류
  bad <- which(futime < 0)
  if (length(bad)) {
    stop("fu_col < accept_col 인 레코드가 있습니다. id: ",
         paste(head(data$id[bad], 5), collapse = ", "),
         if (length(bad) > 5) ", ...", "")
  }
  
  # --- 구간 길이 0 문제가 생기지 않도록 처리 ---
  # fu=acc(0일 구간) 처리 
  zero_fu <- which(futime == 0 | is.na(futime))
  if (length(zero_fu)) {
    if (fix_zero_followup == "bump") { # bump: 아주 작게 늘림
      futime[zero_fu] <- futime[zero_fu] + eps
    } else if (fix_zero_followup == "drop") { # drop: 제외(표본 크기 감소)
      keep <- setdiff(seq_len(nrow(data)), zero_fu)
      data <- data[keep, , drop = FALSE]
      acc  <- acc[keep]; tx <- tx[keep]; fu <- fu[keep]
      st   <- st[keep];  futime <- futime[keep]; txtime <- txtime[keep]
    } else { # error: 중단
      stop(sprintf("추적 0일 사례 %d건 error 발생", length(zero_fu)))
    }
  }
  
  # tx < accept → 이식 없음으로 간주(NA) 
  txtime[!is.na(txtime) & txtime < 0] <- NA_real_
  
  # accept 당일 이식이면 eps만큼 뒤로
  txtime[!is.na(txtime) & txtime == 0] <- eps
  
  # 이식시점이 추적종료와 같거나 이후 일 때(tx >= fu) NA 처리
  idx_ge <- which(!is.na(txtime) & !is.na(futime) & txtime >= futime)
  if (length(idx_ge)) {
    txtime[idx_ge] <- NA_real_        # post-tx 시간 0이므로 이식없음으로 간주
  }
  
  # tmerge: 자료 생성
  base <- data.frame(id = data$id, tstart = 0, tstop = futime)
  
  tm <- tmerge(data1 = base, data2 = base, id = id,
               event = event(tstop, st == 1L)) # 사건 시간 지정(개인별 tstop 시점)
  tm <- tmerge(tm, base, id = id,
               response = tdc(txtime)) # time-varying binary covariate
  
  # 남은 0길이 구간 제거
  tiny <- 1e-12
  tm <- tm[(tm$tstop - tm$tstart) > tiny, , drop = FALSE]
  
  # check
  if (any(tm$tstop <= tm$tstart)) {
    stop("여전히 0-길이 구간이 남아 있습니다. eps를 키우세요.")
  }
  
  tm[, c("id","tstart","tstop","event","response")]
}

