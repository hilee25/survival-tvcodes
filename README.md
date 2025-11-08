# 함수별 상세 사용 설명서 및 예제: survival-tvcodes

> 시간 종속 이진 공변량이 생존에 미치는 영향을 평가하기 위한 R 함수 세트  
> (랜드마크 분석, Mantel–Byar 검정, Cox 모형 적합 및 시각화)

---

## 환경

- **R**: 4.4.1.에서 작성 (2025.11월 기준)
- **필요한 패키지**
  - `survival` (Cox, Surv, tmerge)
  - `dplyr`, `tibble`, `purrr`
  - `ggplot2`, `patchwork`, `survminer` (시각화)
- 설치
  ```r
  install.packages(c("survival","dplyr","tibble","purrr","ggplot2","patchwork","survminer"))
  ```

## 함수 요약
| 함수명                    | 설명                                              |
| :--------------------- | :---------------------------------------------- |
| `make_cp_from_dates()` | Counting process 자료 변환 함수 |
| `fit_for_L()`          | 단일 랜드마크 분석 함수          |
| `.survfit_at_L()`      | L 이후 데이터 추출 + Kaplan–Meier 적합                   |
| `plot_all_landmarks()` | 다중 랜드마크 민감도 분석 및 시각화 함수     |
| `mantel_byar_cp()`     | Mantel-Byar 검정 함수            |

## 함수별 설명
**1) make_cp_from_dates()**

- 목적
  - 베이스라인, 반응일, 추적종료일 정보를 이용해 날짜 기반 자료를 counting-process 형식으로 변환하고, 시간종속 이진 공변량(response)을 생성한다.

- 사용법
  ```r
  make_cp_from_dates(
    data,
    accept_col,
    tx_col,
    fu_col,
    status_col,
    eps = 1e-8,
    fix_zero_followup = c("bump","drop","error")
  )
  ```

- 주요인자

| 인자                  | 설명                                                         |
| ------------------- | ---------------------------------------------------------- |
| `data`              | 원본 데이터프레임                                                  |
| `accept_col`        | 베이스라인(등록일) 변수명                                             |
| `tx_col`            | 반응(이식)일 변수명                                                |
| `fu_col`            | 추적 종료(사건 또는 검열 시점) 변수명                                     |
| `status_col`        | 사건 지표(1=사건, 0=검열)                                          |
| `eps`               | 0-길이 구간 방지용 아주 작은 시간값                                      |
| `fix_zero_followup` | 추적기간 0일(`fu=accept`) 처리 방식 (`"bump"`, `"drop"`, `"error"`) |


- 반환값
  - tstart, tstop, event, response 변수를 포함한 data.frame.
  - 이는 coxph(Surv(tstart, tstop, event) ~ response) 형태의 모형에 바로 사용 가능하다.
    ```r
    data.frame(id, tstart, tstop, event, response)
    ```

- 예시
  ```r
  cp <- make_cp_from_dates(
    data = dat,
    accept_col = "accept.dt",
    tx_col = "tx.date",
    fu_col = "fu.date",
    status_col = "fustat"
  )
  ```

**2) fit_for_L()**

- 목적
  - 단일 랜드마크 시점 L에서 위험집단을 정의하고, L 구간의 response 값을 고정한 그룹(Z_L)으로 Cox 비례위험 모형을 적합한다. Mantel–Byar 접근과 동일한 맥락의 score test p-value를 제공한다.

- 사용법
  ```r
  fit_for_L(cp_df, L, robust = TRUE)
  ```

- 주요인자

| 인자       | 설명                                                                    |
| -------- | --------------------------------------------------------------------- |
| `cp_df`  | `id, tstart, tstop, event, response`를 갖는 counting-process 데이터프레임 |
| `L`      | 랜드마크 시점                                                        |
| `robust` | `cluster(id)` 기반 강건 분산 사용 여부(기본 `TRUE`)                               |

- 반환값: 다음 원소를 갖는 리스트
  - L, n_total, n_events, n_group(Z=0/1 개수)
  - hr, ci_low, ci_high, p
  - fit (survival::coxph 객체)
 
- 예시
  ```r
  res_L120 <- fit_for_L(cp_df, L = 120, robust = TRUE)
  res_L120$hr; res_L120$p
  summary(res_L120$fit)
  ```


**3) .survfit_at_L()**

- 목적
  - 랜드마크 시점 L 이후 데이터(postL)를 구성하고, Z_L 그룹으로 Kaplan–Meier 생존곡선을 적합한다.

- 사용법
  ```r
  .survfit_at_L(cp_df, L)
  ```

- 주요인자

| 인자       | 설명                                                                    |
| -------- | --------------------------------------------------------------------- |
| `cp_df`  | `id, tstart, tstop, event, response`를 갖는 counting-process 데이터프레임 |
| `L`      | 랜드마크 시점                                                        |
  
- 반환값: 다음 원소를 갖는 리스트
  - sf: survfit 객체 (한쪽 그룹만 존재하면 NULL 반환)
  - postL: L 이후 데이터
  - ztab: id 기준 그룹 분포 table
 
- 예시
  ```r
  .survfit_at_L(cp_df, 30)
  ```

**4) plot_all_landmarks()**

- 목적
  - 여러 랜드마크 시점 L들에 대해, 각 L의 랜드마크 분석(HR/CI/p) 결과와 Kaplan–Meier 생존곡선, 요약표를 한 번에 비교·시각화한다.

- 사용법
  ```r
  plot_all_landmarks(
    cp_df = cp,
    landmarks = c(30, 60, 90, 120),
    robust = TRUE,
    conf_int = TRUE,
    tmax = NULL,
    palette = c("#1f77b4", "#d62728"),
    theme_fn = ggplot2::theme_bw,
    base_size = 12,
    legend_position = "top",
    ncol = 2,
    title_prefix = "Landmark L = ",
    show_titles = TRUE,
    show_risktable = TRUE,
    risk_table_height = 0.35,
    risk_table_base_size = 10,
    annotate_hrci = TRUE,
    annotate_mb_p = TRUE,
    annotate_cox_p = FALSE,
    ref_label = "Z=0"
  )
  ```
  
- 주요인자
  
| 인자                     | 설명                                                                                                                    |
| ---------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `cp_df`                | counting-process 형식 데이터 |
| `landmarks`            | 분석할 랜드마크 시점(원 데이터 시간 단위와 동일)                                                                   |
| `robust`               | `fit_for_L()`에서 `cluster(id)` 강건 분산 사용 여부, 추정치는 동일하고 표준오차/검정만 달라짐                                          |
| `conf_int`             | KM 곡선에 신뢰구간 밴드를 표시할지 여부                                                                 |
| `tmax`                 | 모든 패널의 공통 x축 상한으로, `NULL`이면 자동 스케일                                                                |
| `palette`              | 두 그룹(`Z=0`, `Z=1`)에 사용할 라인 색상 팔레트                                                                 |
| `theme_fn`             | KM 패널에 적용할 ggplot2 테마 함수                                                                         |
| `base_size`            | KM 패널의 기본 글자 크기(테마에 전달)                                                                                              |
| `legend_position`      | 범례 위치 (`"top"`, `"bottom"`, `"left"`, `"right"`, `"none"`)                                                           |
| `ncol`                 | patchwork 그리드 열 수                                                                                                |
| `title_prefix`         | 각 패널 제목의 접두어                                                                                                         |
| `show_titles`          | 각 패널의 제목을 표시할지 여부                                                                                                    |
| `show_risktable`       | 패널 하단에 위험표(risk table)를 붙일지 여부                                                                      |
| `risk_table_height`    | 위험표가 차지하는 상대 높이                                                                    |
| `risk_table_base_size` | 위험표 텍스트 크기                                                                                                           |
| `annotate_hrci`        | 패널 내부에 **HR과 95% CI** 주석을 표시할지 여부                                                                                    |
| `annotate_mb_p`        | 패널 내부에 **Mantel–Byar p** 주석을 표시할지 여부                                                        |
| `annotate_cox_p`       | 패널 내부에 **Cox score p** 주석을 표시할지 여부                                                                                   |
| `ref_label`            | Mantel–Byar 계산 시 참조 그룹 레이블으로, 내부에서 `response = as.integer(Z_L != ref_label)`로 변환하여 참조를 0으로 둔다                   |

  
- 반환값: 다음 원소를 갖는 리스트
  - plot: 랜드마크별 KM 패널을 묶은 patchwork 객체
  - table: 랜드마크별 요약표(L, n_total, n_events, n_Z0, n_Z1, HR, CI_low, CI_high, p, note)
  - per_L: 내부 객체 목록(fit_for_L 결과와 .survfit_at_L 결과)

- 예시
  ```r
  # 기본(HR/CI + Mantel–Byar p)
  res1 <- plot_all_landmarks(cp, landmarks = c(30, 60, 90, 120), conf_int = FALSE)
  print(res1$plot)    # patchwork grid
  print(res1$table)   # L별 요약표


  # Mantel–Byar & Cox p 함께 표시
  res2 <- plot_all_landmarks(
    cp, landmarks = c(30, 60, 90, 120),
    annotate_cox_p = TRUE,
    conf_int = FALSE
  )
  print(res2$plot)

  # HR/CI + Cox p만 표시
  res3 <- plot_all_landmarks(
    cp, landmarks = c(30, 60, 90, 120),
    annotate_mb_p = FALSE,
    annotate_cox_p = TRUE,
    conf_int = FALSE
  )
  print(res3$plot)

  # 주석 전부 끄기
  res4 <- plot_all_landmarks(
    cp, landmarks = c(30, 60, 90, 120),
    annotate_hrci = FALSE, annotate_mb_p = FALSE, annotate_cox_p = FALSE,
    conf_int = FALSE
  )
  print(res4$plot)

  # 테이블 미표시
  res5 <- plot_all_landmarks(
    cp, landmarks = c(30, 60, 90, 120),
    show_risktable = FALSE,   
    risk_table_base_size = 9,
    conf_int = FALSE
  )
  print(res5$plot)
  ```


**5) mantel_byar_cp()**

- 목적
  - 시간종속 이진 공변량이 있는 counting-process 형식 데이터를 이용하여 Mantel–Byar 검정을 수행한다.
  - 각 사건 시점마다 위험집단을 재구성하고, 참조군(ref_val) 대비 관측·기대 사건 수를 비교해 로그순위(log-rank)형 검정통계량을 계산한다.

- 사용법
  ```r
  mantel_byar_cp(cp, ref_val = 0)
  ```
- 주요인자

| 인자            | 설명                                                                                                                                                                 |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `cp`      | counting-process 형식 데이터                                                     |
| `ref_val` | 기준(참조) 상태 값으로, 일반적으로 0을 참조군으로 두며, `response`의 이진값(0/1)에 대해 비교를 수행                                         |

  
- 반환값: Mantel–Byar 검정 결과를 담는, 다음 원소를 갖는 리스트
  - chisq: Mantel–Byar 검정통계량 
  - p: 유의확률(p-value), 자유도 1의 카이제곱 분포 기반
  - D_ref: 참조군(`ref_val`)에서 관측된 사건 수의 총합
  - E: 참조군에서 기대된 사건 수의 총합
  - V: 분산의 총합
  - contrib: 각 사건 시점별 기여도 테이블(time, d, N, Nref, dref, OE, V)

- 예시
  ```r
  mb_res <- mantel_byar_cp(cp, ref_val = 0)
  mb_res$chisq
  mb_res$p
  ```
