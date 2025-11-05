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
- 설치 예
  ```r
  install.packages(c("survival","dplyr","tibble","purrr","ggplot2","patchwork","survminer"))

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
  - 베이스라인, 반응(이식)일, 추적종료일 정보를 이용해 날짜 기반 자료를 Cox 모형용 counting-process(long) 형식으로 변환하고, 시간종속 이진 공변량(response)을 생성한다.
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
tstart, tstop, event, response 변수를 포함한 long-format data.frame을 반환한다.
이는 coxph(Surv(tstart, tstop, event) ~ response) 형태의 모형에 바로 사용 가능하다.


  ```r
  data.frame(id, tstart, tstop, event, response)

- 예시
  ```r
  cp <- make_cp_from_dates(
  data = dat,
  accept_col = "accept.dt",
  tx_col = "tx.date",
  fu_col = "fu.date",
  status_col = "fustat"
)


