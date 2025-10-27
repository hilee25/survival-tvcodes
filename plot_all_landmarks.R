# =======================
# Step 3 (plot 그리기)
# =======================
.survfit_at_L <- function(cp_df, L) {
  stopifnot(all(c("id","tstart","tstop","event","response") %in% names(cp_df)))
  had_event_before_L <- cp_df %>%
    group_by(id) %>% summarize(ev_before = any(event == 1 & tstop <= L), .groups="drop")
  at_risk_ids <- cp_df %>%
    group_by(id) %>% summarize(has_followup = any(tstop > L), .groups="drop")
  keep_ids <- had_event_before_L %>%
    inner_join(at_risk_ids, by = "id") %>% filter(!ev_before & has_followup) %>% pull(id)
  df_risk <- cp_df %>% filter(id %in% keep_ids)
  zL <- df_risk %>%
    filter(tstart <= L, tstop > L) %>% select(id, response) %>% distinct() %>% rename(Z_L = response)
  postL <- df_risk %>%
    filter(tstop > L) %>%
    mutate(tstart = pmax(tstart, L)) %>%
    filter(tstop > tstart) %>%
    left_join(zL, by = "id") %>%
    mutate(Z_L = factor(Z_L, levels = c(0,1), labels = c("Z=0","Z=1")))
  gtbl <- table(postL$Z_L)
  if (length(gtbl) < 2 || any(gtbl == 0)) {
    return(list(sf = NULL, postL = postL, ztab = gtbl))
  }
  sf <- survfit(Surv(tstart, tstop, event) ~ Z_L, data = postL)
  list(sf = sf, postL = postL, ztab = gtbl)
}

#' 분석결과 시각화 및 비교: 여러 Landmark를 한 번에
#'
#' @param cp_df data.frame(id, tstart, tstop, event, response)
#' @param landmarks numeric 또는 길이 1 이상의 벡터
#' @param robust 논문에서처럼 fit_for_L의 robust 옵션(기본 TRUE)
#' @param conf_int KM 곡선의 신뢰구간 표시 여부
#' @param tmax 모든 패널에서 공통 x축 상한(미지정 시 자동)
#' @param palette 라인 팔레트 (길이 2 권장)
#' @param annotate_p HR과 p-value 주석을 곡선 위에 표시할지 여부
#' @return list(plot = patchwork plot, table = tibble 요약표, per_L = 각 L의 객체 리스트)
plot_all_landmarks <- function(cp_df,
                               landmarks,
                               robust = TRUE,
                               conf_int = TRUE,
                               tmax = NULL,
                               palette = c("#1f77b4", "#d62728"),
                               annotate_p = TRUE,
                               theme_fn = ggplot2::theme_bw,
                               base_size = 12,
                               legend_position = "top",
                               ncol = 2,
                               title_prefix = "Landmark L = ",
                               show_titles = TRUE,
                               # ▼ 추가: 위험표 옵션
                               show_risktable = TRUE,
                               risk_table_height = 0.38,
                               risk_table_base_size = 10) {
  if (length(palette) < 2) palette <- rep_len(palette, 2)
  L_vec <- sort(unique(as.numeric(landmarks)))
  if (any(!is.finite(L_vec))) stop("landmarks must be finite numeric.")
  
  fits <- purrr::map(L_vec, ~tryCatch(fit_for_L(cp_df, .x, robust = robust),
                                      error = function(e) list(error = e, L = .x)))
  sfs  <- purrr::map(L_vec, ~.survfit_at_L(cp_df, .x))
  
  tab <- purrr::map2_dfr(fits, sfs, function(fit, sflist) {
    if (!is.null(fit$error)) {
      tibble(L = fit$L, n_total = NA_integer_, n_events = NA_integer_,
             n_Z0 = NA_integer_, n_Z1 = NA_integer_,
             HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
             note = paste0("fit error: ", conditionMessage(fit$error)))
    } else {
      z0 <- unname(if ("Z=0" %in% names(fit$n_group)) fit$n_group[["Z=0"]] else NA_integer_)
      z1 <- unname(if ("Z=1" %in% names(fit$n_group)) fit$n_group[["Z=1"]] else NA_integer_)
      note <- if (is.null(sflist$sf)) "KM unavailable (single/empty group)" else NA_character_
      tibble(L = fit$L, n_total = fit$n_total, n_events = fit$n_events,
             n_Z0 = z0, n_Z1 = z1,
             HR = fit$hr, CI_low = fit$ci_low, CI_high = fit$ci_high, p = fit$p,
             note = note)
    }
  }) %>% arrange(L)
  
  km_plots <- purrr::map2(L_vec, sfs, function(Li, sflist) {
    ttl <- paste0(title_prefix, format(round(Li, 3), trim = TRUE))
    if (is.null(sflist$sf)) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste0(ttl, "\nKM not available (group size issue)"),
                   size = 4) +
          theme_void() +
          labs(title = if (show_titles) ttl else NULL)
      )
    }
    
    g <- ggsurvplot(
      sflist$sf,
      data = sflist$postL,                
      conf.int = conf_int,
      risk.table = show_risktable,
      risk.table.col = "strata",
      risk.table.y.text = TRUE,           
      ggtheme = theme_fn(base_size = base_size),
      palette = palette,
      censor = TRUE,
      legend.title = NULL,      
      legend.labs = c("response=0", "response=1") 
    )
    
    p <- g$plot +
      labs(title = if (show_titles) ttl else NULL,
           x = "Time since L", y = "Survival") +
      theme(legend.position = legend_position,
            plot.title = element_text(face = "bold"),
            legend.title = element_blank()) + 
      guides(color   = guide_legend(title = NULL),  
             linetype= guide_legend(title = NULL),
             fill    = guide_legend(title = NULL))
    if (!is.null(tmax) && is.finite(tmax)) {
      p <- p + coord_cartesian(xlim = c(0, tmax))
    }
    if (annotate_p) {
      row <- tab %>% filter(L == Li)
      if (nrow(row) == 1 && is.finite(row$HR)) {
        lab <- sprintf("HR=%.2f (%.2f–%.2f),\nMantel-Byar test p=%s",
                       row$HR, row$CI_low, row$CI_high,
                       ifelse(row$p < 1e-4, "<1e-4",
                              formatC(row$p, format = "f", digits = 3)))
        p <- p + annotate("label", x = 0, y = 0, #위치 조정
                          label = lab, vjust = -0.2, hjust = 0,
                          size = 3, label.size = 0.2,
                          lineheight = 0.9)
      }
    }
    
    if (isTRUE(show_risktable) && !is.null(g$table)) {
      tbl <- g$table +
        theme_minimal(base_size = risk_table_base_size) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
        )
      p <- p / tbl + plot_layout(heights = c(1, risk_table_height))
    }
    p
  })
  
  grid_plot <- wrap_plots(km_plots, ncol = ncol)
  
  list(
    plot = grid_plot,
    table = tab,
    per_L = list(L = L_vec, fits = fits, survfits = sfs)
  )
}
