#' @import patchwork
utils::globalVariables(c("fstrata", "lower", "n_risk", "surv", "strata", "upper"))

ggKM.BPCP <- function(data_input, method) {
  if (!requireNamespace("bpcp", quietly = TRUE)) stop("[ggKM.BPCP] requires package 'bpcp'")
  factors <- unique(data_input$strata)
  output <- list()
  for (f in factors) {
    subdata_input <- data_input[data_input$strata == f,]
    fit <- bpcp::bpcpfit(subdata_input$time, subdata_input$status)
    output_data <- data.frame("time" = fit[[1]]$L,
                              "surv" = fit[[1]]$surv,
                              "lower" = fit[[1]]$lower,
                              "upper" = fit[[1]]$upper,
                              "strata" = f) 
    output[[as.character(f)]] <- output_data
  }
  return(do.call(rbind, output))
}

ggKM.censor <- function(data_summary, data_input) {
  data_input <- data_input[data_input$status == 0,]
  factors <- unique(data_input$strata)
  output <- list()
  for (f in factors) {
    subdata_input <- data_input[data_input$strata == f,]
    subdata_summary <- data_summary[data_summary$strata == f,]
    if (subdata_summary$time[1] != 0) {
      startrow <- data.frame("time" = 0, "surv" = 1, "lower" = 1, "upper" = 1, "strata" = f)
      subdata_summary <- rbind(startrow, subdata_summary)
    }
    output[[as.character(f)]] <- data.frame(
      "time" = subdata_input$time,
      "surv" = subdata_summary$surv[findInterval(subdata_input$time, subdata_summary$time)],
      "strata" = f)
  }
  return(do.call(rbind, output))
}

ggKM.custom <- function(data_input, f_custom) {
  factors <- unique(data_input$strata)
  output <- list()
  for (f in factors) {
    sub_data <- data_input[data_input$strata == f,]
    .time <- sub_data$time
    .status <- sub_data$status
    .weights <- sub_data$weights
    environment(f_custom) <- environment()
    output.custom <- f_custom()
    #output.custom$strata <- f
    #output[[as.character(f)]] <- output.custom
    # output.custom <- f_custom()
    output[[as.character(f)]] <- data.frame(
      "time" = output.custom$time,
      "surv" = output.custom$surv,
      "lower" = output.custom$lower,
      "upper" = output.custom$upper,
      "strata" = f)
  }
  return(do.call(rbind, output))
}

ggKM.LOCF <- function(x) {
  n <- length(x)
  if (n == 0) return(x)
  if (is.na(x[1])) return(x)
  for (i in 2:n) {
    if (is.na(x[i])) {
      x[i] <- x[i - 1]
      break
    }
  }
  return(x)
}

ggKM.step <- function(data_summary, data_input, t_max = Inf) {
  factors <- unique(data_summary$strata)
  output <- list()
  i <- 0
  for (f in factors) {
    subdata_input <- data_input[data_input$strata == f,]
    subdata_summary <- data_summary[data_summary$strata == f,]
    # deal with NA confidence limits when surv == 0
    subdata_summary$lower[is.na(subdata_summary$lower) & subdata_summary$surv == 0] <- 0
    subdata_summary$upper[is.na(subdata_summary$upper) & subdata_summary$surv == 0] <- 0
    if (subdata_summary$time[1] != 0) {
      startrow <- data.frame("time" = 0, "surv" = 1, "lower" = 1, "upper" = 1, "strata" = f)
      subdata_summary <- rbind(startrow, subdata_summary)
    }
    end <- 2 * nrow(subdata_summary)
    newdata <- data.frame(
      "time" = rep(subdata_summary$time, each = 2)[-1],
      "surv" = rep(subdata_summary$surv, each = 2)[-end],
      "lower" = rep(subdata_summary$lower, each = 2)[-end],
      "upper" = rep(subdata_summary$upper, each = 2)[-end],
      "strata" = f)
    newdata <- newdata[newdata$time <= t_max,]
    subdata_t_max <- min(max(subdata_input$time), t_max)
    if (max(newdata$time) < subdata_t_max) {
      endrow <- newdata[nrow(newdata),]
      endrow$time <- subdata_t_max
      output[[as.character(f)]] <- rbind(newdata, endrow)
    } else {
      output[[as.character(f)]] <- newdata
    }
  }
  return(do.call(rbind, output))
}

ggKM.WH <- function(data_input, method) {
  if (!requireNamespace("WHKMconf", quietly = TRUE)) stop("[ggKM.WH] requires package 'WHKMconf'")
  factors <- unique(data_input$strata)
  output <- list()
  for (f in factors) {
    subdata_input <- data_input[data_input$strata == f,]
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = subdata_input)
    summary <- summary(fit)
    if (method == "Rothman") {
      CI <- WHKMconf::WH_Rothman(summary$surv, summary$n.risk, summary$n.event)
    } else if (method == "TG") {
      CI <- WHKMconf::WH_ThomasGrunkemeier(summary$time, summary$n.risk, summary$n.event)
    } else if (method == "Nair") {
      CI <- WHKMconf::WH_Nair(summary$time, summary$surv, summary$std.err, summary$n.risk, summary$n.event)
    } else if (method == "HM") {
      CI <- WHKMconf::WH_HollanderMcKeague(summary$time, summary$n.risk, summary$n.event, adapt = FALSE)
    } else stop("[ggKM.WH] invalid method")
    output[[as.character(f)]] <- data.frame(
      "time" = summary$time,
      "surv" = summary$surv,
      "lower" = CI[, 2],
      "upper" = CI[, 3],
      "strata" = f)
  }
  return(do.call(rbind, output))
}

################################################################################

#' Kaplan–Meier plot
#'
#' Generates a Kaplan–Meier plot with optional confidence intervals and risk 
#' table.
#' @param time Numeric vector of follow-up times.
#' @param status Integer vector event indicator (`1` = event, `0` = censored).
#' @param group Optional integer vector grouping variable.
#' @param breaks.s Y-axis (survival) tick marks. Default = `seq(0, 1, 0.25)`.
#' @param breaks.t X-axis (time) tick marks. Default = 
#' `seq(0, max(time), by = 12)`.
#' @param CI String indicating the CI type (or a custom function; see Details):
#'   \itemize{
#'     \item `"none"`: none
#'     \item Pointwise confidence intervals:
#'     \itemize{
#'       \item `"cloglog"` (default): Greenwood's variance¹ with complementary 
#'       log–log transformation² 
#'       \item `"modcloglog"`: `"cloglog"` with the lower confidence limit 
#'       modified according to the effective sample size at each censored 
#'       observation³
#'       \item `"Rothman"`: Rothman's binomial method⁴ (via the `WHKMconf` 
#'       package)
#'       \item `"TG"`: Thomas–Grunkemeier likelihood-ratio method⁵ (via the 
#'       `WHKMconf` package)
#'       \item `"BPCP"`: beta product confidence procedure⁶ (via the `bpcp` 
#'       package)
#'     }
#'     \item Simultaneous confidence bands:
#'     \itemize{
#'       \item `"Nair"`: Nair's log-transformed equal precision method⁷ (via the 
#'       `WHKMconf` package)
#'       \item `"HM"`: Hollander–McKeague likelihood-ratio method⁸ (via the 
#'       `WHKMconf` package)
#'     }
#'   }
#' @param CI.alpha Alpha transparency of the confidence intervals. Default = 
#' `0.2`.
#' @param colors Vector of colors of the survival curves. Default = `c("#0087BD", 
#' "#C40233", "#009F6B", "#FEA500", "#9755AE", "#33B1BD", "#FFD300", "#D8668B")`.
#' @param digits.fixed Number of decimal places for the risk table if `weights` 
#' are used. Default = `1`.
#' @param grid.color Gridline color. Default = `grDevices::rgb(0.95, 0.95, 0.95)`.
#' @param grid.s Horizontal gridline positions. Default = `seq(0, 1, 0.25)`.
#' @param grid.t Vertical gridline positions. Default = `NULL`.
#' @param grid.width Gridline thickness. Default = `0.5`.
#' @param legend.direction Legend orientation; either `"vertical"` (default) or 
#' `"horizontal"`.
#' @param legend.just Alignment anchor for the legend relative to its position. 
#' Can be a single keyword (e.g., `"center"`), a keyword pair specifying 
#' horizontal and vertical justification (e.g., `c("left", "top")`, 
#' `c("right", "bottom")`) respectively, or a numeric vector of length 2 giving 
#' relative coordinates within the plot area. Default = `c("right", "top")`.
#' @param legend.label.pos Legend label position (`"left"` or `"right"`) 
#' relative to the legend symbol. Default = `"left"`.
#' @param legend.labels Character vector of group labels.
#' @param legend.ncol Integer specifying the number of columns in the legend. 
#' Default = `NULL`.
#' @param legend.nrow Integer specifying the number of rows in the legend. 
#' Default = `NULL`.
#' @param legend.pos Position of the legend. Can be a keyword such as `"none"`, 
#' `"left"`, `"right"`, `"bottom"`, or `"top"`, or a numeric vector of length 2 
#' giving relative coordinates within the plot area. Default = `c(1, 1)`. 
#' Set to `"none"` if `group` is `NULL`.
#' @param legend.text.align Legend text alignment: `0` = left, `0.5` = center, 
#' `1` = right. Default = `1`.
#' @param line.width Line width survival curves and censor marks. Default = `0.5`.
#' @param line.height Height of censor marks. Default = `0.025`.
#' @param margin.KM.b Numeric value specifying the size of the Kaplan–Meier plot
#' bottom margin in points. Default = `0`.
#' @param margin.KM.l Numeric value specifying the size of the Kaplan–Meier plot
#' left margin in points. Default = `0`.
#' @param margin.KM.r Numeric value specifying the size of the Kaplan–Meier plot
#' right margin in points. Default = `0`.
#' @param margin.KM.t Numeric value specifying the size of the Kaplan–Meier plot
#' top margin in points. Default = `4`.
#' @param margin.risk Numeric value specifying the horizontal spacing between 
#' the risk table labels and the risk table. Default = `16`. 
#' @param max.t Optional numeric value specifying the maximum time up to which 
#' the Kaplan–Meier plot and risk table are displayed. 
#' @param risk.table Integer specifying the risk table mode:
#' \itemize{
#'   \item `0`: No risk table 
#'   \item `1`: Labelled risk table
#'   \item `2` (default): Color-coded risk table 
#' }
#' @param risk.table.prop Relative height of the risk table. Default = `0.1`.
#' @param t.labels Optional string vector specifying custom labels for the 
#' time-axis breaks defined by `breaks.t`.
#' @param textsize.axis Axis text size. Default = `12`.
#' @param textsize.legend Legend text size. Default = `12`.
#' @param textsize.risk Risk table text size. Default = `12`.
#' @param title.s Y-axis (survival) title. Default = `"Survival"`.
#' @param title.t X-axis (time) title. Default = `"Time"`.
#' @param weights Optional numeric vector of observation weights.
#' @return A `ggplot` (+ `patchwork` if `risk.table` ≥ `1`) object. The returned 
#' object also stores the processed input `time`, `status`, and `group` as attributes.
#' 
#' @return A `ggplot` object (combined with `patchwork` if `risk.table ≥ 1`). 
#' The returned object stores additional data as attributes:
#' \itemize{
#'   \item Processed input data:
#'   \itemize{
#'     \item `attr(x, "time")`: follow-up times
#'     \item `attr(x, "status")`: event indicator
#'     \item `attr(x, "group")`: grouping variable
#'   }
#'   \item Kaplan–Meier estimates (truncated at `max.t`):
#'   \itemize{
#'     \item `attr(x, "KM_time")`: step times
#'     \item `attr(x, "KM_surv")`: survival estimates
#'     \item `attr(x, "KM_lower")`: lower confidence limits
#'     \item `attr(x, "KM_upper")`: upper confidence limits
#'     \item `attr(x, "KM_group")`: grouping variable
#'   }
#' }
#' @details
#' Observations with invalid `time` or missing `status` are automatically 
#' excluded.
#' 
#' A custom function passed to `CI` is called once per `group` in an environment 
#' containing the following variables:
#' \itemize{
#'   \item `.time`: Numeric vector of follow-up times for the current `group` 
#'   \item `.status`: Integer vector event indicator (`1` = event, `0` = 
#'   censored) for the current `group`
#'   \item `.weights`: Numeric vector of weights for the current `group` 
#' }
#' The custom function must return a `data.frame` with the following columns:
#' \itemize{
#'   \item `$time`: times
#'   \item `$surv`: survival estimate 
#'   \item `$lower`: lower confidence limit
#'   \item `$upper`: upper confidence limit
#' }
#' @references
#' 1. Greenwood, M., 1926. A report on the natural duration of cancer. In: 
#' *Reports on Public Health and Medical Subjects*, 33, pp. 1–26. London: 
#' Her Majesty’s Stationery Office, Ministry of Health.
#' 2. Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K., 2007. Analyzing 
#' survival curves at a fixed point in time. *Statistics in Medicine*, 26(24), 
#' pp. 4505–4519.
#' 3. Dorey, F.J. and Korn, E.L., 1987. Effective sample sizes for confidence 
#' intervals for survival probabilities. *Statistics in Medicine*, 6(6), pp. 
#' 679–687.
#' 4. Rothman, K.J., 1978. Estimation of confidence limits for the cumulative 
#' probability of survival in life table analysis. *Journal of Chronic Diseases*,
#' 31(8), pp. 557–560.
#' 5. Thomas, D.R. and Grunkemeier, G.L., 1975. Confidence interval estimation 
#' of survival probabilities for censored data. *Journal of the American 
#' Statistical Association*, 70(352), pp. 865–871.
#' 6. Fay, M.P., Brittain, E.H. and Proschan, M.A., 2013. Pointwise confidence 
#' intervals for a survival distribution with small samples or heavy censoring. 
#' *Biostatistics*, 14(4), pp. 723–736.
#' 7. Nair, V.N., 1984. Conﬁdence bands for survival functions with censored
#' data: a comparative study. *Technometrics*, 26, pp. 265–275.
#' 8. Hollander, M. and McKeague, I.W., 1997. Likelihood ratio-based confidence 
#' bands for survival functions. *Journal of the American Statistical 
#' Association*, 92(437), pp. 215–226.
#' @examples
#' data <- survival::lung
#' fig1a <- ggKM(data$time * 12 / 365.2425, data$status - 1, data$sex, 
#'               breaks.t = seq(0, 30, 6), legend.labels = c("Male", "Female"),
#'               title.s = "Overall survival", title.t = "Time (months)")
#' print(fig1a)
#' f.custom <- function() {
#'   fit <- survival::survfit(survival::Surv(.time, .status) ~ 1,  
#'                            conf.type = "arcsin")
#'   return(data.frame("time" = fit$time, "surv" = fit$surv, 
#'                     "lower" = fit$lower, "upper" = fit$upper))
#' }
#' fig1b <- ggKM(data$time * 12 / 365.2425, data$status - 1, data$sex, 
#'               breaks.t = seq(0, 30, 6), legend.labels = c("Male", "Female"), 
#'               title.s = "Overall survival", title.t = "Time (months)", 
#'               CI = f.custom)
#' print(fig1b)
#' @export
ggKM <- function(time, status, group = NULL,
                 breaks.s = seq(0, 1, 0.25), breaks.t = NULL,
                 CI = "cloglog", CI.alpha = 0.2,
                 colors = c("#0087BD", "#C40233", "#009F6B", "#FEA500", 
                            "#9755AE", "#33B1BD", "#FFD300", "#D8668B"),
                 digits.fixed = 1,
                 grid.color = grDevices::rgb(0.95, 0.95, 0.95),
                 grid.s = seq(0, 1, 0.25), grid.t = NULL, grid.width = 0.5,
                 legend.direction = "vertical", legend.just = c("right", "top"),
                 legend.label.pos = "left", 
                 legend.labels = NULL, legend.ncol = NULL, legend.nrow = NULL,
                 legend.pos = c(1, 1), legend.text.align = 1,
                 margin.KM.b = 0, margin.KM.l = 0, margin.KM.r = 0, margin.KM.t = 4,
                 margin.risk = 16, line.width = 0.5, line.height = 0.025, max.t = NULL,
                 risk.table = 2, risk.table.prop = 0.1, t.labels = NULL,
                 textsize.axis = 12, textsize.legend = 12, textsize.risk = 12,
                 title.s = "Survival", title.t = "Time", weights = NULL) {
  keep <- !is.na(time) & !is.na(status) & time > 0
  colors <- c(colors, "#00000000") # need an extra dummy color
  if (is.null(weights)) {
    w <- rep(1, length(time))
    digits.fixed <- 0
  } else {
    w <- weights
  }
  time <- time[keep]
  status <- status[keep]
  if (!is.null(group)) group <- group[keep]
  w <- w[keep]
  if (inherits(CI, "function")) {
    f_custom <- CI
    CI <- "custom"
  }
  if (is.null(group)) group <- rep(1, length(time))
  group <- as.numeric(group)
  n_group <- length(unique(group))
  if (n_group == 1) legend.pos = "none"
  if (is.null(breaks.t)) breaks.t <- seq(0, max(time), by = 12)
  if (is.null(legend.labels)) legend.labels <- 0:(n_group - 1)
  if (is.null(t.labels)) t.labels <- breaks.t
  data_input <- data.frame("time" = time, "status" = status, "strata" = group, "w" = w)
  data_censor <- data_input[data_input$status == 0,]
  if (F) {
    fit <- survival::survfit(survival::Surv(time, status) ~ strata, 
                             data = data_input, conf.type = "log-log", 
                             conf.lower = "modified", weights = w)
    s <- data.frame("time" = fit$time, "surv" = fit$surv, "lower" = fit$lower, "upper" = fit$upper)
    if (n_group == 1) {
      s$strata <- 1
    } else {
      x <- as.numeric(fit$strata)
      s$strata <- rep(seq_along(x), times = x)
    }
  } else {
    fit <- survival::survfit(survival::Surv(time, status) ~ strata, 
                             data = data_input, conf.type = "log-log", 
                             conf.lower = "usual", weights = w)
    s <- summary(fit)
    if (n_group == 1) {
      s$strata <- rep(1, length(s$time))
    } else {
      s$strata <- as.numeric(sub("strata=", "", s$strata))
    }
  }
  if (CI == "cloglog" | CI == "modcloglog" | CI == "none") {
    data_summary <- data.frame(s[c("time", "surv", "lower", "upper", "strata")])
  } else if (CI == "HM" | CI == "Nair" | CI == "TG" | CI == "Rothman") {
    if (!is.null(weights)) stop(paste0("[ggKM] weights not implemented for CI = ", CI))
    data_summary <- ggKM.WH(data_input, CI)
  } else if (CI == "BPCP") {
    if (!is.null(weights)) stop(paste0("[ggKM] weights not implemented for CI = ", CI))
    data_summary <- ggKM.BPCP(data_input)
  } else if (CI == "custom") {
    data_summary <- ggKM.custom(data_input, f_custom)
  } else stop("[ggKM] invalid CI")
  if (is.null(max.t)) max.t <- max(time)
  x_lim <- c(0, max.t)
  data_step <- ggKM.step(data_summary, data_input, max.t)
  levels_order <- sort(unique(group)) # enforce factor level order consistency
  data_step$fstrata <- factor(data_step$strata, levels = levels_order)
  data_pruned <- data_step[!is.na(data_step$lower),]
  g_KM <- ggplot2::ggplot(data_step, ggplot2::aes(time, surv, color = fstrata, fill = fstrata)) +
    ggplot2::coord_cartesian(ylim = c(0, 1), xlim = x_lim, clip = "off") +
    ggplot2::geom_hline(yintercept = grid.s, color = grid.color, linewidth = grid.width) +
    ggplot2::geom_vline(xintercept = grid.t, color = grid.color, linewidth = grid.width) +
    ggplot2::guides(color = ggplot2::guide_legend(label.position = legend.label.pos, ncol = legend.ncol, nrow = legend.nrow),
                    fill  = ggplot2::guide_legend(label.position = legend.label.pos, ncol = legend.ncol, nrow = legend.nrow)) +
    ggplot2::labs(x = title.t, y = title.s, color = NULL) +
    ggplot2::scale_color_manual(
      name = "",
      values = colors,
      labels = legend.labels,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::scale_x_continuous(breaks = breaks.t, labels = t.labels, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = breaks.s, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = textsize.axis),
                   axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 1)),
                   axis.title.x = ggplot2::element_text(size = textsize.axis),
                   legend.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.box.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.direction = legend.direction, 
                   legend.justification = legend.just,
                   legend.position = legend.pos,
                   legend.text = ggplot2::element_text(size = textsize.legend),
                   legend.text.align = legend.text.align, 
                   legend.title = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = margin.KM.t, r = margin.KM.r, 
                                                 b = margin.KM.b, l = margin.KM.l))
  if (CI != "none") {
    g_KM <- g_KM +
      ggplot2::geom_ribbon(data = data_pruned, 
                           ggplot2::aes(ymin = lower, ymax = upper, fill = fstrata),
                           alpha = CI.alpha, colour = NA)
  }
  data_censor <- ggKM.censor(data_summary, data_input)
  data_censor$fstrata <- factor(data_censor$strata, levels = levels_order)
  if (!is.null(data_censor$time)) {
    g_KM <- g_KM + ggplot2::geom_segment(
      data = data_censor[data_censor$time <= max.t,],
      ggplot2::aes(x = time, xend = time,
                   y = surv - line.height/2, yend = surv + line.height/2,
                   color = fstrata),
      linewidth = line.width
    )
  }
  g_KM <- g_KM + ggplot2::geom_line(linewidth = line.width)
  if (risk.table >= 1) { 
    s <- summary(fit, times = breaks.t)
    data_risk <- data.frame(
      time = s$time,
      strata = NA,
      n_risk = s$n.risk
    )
    data_risk$n_risk <- formatC(data_risk$n_risk, format = "f", digits = digits.fixed)
    if (risk.table == 1) {
      if (n_group > 1) {
        data_risk$strata <- as.factor(as.numeric(sub("strata=", "", s$strata)))      
        data_risk <- rbind(data.frame(time = 0, strata = "At risk:", n_risk = ""), data_risk)
        data_risk$strata <- as.factor(data_risk$strata)
        data_risk$strata <- stats::relevel(data_risk$strata, "At risk:")
        levels(data_risk$strata) <- c("At risk:", legend.labels)
      } else {
        data_risk$strata <- "At risk:"
      }
      g_risk <- ggplot2::ggplot(data_risk[data_risk$time <= max.t,], 
                                ggplot2::aes(time, strata, label = n_risk)) +
        ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
        ggplot2::geom_text(size = textsize.risk * 127/360) +
        ggplot2::scale_x_continuous(breaks = breaks.t, expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          axis.line = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = margin.risk), 
                                              size = textsize.risk),
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank()
        )
    } else if (risk.table == 2) {
      if (n_group > 1) {
        data_risk$strata <- factor(as.numeric(sub("strata=", "", s$strata)), levels = levels_order) 
        colors_risk <- colors[seq_len(n_group)]
      } else {
        stop("[ggKM] risk.table == 2 not implemented for n_group == 1")
      }
      g_risk <- ggplot2::ggplot(
        data_risk[data_risk$time <= max.t,],
        ggplot2::aes(time, strata, label = n_risk)) +
        ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
        ggplot2::geom_text(ggplot2::aes(color = strata), size = textsize.risk * 127/360) +
        ggplot2::scale_color_manual(values = colors_risk, guide = "none") +
        ggplot2::scale_x_continuous(breaks = breaks.t, expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev(levels(data_risk$strata))) +
        ggplot2::theme_classic() + ggplot2::ggtitle(" ") +
        ggplot2::theme(
          axis.line   = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks  = ggplot2::element_blank(),
          axis.title  = ggplot2::element_blank(),
          panel.grid  = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0, size = textsize.axis, 
                                             margin = ggplot2::margin(b = 0)),
          plot.title.position = "plot"
        ) + ggplot2::labs(y = NULL) +
        ggplot2::annotation_custom(
          grid::textGrob(
            "At risk", rot = 90, x = grid::unit(-2.75, "lines"),
            gp = grid::gpar(fontsize = textsize.axis)
          )
        )
    } else stop("[ggKM] invalid risk.table")
    if (risk.table == 2) {
      y_lab <- "" # dummy title for spacing
    } else {
      y_lab <- NULL
    }
    g_KM <- g_KM + ggplot2::labs(y = y_lab) + 
      ggplot2::annotation_custom(
        grid::textGrob(
          title.s, rot = 90, x = grid::unit(-2.75, "lines"),
          gp = grid::gpar(fontsize = textsize.axis)
        )
      )
    g <- g_KM / g_risk +
      patchwork::plot_layout(heights = c(((1 - risk.table.prop) / risk.table.prop), 1))
  } else {
    g <- g_KM
  }
  attr(g, "time") <- time
  attr(g, "status") <- status
  attr(g, "group") <- group
  attr(g, "KM_time") <- data_pruned$time
  attr(g, "KM_surv") <- data_pruned$surv
  attr(g, "KM_lower") <- data_pruned$lower
  attr(g, "KM_upper") <- data_pruned$upper
  attr(g, "KM_group") <- data_pruned$strata 
  return(g)
}

#' Add extra label to a ggKM legend
#'
#' Adds an extra text label to an existing legend of a `ggKM` plot.
#' @param input A `ggKM` plot object.
#' @param label String to append as an extra legend label.
#' @return A deep copy of the input `ggKM` object with the added legend label.
#' @export
ggKM.legend.extra <- function(input, label) {
  # R CMD check behavior
  time <- surv <- fstrata <- NULL
  output <- unserialize(serialize(input, NULL))
  output[[1]]$scales$scales[[1]]$labels <- c(input[[1]]$scales$scales[[1]]$labels, label)
  output[[1]] <- output[[1]] +
    ggplot2::geom_point(
      data = data.frame(time = 0, surv = 0, fstrata = "dummy"),
      ggplot2::aes(time, surv, colour = fstrata, fill = fstrata),
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::theme(
      legend.key = ggplot2::element_rect(fill = NA, colour = NA)
    )
  return(output)  
}

ggKM.quantile <- function(KM_time, KM_surv, KM_lower, KM_upper, KM_group, 
                          p = 0.5, digits.fixed = 1) {
  groups <- sort(unique(KM_group))
  output <- matrix("", length(groups), 3)
  colnames(output) <- c("surv", "lower", "upper")
  rownames(output) <- groups
  format <- paste0("%.", digits.fixed, "f")
  for (g in seq_along(groups)) {
    current_i <- KM_group == groups[g]
    current_t  <- KM_time[current_i]
    current_S  <- KM_surv[current_i]
    current_LL <- KM_lower[current_i]
    current_UL <- KM_upper[current_i]
    targeti_S  <- which(current_S  <= p)[1]
    targeti_LL <- which(current_UL <= p)[1]
    targeti_UL <- which(current_LL <= p)[1]
    if (!is.na(targeti_S))  output[g, 1] <- sprintf(format, current_t[targeti_S]) 
    if (!is.na(targeti_UL)) output[g, 2] <- sprintf(format, current_t[targeti_UL]) 
    if (!is.na(targeti_LL)) output[g, 3] <- sprintf(format, current_t[targeti_LL]) 
  }
  output[output == ""] <- "NR"
  return(output)
}

#' Add survival summary information to a ggKM legend
#'
#' Modifies the legend of a `ggKM` plot by appending the median survival time
#' and its 95% confidence interval¹ to each group label. If exactly two groups
#' are present, the legend is further extended with the hazard ratio, 95% 
#' confidence interval, and *p*-value from a Cox regression model² comparing the 
#' two groups.
#' @param input A `ggKM` plot object.
#' @return A deep copy of the input `ggKM` object with the added legend label.
#' @references
#' 1. Brookmeyer, R. and Crowley, J., 1982. A confidence interval for the median
#' survival time. *Biometrics*, pp. 29–41.
#' 2. Cox, D.R., 1972. Regression Models and Life-Tables. *Journal of the Royal 
#' Statistical Society*, Series B. 34(2), pp. 187–220.
#' @export
ggKM.surv.extra <- function(input) {
  KM_time <- attr(input, "KM_time")
  KM_surv <- attr(input, "KM_surv")
  KM_lower <- attr(input, "KM_lower")
  KM_upper <- attr(input, "KM_upper")
  KM_group <- attr(input, "KM_group")
  n_group <- length(unique(KM_group))
  output <- unserialize(serialize(input, NULL))
  labels <- input[[1]]$scales$scales[[1]]$labels
  medians <- ggKM.quantile(KM_time, KM_surv, KM_lower, KM_upper, KM_group, p = 0.5)
  for (i in 1:n_group) {
    labels[i] <- paste0(labels[i], ": median ", medians[i, 1],  
                        " (95%CI ", medians[i, 2], "\u2212", medians[i, 3], ")")
  }
  output[[1]]$scales$scales[[1]]$labels <- labels
  if (n_group == 2) {
    time <- attr(input, "time")
    status <- attr(input, "status")
    group <- attr(input, "group")
    Cox <- summary(survival::coxph(survival::Surv(time, status) ~ as.factor(group)))$coefficients
    HR <- sprintf("%.2f", Cox[1, 2])
    LL <- sprintf("%.2f", exp(Cox[1, 1] - stats::qnorm(0.975) * Cox[1, 3]))
    UL <- sprintf("%.2f", exp(Cox[1, 1] + stats::qnorm(0.975) * Cox[1, 3]))
    if (Cox[1, 5] < 0.001) {
      p <- "\U0001D45D < 0.001"
    } else if (Cox[1, 5] > 0.999) {
      p <- "\U0001D45D > 0.999"
    } else {
      p <- paste0("\U0001D45D = ", sprintf("%.2g", Cox[1, 5]))
    }
    Coxlabel <- paste0("\nHR ", HR, " (95%CI ", LL, "\u2212", UL, ")\n", p)
    output <- ggKM.legend.extra(output, Coxlabel)
  }
  return(output)
}

#' Kaplan–Meier plot: supplied survival curve versus individual-level data
#'
#' Arm 0 is a survival curve supplied as vectors (for example, pooled across
#' multiply imputed datasets); arm 1 is individual-level survival data.
#' @param KMtime0 Numeric vector of time points for arm 0, sorted in ascending order.
#' @param KMlower0 Numeric vector of lower confidence limits for arm 0.
#' @param KMsurv0 Numeric vector of survival estimates for arm 0.
#' @param KMupper0 Numeric vector of upper confidence limits for arm 0.
#' @param KMatrisk0 Numeric vector of numbers at risk for arm 0 (may be non-integer).
#' @param survtime1 Numeric vector of follow-up times for arm 1.
#' @param survstatus1 Integer vector event indicator for arm 1 (`1` = event, `0` = censored).
#' @param CI Logical; draw confidence intervals. Default = `TRUE`. Arm 0 uses the
#' supplied limits and arm 1 uses Greenwood's variance with the complementary
#' log–log transformation.
#' @param digits.fixed Decimal places for arm 0's number at risk. Default = `1`.
#' @details Other arguments are as in `ggKM`. No censor marks are drawn for arm 0.
#' @examples
#' data <- survival::lung
#' KMtime0 <- c(0.164, 0.361, 0.394, 0.427, 0.493, 0.854, 0.986, 1.019,
#'              1.741, 1.774, 1.938, 1.971, 2.004, 2.037, 2.136, 2.333,
#'              2.596, 2.661, 2.891, 3.023, 3.056, 3.121, 3.45, 3.515,
#'              3.614, 3.811, 3.877, 4.008, 4.304, 4.337, 4.435, 4.665,
#'              4.731, 4.764, 4.83, 5.027, 5.125, 5.355, 5.454, 5.487,
#'              5.585, 5.684, 5.717, 5.75, 5.782, 5.815, 5.881, 5.914,
#'              5.947, 5.98, 6.012, 6.078, 6.111, 6.177, 6.21, 6.275,
#'              6.308, 6.374, 6.44, 6.472, 6.538, 6.604, 6.637, 6.67,
#'              6.801, 6.834, 6.9, 6.932, 6.965, 7.162, 7.261, 7.294,
#'              7.327, 7.359, 7.392, 7.425, 7.524, 7.557, 7.721, 7.787,
#'              7.852, 7.885, 7.984, 8.049, 8.082, 8.279, 8.509, 8.739,
#'              8.772, 8.805, 8.838, 8.871, 8.937, 9.068, 9.167, 9.298,
#'              9.331, 9.364, 9.396, 9.462, 9.561, 9.594, 9.626, 9.725,
#'              9.856, 9.889, 9.955, 10.021, 10.054, 10.185, 10.349, 10.514,
#'              10.809, 10.908, 11.072, 11.171, 11.335, 11.433, 11.499, 11.532,
#'              11.598, 11.696, 11.861, 11.926, 11.959, 12.189, 12.353, 12.551,
#'              12.616, 12.715, 12.813, 12.945, 13.273, 13.569, 13.996, 14.062,
#'              14.095, 14.226, 14.522, 14.588, 14.785, 14.949, 15.015, 15.048,
#'              15.113, 15.54, 15.672, 16.789, 17.052, 17.085, 17.216, 17.38,
#'              17.512, 17.84, 18.07, 18.103, 18.333, 18.366, 18.629, 18.859,
#'              19.154, 19.319, 20.14, 20.501, 21.06, 21.126, 21.487, 21.52,
#'              22.571, 22.637, 23.163, 23.228, 23.918, 24.017, 24.148, 24.313,
#'              25.134, 25.988, 26.481, 26.744, 26.974, 27.598, 29.011, 31.705,
#'              33.183, 33.578)
#' KMlower0 <- c(1, 0.983, 0.958, 0.951, 0.931, 0.931, 0.918, 0.912,
#'               0.912, 0.9, 0.894, 0.888, 0.876, 0.871, 0.854, 0.849,
#'               0.843, 0.843, 0.832, 0.817, 0.812, 0.812, 0.795, 0.795,
#'               0.784, 0.779, 0.768, 0.768, 0.758, 0.747, 0.747, 0.743,
#'               0.737, 0.722, 0.717, 0.711, 0.711, 0.701, 0.677, 0.672,
#'               0.672, 0.667, 0.667, 0.662, 0.662, 0.657, 0.652, 0.637,
#'               0.627, 0.622, 0.622, 0.617, 0.617, 0.612, 0.607, 0.607,
#'               0.607, 0.602, 0.602, 0.602, 0.598, 0.582, 0.577, 0.577,
#'               0.572, 0.567, 0.561, 0.561, 0.561, 0.556, 0.551, 0.546,
#'               0.54, 0.54, 0.54, 0.54, 0.529, 0.524, 0.524, 0.524,
#'               0.524, 0.513, 0.513, 0.513, 0.507, 0.502, 0.502, 0.502,
#'               0.502, 0.496, 0.485, 0.479, 0.479, 0.479, 0.479, 0.472,
#'               0.467, 0.455, 0.455, 0.45, 0.438, 0.438, 0.438, 0.431,
#'               0.431, 0.431, 0.425, 0.413, 0.407, 0.407, 0.394, 0.388,
#'               0.388, 0.382, 0.382, 0.369, 0.363, 0.363, 0.357, 0.351,
#'               0.332, 0.332, 0.326, 0.326, 0.314, 0.308, 0.295, 0.295,
#'               0.295, 0.288, 0.288, 0.276, 0.276, 0.276, 0.276, 0.263,
#'               0.257, 0.257, 0.243, 0.237, 0.23, 0.224, 0.216, 0.216,
#'               0.216, 0.209, 0.195, 0.195, 0.187, 0.18, 0.166, 0.166,
#'               0.159, 0.159, 0.159, 0.152, 0.152, 0.144, 0.136, 0.128,
#'               0.128, 0.119, 0.119, 0.111, 0.095, 0.088, 0.088, 0.073,
#'               0.073, 0.066, 0.051, 0.051, 0.044, 0.031, 0.031, 0.024,
#'               0.017, 0.017, 0.01, 0.002, 0.002, 0.002, 0, 0,
#'               0, 0)
#' KMsurv0 <- c(1, 0.996, 0.982, 0.978, 0.964, 0.964, 0.956, 0.951,
#'              0.951, 0.942, 0.938, 0.933, 0.925, 0.92, 0.907, 0.903,
#'              0.899, 0.899, 0.89, 0.877, 0.873, 0.873, 0.859, 0.859,
#'              0.85, 0.846, 0.837, 0.837, 0.828, 0.819, 0.819, 0.815,
#'              0.81, 0.797, 0.792, 0.788, 0.788, 0.779, 0.757, 0.752,
#'              0.752, 0.748, 0.748, 0.743, 0.743, 0.739, 0.734, 0.72,
#'              0.711, 0.707, 0.707, 0.702, 0.702, 0.697, 0.693, 0.693,
#'              0.693, 0.688, 0.688, 0.688, 0.684, 0.67, 0.665, 0.665,
#'              0.66, 0.655, 0.65, 0.65, 0.65, 0.645, 0.64, 0.635,
#'              0.63, 0.63, 0.63, 0.63, 0.62, 0.615, 0.615, 0.615,
#'              0.615, 0.605, 0.605, 0.605, 0.6, 0.595, 0.595, 0.595,
#'              0.595, 0.589, 0.579, 0.573, 0.573, 0.573, 0.573, 0.568,
#'              0.562, 0.551, 0.551, 0.546, 0.534, 0.534, 0.534, 0.528,
#'              0.528, 0.528, 0.522, 0.511, 0.505, 0.505, 0.493, 0.486,
#'              0.486, 0.48, 0.48, 0.468, 0.462, 0.462, 0.456, 0.45,
#'              0.431, 0.431, 0.425, 0.425, 0.413, 0.407, 0.394, 0.394,
#'              0.394, 0.388, 0.388, 0.375, 0.375, 0.375, 0.375, 0.361,
#'              0.355, 0.355, 0.341, 0.335, 0.328, 0.321, 0.314, 0.314,
#'              0.314, 0.306, 0.292, 0.292, 0.284, 0.277, 0.262, 0.262,
#'              0.254, 0.254, 0.254, 0.246, 0.246, 0.238, 0.23, 0.221,
#'              0.221, 0.213, 0.213, 0.204, 0.187, 0.178, 0.178, 0.16,
#'              0.16, 0.151, 0.133, 0.133, 0.125, 0.106, 0.106, 0.097,
#'              0.087, 0.087, 0.078, 0.067, 0.067, 0.067, 0.05, 0.05,
#'              0.05, 0)
#' KMupper0 <- c(1, 1, 1, 1, 0.998, 0.998, 0.993, 0.99,
#'               0.99, 0.985, 0.982, 0.979, 0.973, 0.97, 0.961, 0.958,
#'               0.954, 0.954, 0.948, 0.937, 0.934, 0.934, 0.923, 0.923,
#'               0.917, 0.913, 0.906, 0.906, 0.898, 0.891, 0.891, 0.887,
#'               0.883, 0.872, 0.868, 0.864, 0.864, 0.857, 0.837, 0.833,
#'               0.833, 0.829, 0.829, 0.824, 0.824, 0.82, 0.816, 0.803,
#'               0.795, 0.791, 0.791, 0.787, 0.787, 0.783, 0.779, 0.779,
#'               0.779, 0.774, 0.774, 0.774, 0.77, 0.757, 0.753, 0.753,
#'               0.748, 0.744, 0.739, 0.739, 0.739, 0.735, 0.73, 0.725,
#'               0.72, 0.72, 0.72, 0.72, 0.711, 0.706, 0.706, 0.706,
#'               0.706, 0.697, 0.697, 0.697, 0.692, 0.687, 0.687, 0.687,
#'               0.687, 0.683, 0.673, 0.668, 0.668, 0.668, 0.668, 0.663,
#'               0.658, 0.647, 0.647, 0.642, 0.631, 0.631, 0.631, 0.625,
#'               0.625, 0.625, 0.62, 0.608, 0.603, 0.603, 0.591, 0.585,
#'               0.585, 0.579, 0.579, 0.567, 0.561, 0.561, 0.555, 0.549,
#'               0.531, 0.531, 0.524, 0.524, 0.512, 0.506, 0.494, 0.494,
#'               0.494, 0.488, 0.488, 0.474, 0.474, 0.474, 0.474, 0.46,
#'               0.453, 0.453, 0.44, 0.433, 0.426, 0.419, 0.411, 0.411,
#'               0.411, 0.403, 0.389, 0.389, 0.382, 0.374, 0.358, 0.358,
#'               0.349, 0.349, 0.349, 0.341, 0.341, 0.332, 0.323, 0.315,
#'               0.315, 0.306, 0.306, 0.298, 0.278, 0.268, 0.268, 0.247,
#'               0.247, 0.236, 0.216, 0.216, 0.205, 0.182, 0.182, 0.171,
#'               0.158, 0.158, 0.145, 0.131, 0.131, 0.131, 0.113, 0.113,
#'               0.113, 0)
#' KMatrisk0 <- c(228, 227, 223.9, 222.9, 219.9, 219.9, 217.9, 216.8,
#'                216.8, 214.8, 213.8, 212.8, 210.8, 209.8, 206.9, 205.9,
#'                204.9, 204.9, 202.9, 199, 198, 198, 193.9, 193.9,
#'                191.9, 190.9, 188.9, 188.9, 186.9, 184.8, 184.8, 183.9,
#'                182.9, 179.8, 178.8, 177.8, 177.8, 175.8, 170.8, 169.7,
#'                169.7, 167.7, 166.7, 164.7, 164.7, 163.7, 161.7, 158.6,
#'                156.6, 155.6, 155.6, 154.6, 153.6, 151.6, 150.7, 150.7,
#'                149.7, 147.8, 146.9, 146.9, 144.8, 141.8, 139.8, 138.8,
#'                137.8, 136.8, 135.8, 135.8, 134.8, 133.9, 131.8, 129.7,
#'                128.7, 128.7, 127.7, 125.6, 123.6, 122.6, 121.6, 120.7,
#'                120.7, 118.7, 116.7, 116.7, 115.7, 114.7, 113.6, 112.6,
#'                111.6, 110.6, 107.6, 106.6, 105.6, 104.6, 103.6, 102.6,
#'                100.6, 98.7, 98.7, 97.7, 95.7, 93.6, 93.6, 92.5,
#'                91.5, 90.5, 88.5, 85.5, 84.5, 84.5, 82.5, 80.6,
#'                80.6, 78.5, 78.5, 76.5, 75.5, 75.5, 74.6, 73.6,
#'                70.5, 70.5, 68.5, 68.5, 66.6, 64.5, 62.6, 60.5,
#'                60.5, 58.4, 58.4, 56.5, 56.5, 55.5, 54.6, 52.6,
#'                51.6, 51.6, 49.7, 47.7, 46.8, 45.8, 44.7, 43.7,
#'                43.7, 42.7, 40.7, 38.8, 37.8, 36.8, 34.8, 34.8,
#'                32.8, 32.8, 31.8, 30.8, 29.8, 27.9, 26.9, 25.9,
#'                25.9, 24, 24, 23, 21, 20, 20, 18,
#'                18, 17, 15, 15, 14, 12, 12, 10,
#'                9, 9, 8, 6, 5, 5, 3, 2,
#'                2, 0)
#' fig1 <- ggKM.Xvs1(KMtime0, KMlower0, KMsurv0, KMupper0, KMatrisk0,
#'                   survtime1 = data$time * 12 / 365.2425, 
#'                   survstatus1 = data$status - 1, breaks.t = seq(0, 30, 6), 
#'                   legend.labels = c("Pooled bootstrap", "Observed"),
#'                   title.s = "Overall survival", title.t = "Time (months)")
#' print(fig1)
#' @export
ggKM.Xvs1 <- function(KMtime0, KMlower0, KMsurv0, KMupper0, KMatrisk0,
                      survtime1, survstatus1,
                      breaks.s = seq(0, 1, 0.25), breaks.t = NULL,
                      CI = TRUE, CI.alpha = 0.2,
                      colors = c("#C40233", "#0087BD"),
                      digits.fixed = 1,
                      grid.color = grDevices::rgb(0.95, 0.95, 0.95),
                      grid.s = seq(0, 1, 0.25), grid.t = NULL, grid.width = 0.5,
                      legend.direction = "vertical", legend.just = c("right", "top"),
                      legend.label.pos = "left",
                      legend.labels = NULL, legend.ncol = NULL, legend.nrow = NULL,
                      legend.pos = c(1, 1), legend.text.align = 1,
                      margin.KM.b = 0, margin.KM.l = 0, margin.KM.r = 0, margin.KM.t = 4,
                      margin.risk = 16, line.width = 0.5, line.height = 0.025, max.t = NULL,
                      risk.table = 2, risk.table.prop = 0.1, t.labels = NULL,
                      textsize.axis = 12, textsize.legend = 12, textsize.risk = 12,
                      title.s = "Survival", title.t = "Time") {
  colors <- c(colors, "#00000000") # need an extra dummy color
  # ---- arm 1: ordinary KM with cloglog CI ----
  keep <- !is.na(survtime1) & !is.na(survstatus1) & survtime1 > 0
  survtime1 <- survtime1[keep]
  survstatus1 <- survstatus1[keep]
  input1 <- data.frame("time" = survtime1, "status" = survstatus1, "strata" = 1)
  fit1 <- survival::survfit(survival::Surv(time, status) ~ 1, data = input1, conf.type = "log-log")
  s1 <- summary(fit1)
  summary1 <- data.frame("time" = s1$time, "surv" = s1$surv,
                         "lower" = s1$lower, "upper" = s1$upper, "strata" = 1)
  # ---- arm 0: pooled KM ----
  summary0 <- data.frame("time" = KMtime0, "surv" = KMsurv0,
                         "lower" = KMlower0, "upper" = KMupper0, "strata" = 0)
  input0 <- data.frame("time" = max(KMtime0), "status" = NA, "strata" = 0) # ggKM.step only uses its max time
  # ---- combine ----
  data_summary <- rbind(summary0, summary1)
  data_input <- rbind(input0, input1)
  if (is.null(max.t)) max.t <- max(data_input$time)
  if (is.null(breaks.t)) breaks.t <- seq(0, max.t, by = 12)
  if (is.null(legend.labels)) legend.labels <- 0:1
  if (is.null(t.labels)) t.labels <- breaks.t
  x_lim <- c(0, max.t)
  levels_order <- c(0, 1)
  data_step <- ggKaplanMeier:::ggKM.step(data_summary, data_input, max.t)
  data_step$fstrata <- factor(data_step$strata, levels = levels_order)
  data_pruned <- data_step[!is.na(data_step$lower),]
  # ---- KM plot ----
  g_KM <- ggplot2::ggplot(data_step, ggplot2::aes(time, surv, color = fstrata, fill = fstrata)) +
    ggplot2::coord_cartesian(ylim = c(0, 1), xlim = x_lim, clip = "off") +
    ggplot2::geom_hline(yintercept = grid.s, color = grid.color, linewidth = grid.width) +
    ggplot2::geom_vline(xintercept = grid.t, color = grid.color, linewidth = grid.width) +
    ggplot2::guides(color = ggplot2::guide_legend(label.position = legend.label.pos, ncol = legend.ncol, nrow = legend.nrow),
                    fill  = ggplot2::guide_legend(label.position = legend.label.pos, ncol = legend.ncol, nrow = legend.nrow)) +
    ggplot2::labs(x = title.t, y = title.s, color = NULL) +
    ggplot2::scale_color_manual(name = "", values = colors, labels = legend.labels,
                                aesthetics = c("color", "fill")) +
    ggplot2::scale_x_continuous(breaks = breaks.t, labels = t.labels, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = breaks.s, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = textsize.axis),
                   axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 1)),
                   axis.title.x = ggplot2::element_text(size = textsize.axis),
                   legend.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.box.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.direction = legend.direction,
                   legend.justification = legend.just,
                   legend.position = legend.pos,
                   legend.text = ggplot2::element_text(size = textsize.legend),
                   legend.text.align = legend.text.align,
                   legend.title = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = margin.KM.t, r = margin.KM.r,
                                                 b = margin.KM.b, l = margin.KM.l))
  if (CI) {
    g_KM <- g_KM +
      ggplot2::geom_ribbon(data = data_pruned,
                           ggplot2::aes(ymin = lower, ymax = upper, fill = fstrata),
                           alpha = CI.alpha, colour = NA)
  }
  # censor marks: arm 1 only
  data_censor <- ggKaplanMeier:::ggKM.censor(summary1, input1)
  if (!is.null(data_censor$time)) {
    data_censor$fstrata <- factor(data_censor$strata, levels = levels_order)
    g_KM <- g_KM + ggplot2::geom_segment(
      data = data_censor[data_censor$time <= max.t,],
      ggplot2::aes(x = time, xend = time,
                   y = surv - line.height/2, yend = surv + line.height/2,
                   color = fstrata),
      linewidth = line.width
    )
  }
  g_KM <- g_KM + ggplot2::geom_line(linewidth = line.width)
  # ---- risk table ----
  if (risk.table >= 1) {
    # arm 0: mean number at risk at each break = value at the first time point >= break
    j <- findInterval(breaks.t, KMtime0, left.open = TRUE) + 1
    risk0 <- rep(0, length(breaks.t))
    inside <- j <= length(KMtime0)
    risk0[inside] <- KMatrisk0[j[inside]]
    # arm 1: ordinary number at risk
    risk1 <- summary(fit1, times = breaks.t, extend = TRUE)$n.risk
    data_risk <- data.frame(
      time = rep(breaks.t, 2),
      strata = rep(c("0", "1"), each = length(breaks.t)),
      n_risk = c(formatC(risk0, format = "f", digits = digits.fixed),
                 formatC(risk1, format = "f", digits = 0))
    )
    if (risk.table == 1) {
      data_risk <- rbind(data.frame(time = 0, strata = "At risk:", n_risk = ""), data_risk)
      data_risk$strata <- factor(data_risk$strata, levels = c("At risk:", "0", "1"),
                                 labels = c("At risk:", legend.labels))
      g_risk <- ggplot2::ggplot(data_risk[data_risk$time <= max.t,],
                                ggplot2::aes(time, strata, label = n_risk)) +
        ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
        ggplot2::geom_text(size = textsize.risk * 127/360) +
        ggplot2::scale_x_continuous(breaks = breaks.t, expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          axis.line = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = margin.risk),
                                              size = textsize.risk),
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank()
        )
    } else if (risk.table == 2) {
      data_risk$strata <- factor(data_risk$strata, levels = levels_order)
      colors_risk <- colors[1:2]
      g_risk <- ggplot2::ggplot(
        data_risk[data_risk$time <= max.t,],
        ggplot2::aes(time, strata, label = n_risk)) +
        ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
        ggplot2::geom_text(ggplot2::aes(color = strata), size = textsize.risk * 127/360) +
        ggplot2::scale_color_manual(values = colors_risk, guide = "none") +
        ggplot2::scale_x_continuous(breaks = breaks.t, expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev(levels(data_risk$strata))) +
        ggplot2::theme_classic() + ggplot2::ggtitle(" ") +
        ggplot2::theme(
          axis.line   = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks  = ggplot2::element_blank(),
          axis.title  = ggplot2::element_blank(),
          panel.grid  = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0, size = textsize.axis,
                                             margin = ggplot2::margin(b = 0)),
          plot.title.position = "plot"
        ) + ggplot2::labs(y = NULL) +
        ggplot2::annotation_custom(
          grid::textGrob(
            "At risk", rot = 90, x = grid::unit(-2.75, "lines"),
            gp = grid::gpar(fontsize = textsize.axis)
          )
        )
    } else stop("[ggKM.Xvs1] invalid risk.table")
    if (risk.table == 2) {
      y_lab <- "" # dummy title for spacing
    } else {
      y_lab <- NULL
    }
    g_KM <- g_KM + ggplot2::labs(y = y_lab) +
      ggplot2::annotation_custom(
        grid::textGrob(
          title.s, rot = 90, x = grid::unit(-2.75, "lines"),
          gp = grid::gpar(fontsize = textsize.axis)
        )
      )
    g <- g_KM / g_risk +
      patchwork::plot_layout(heights = c(((1 - risk.table.prop) / risk.table.prop), 1))
  } else {
    g <- g_KM
  }
  attr(g, "KM_time") <- data_pruned$time
  attr(g, "KM_surv") <- data_pruned$surv
  attr(g, "KM_lower") <- data_pruned$lower
  attr(g, "KM_upper") <- data_pruned$upper
  attr(g, "KM_group") <- data_pruned$strata
  return(g)
}