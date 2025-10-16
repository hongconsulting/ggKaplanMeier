#' @import patchwork
utils::globalVariables(c("fstrata", "lower", "n_risk", "surv", "strata", "upper"))

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

ggKM.step <- function(data_summary, data_input) {
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
    endrow <- newdata[nrow(newdata),]
    endrow$time <- max(subdata_input$time)
    output[[as.character(f)]] <- rbind(newdata, endrow)
  }
  return(do.call(rbind, output))
}

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

ggKM.CI <- function(data_input, method) {
  if (!requireNamespace("km.ci", quietly = TRUE)) stop("[ggKM.CI] requires package 'km.ci'")
  factors <- unique(data_input$strata)
  output <- list()
  for (f in factors) {
    subdata_input <- data_input[data_input$strata == f,]
    fit0 <- survival::survfit(survival::Surv(time, status) ~ 1, data = subdata_input)
    fit <- km.ci::km.ci(fit0, conf.level = 0.95, tl = NA, tu = NA, method = method)
    summary <- summary(fit)
    output[[as.character(f)]] <- data.frame(
      "time" = summary$time,
      "surv" = summary$surv,
      "lower" = summary$lower,
      "upper" = summary$upper,
      "strata" = f)
  }
  return(do.call(rbind, output))
}

#' Kaplan–Meier plot
#'
#' Generates a Kaplan–Meier plot with optional confidence intervals and risk table.
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator (`1` = event, `0` = censored).
#' @param group Integer vector grouping variable.
#' @param breaks.s Y-axis (survival) tick marks. Default = `seq(0, 1, 0.25)`.
#' @param breaks.t X-axis (time) tick marks. Default = `seq(0, max(time), by = 12)`.
#' @param CI Integer indicating the CI type (options `2`, `3` and `5` are via 
#' the `km.ci` package, option `4` is via the `bpcp` package):
#'   \itemize{
#'     \item \code{0}: none
#'     \item \code{1}: pointwise CI using Greenwood's variance¹ with complementary log–log transformation² (default)
#'     \item \code{2}: pointwise CI using Rothman's binomial method³ 
#'     \item \code{3}: pointwise CI using the Thomas–Grunkemeier likelihood-ratio method⁴
#'     \item \code{4}: pointwise CI using the Fay–Brittain beta product confidence procedure⁵
#'     \item \code{5}: simultaneous confidence bands using Nair's equal precision method⁶ with log transformation 
#'   }
#' @param CI.alpha Alpha transparency of the confidence intervals. Default = `0.2`.
#' @param colors Vector of colors of the survival curves. Default = `ggsci::pal_nejm()(8)`.
#' @param grid.color Gridline color. Default = `grDevices::rgb(0.95, 0.95, 0.95)`.
#' @param grid.s Horizontal gridline positions. Default = `seq(0, 1, 0.25)`.
#' @param grid.t Vertical gridline positions. Default = `NULL`.
#' @param grid.width Gridline thickness. Default = `0.5`.
#' @param legend.direction Legend orientation; either `"vertical"` (default) or `"horizontal"`.
#' @param legend.justification Alignment anchor for the legend relative to its
#'   position. Can be a single keyword (e.g., `"center"`), a keyword pair specifying
#'   horizontal and vertical justification (e.g., `c("left", "top")`, `c("right", "bottom")`) 
#'   respectively, or a numeric vector of length 2 giving relative coordinates 
#'   within the plot area. Default = `"center"`.
#' @param legend.label.position Legend label position (`"left"` or `"right"`) relative to the legend symbol. Default = `"left"`.
#' @param legend.labels Character vector of group labels.
#' @param legend.ncol Integer specifying the number of columns in the legend. Default = `NULL`.
#' @param legend.nrow Integer specifying the number of rows in the legend. Default = `NULL`.
#' @param legend.position Position of the legend. Can be a keyword such as 
#'   `"none"`, `"left"`, `"right"`, `"bottom"`, or `"top"`, or a numeric vector 
#'   of length 2 giving relative coordinates within the plot area. Default = `c(0.9, 0.9)`.
#' @param legend.text.align Legend text alignment: `0` = left, `0.5` = center, `1` = right. Default = `1`.
#' @param line.width Line width survival curves and censor marks. Default = `0.5`.
#' @param line.height Height of censor marks. Default = `0.025`.
#' @param risk.table Logical; if `TRUE`, show risk table. Default = `TRUE`.
#' @param risk.table.proportion Relative height of the risk table. Default = `0.2`.
#' @param textsize.axis Axis text size. Default = `12`.
#' @param textsize.legend Legend text size. Default = `12`.
#' @param textsize.risk Risk table text size. Default = `12`.
#' @param title.s Y-axis (survival) title. Default = `"Survival"`.
#' @param title.t X-axis (time) title. Default = `"Time"`.
#' @return A `ggplot2` object (or `patchwork` composite if `risk.table = TRUE`).
#' @references
#' 1. Greenwood, M., 1926. A report on the natural duration of cancer. In: 
#' \emph{Reports on Public Health and Medical Subjects}, 33, pp. 1–26. London: 
#' Her Majesty’s Stationery Office, Ministry of Health.
#' 2. Klein, J.P., Logan, B., Harhoff, M. and Andersen, P.K., 2007. Analyzing 
#' survival curves at a fixed point in time. \emph{Statistics in Medicine},
#' 26(24), pp. 4505–4519.
#' 3. Rothman, K.J., 1978. Estimation of confidence limits for the cumulative 
#' probability of survival in life table analysis. \emph{Journal of Chronic 
#' Diseases}, 31(8), pp. 557–560.
#' 4. Thomas, D.R. and Grunkemeier, G.L., 1975. Confidence interval estimation 
#' of survival probabilities for censored data. \emph{Journal of the American 
#' Statistical Association}, 70(352), pp. 865–871.
#' 5. Fay, M.P. and Brittain, E.H., 2016. Finite sample pointwise confidence 
#' intervals for a survival distribution with right‐censored data.
#' \emph{Statistics in Medicine}, 35(16), pp. 2726–2740.
#' 6. Nair, V.N., 1984. Conﬁdence bands for survival functions with censored
#' data: a comparative study. \emph{Technometrics}, 26, pp. 265–275.
#' @examples
#' data <- survival::lung
#' g <- ggKM(data$time * 12 / 365.2425, data$status - 1, data$sex, 
#'           breaks.t = seq(0, 30, 6), legend.labels = c("Male", "Female"),
#'           title.s = "Overall survival", title.t = "Time (months)")
#' print(g)
#' @export
ggKM <- function(time, status, group,
                 breaks.s = seq(0, 1, 0.25), breaks.t = NULL,
                 CI = 1, CI.alpha = 0.2,
                 colors = ggsci::pal_nejm()(8),
                 grid.color = grDevices::rgb(0.95, 0.95, 0.95),
                 grid.s = seq(0, 1, 0.25), grid.t = NULL, grid.width = 0.5,
                 legend.direction = "vertical", legend.justification = "center",
                 legend.label.position = "left", 
                 legend.labels = NULL, legend.ncol = NULL, legend.nrow = NULL,
                 legend.position = c(0.9, 0.9), legend.text.align = 1,
                 line.width = 0.5, line.height = 0.025,
                 risk.table = TRUE, risk.table.proportion = 0.2,
                 textsize.axis = 12, textsize.legend = 12, textsize.risk = 12,
                 title.s = "Survival", title.t = "Time") {
  group <- as.numeric(group)
  n_group <- length(unique(group))
  if (is.null(breaks.t)) breaks.t <- seq(0, max(time), by = 12)
  if (is.null(legend.labels)) legend.labels <- 0:(n_group - 1)
  data_input <- data.frame("time" = time, "status" = status, "strata" = group)
  data_censor <- data_input[data_input$status == 0,]
  fit <- survival::survfit(survival::Surv(time, status) ~ strata, data = data_input, conf.type = "log-log")
  s <- summary(fit)
  if (is.null(s$strata)) {
    s$strata <- rep(0, length(s$time))
  } else {
    s$strata <- as.numeric(sub("strata=", "", s$strata))
  }
  if (CI == 0 | CI == 1) {
    data_summary <- data.frame(s[c("time", "surv", "lower", "upper", "strata")])
  } else if (CI == 2) {
    data_summary <- ggKM.CI(data_input, "rothman")
  } else if (CI == 3) {
    data_summary <- ggKM.CI(data_input, "grunkemeier")
  } else if (CI == 4) {
    data_summary <- ggKM.BPCP(data_input)
  } else if (CI == 5) {
    data_summary <- ggKM.CI(data_input, "logep")
  } else {
    stop("[ggKM] invalid CI")
  }
  data_step <- ggKM.step(data_summary, data_input)
  levels_order <- sort(unique(group)) # enforce factor level order consistency
  data_step$fstrata <- factor(data_step$strata, levels = levels_order)
  x_lim <- c(0, max(time))
  g_KM <- ggplot2::ggplot(data_step, ggplot2::aes(time, surv, color = fstrata, fill = fstrata)) +
    ggplot2::coord_cartesian(ylim = c(0, 1), xlim = x_lim, clip = "off") +
    ggplot2::geom_hline(yintercept = grid.s, color = grid.color, linewidth = grid.width) +
    ggplot2::geom_vline(xintercept = grid.t, color = grid.color, linewidth = grid.width) +
    ggplot2::guides(color = ggplot2::guide_legend(label.position = legend.label.position, ncol = legend.ncol, nrow = legend.nrow),
                    fill  = ggplot2::guide_legend(label.position = legend.label.position, ncol = legend.ncol, nrow = legend.nrow)) +
    ggplot2::labs(x = title.t, y = title.s, color = NULL) +
    ggplot2::scale_color_manual(
      name = "",
      values = colors,
      labels = legend.labels,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::scale_x_continuous(breaks = breaks.t, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = breaks.s, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = textsize.axis),
                   axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 1)),
                   axis.title.x = ggplot2::element_text(size = textsize.axis),
                   legend.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.box.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.direction = legend.direction, 
                   legend.justification = legend.justification,
                   legend.position = legend.position,
                   legend.text = ggplot2::element_text(size = textsize.legend),
                   legend.text.align = legend.text.align)
  if (CI > 0) {
    g_KM <- g_KM +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = fstrata),
                           alpha = CI.alpha, colour = NA)
  }
  data_censor <- ggKM.censor(data_summary, data_input)
  data_censor$fstrata <- factor(data_censor$strata, levels = levels_order)
  g_KM <- g_KM + ggplot2::geom_segment(
    data = data_censor,
    ggplot2::aes(x = time, xend = time,
                 y = surv - line.height/2, yend = surv + line.height/2,
                 color = fstrata),
    linewidth = line.width
  )
  g_KM <- g_KM + ggplot2::geom_line(linewidth = line.width)
  if (risk.table) { 
    s <- summary(fit, times = breaks.t)
    data_risk <- data.frame(
      time = s$time,
      strata = NA,
      n_risk = s$n.risk
    )
    if (length(unique(s$strata)) > 1) {
      data_risk$strata <- as.factor(as.numeric(sub("strata=", "", s$strata)))      
    }
    data_risk <- rbind(data.frame(time = 0, strata = "At risk:", n_risk = ""), data_risk)
    data_risk$strata <- as.factor(data_risk$strata)
    data_risk$strata <- stats::relevel(data_risk$strata, "At risk:")
    levels(data_risk$strata) <- c("At risk:", legend.labels)
    g_risk <- ggplot2::ggplot(data_risk, ggplot2::aes(time, strata, label = n_risk)) +
      ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
      ggplot2::geom_text(size = textsize.risk * 127/360) +
      ggplot2::scale_x_continuous(breaks = breaks.t, expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0), limits = rev) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 16), size = textsize.risk),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          hjust = 0,
          size = textsize.axis,
          margin = ggplot2::margin(b = 32)
        ),
        plot.title.position = "plot"
      )
    g_KM <- g_KM + ggplot2::labs(y = NULL) +
      ggplot2::annotation_custom(
        grid::textGrob(
          title.s, rot = 90, x = grid::unit(-2.75, "lines"),
          gp = grid::gpar(fontsize = textsize.axis)
        )
      )
    g <- g_KM / g_risk +
      patchwork::plot_layout(heights = c(((1 - risk.table.proportion) / risk.table.proportion), 1))
  }
  return(g)
}
