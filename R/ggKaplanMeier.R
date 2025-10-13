#' @import patchwork
utils::globalVariables(c("surv", "strata", "lower", "upper", "n_risk"))

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
    output[[f]] <- data.frame(
      "time" = subdata_input$time,
      "surv" = subdata_summary$surv[findInterval(subdata_input$time, subdata_summary$time)],
      "strata" = f)
  }
  return(do.call(rbind, output))
}

ggKM.step <- function(summary, fit) {
  factors <- unique(summary$strata)
  output <- list()
  i <- 0
  for (f in factors) {
    subdata <- summary[summary$strata == f,]
    if (subdata$time[1] != 0) {
      startrow <- data.frame("time" = 0, "surv" = 1, "lower" = 1, "upper" = 1, "strata" = f)
      subdata <- rbind(startrow, subdata)
    }
    end <- 2 * nrow(subdata)
    newdata <- data.frame(
      "time" = rep(subdata$time, each = 2)[-1],
      "surv" = rep(subdata$surv, each = 2)[-end],
      "lower" = rep(subdata$lower, each = 2)[-end],
      "upper" = rep(subdata$upper, each = 2)[-end],
      "strata" = f)
    endrow <- newdata[nrow(newdata),]
    i <- i + 1
    j <- seq_len(fit$strata[i]) + sum(fit$strata[seq_len(i - 1)])
    endrow$time <- max(fit$time[j])
    output[[f]] <- rbind(newdata, endrow)
  }
  return(do.call(rbind, output))
}

#' Kaplan–Meier plot
#'
#' Generates a Kaplan–Meier plot with optional confidence intervals and risk table.
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator (1 = event, 0 = censored).
#' @param group Grouping variable.
#' @param breaks.s Y-axis (survival) tick marks. Default = `seq(0, 1, 0.25)`.
#' @param breaks.t X-axis (time) tick marks. Default = `seq(0, max(time), by = 6)`.
#' @param CI Logical; if `TRUE`, show confidence intervals. Default = `TRUE`.
#' @param CI.alpha Alpha transparency of confidence intervals. Default = 0.2.
#' @param colors Vector of colors. Default = `ggsci::pal_nejm()(8)`.
#' @param grid.color Gridline color. Default = `rgb(0.95, 0.95, 0.95)`.
#' @param grid.s Horizontal gridline positions. Default = `seq(0, 1, 0.25)`.
#' @param grid.t Vertical gridline positions. Default = `NULL`.
#' @param grid.width Gridline thickness. Default = 0.5.
#' @param legend.label.position Legend label position ("left" or "right") relative to the legend symbol. Default = "left".
#' @param legend.labels Character vector of group labels.
#' @param legend.position Legend position coordinates. Default = `c(0.9, 0.9)`.
#' @param legend.text.align Legend text alignment: 0 = left, 0.5 = center, 1 = right. Default = 1.
#' @param line.width Line width survival curves and censor marks. Default = 0.5.
#' @param line.height Height of censor marks. Default = 0.025.
#' @param risk.table Logical; if `TRUE`, show risk table. Default = `TRUE`.
#' @param risk.table.proportion Relative height of the risk table. Default = 0.2.
#' @param textsize.axis Axis text size. Default = 12.
#' @param textsize.legend Legend text size. Default = 12.
#' @param textsize.risk Risk table text size. Default = 12.
#' @param title.s Y-axis (survival) title. Default = "Survival".
#' @param title.t X-axis (time) title. Default = "Time".
#' @return A `ggplot2` object (or `patchwork` composite if `risk.table = TRUE`).
#' @examples
#' data <- survival::lung
#' g <- ggKM(data$time * 12 / 365.2425, data$status - 1, data$sex,
#'           legend.labels = c("Male", "Female"),
#'           title.s = "Overall survival", title.t = "Time (months)")
#' print(g)
#' @export
ggKM <- function(time, status, group,
                 breaks.s = seq(0, 1, 0.25), breaks.t = NULL,
                 CI = TRUE, CI.alpha = 0.2,
                 colors = ggsci::pal_nejm()(8),
                 grid.color = grDevices::rgb(0.95, 0.95, 0.95),
                 grid.s = seq(0, 1, 0.25), grid.t = NULL, grid.width = 0.5,
                 legend.label.position = "left", legend.labels = NULL,
                 legend.position = c(0.9, 0.9), legend.text.align = 1,
                 line.width = 0.5, line.height = 0.025,
                 risk.table = TRUE, risk.table.proportion = 0.2,
                 textsize.axis = 12, textsize.legend = 12, textsize.risk = 12,
                 title.s = "Survival", title.t = "Time") {
  n_group <- length(unique(group))
  if (is.null(breaks.t)) breaks.t <- seq(0, max(time), by = 6)
  if (is.null(legend.labels)) legend.labels <- 0:(n_group - 1)
  data_input <- data.frame("time" = time, "status" = status, "strata" = as.factor(group))
  data_censor <- data_input[data_input$status == 0,]
  fit <- survival::survfit(survival::Surv(time, status) ~ strata, data = data_input)
  s <- summary(fit)
  if (is.null(s$strata)) {
    s$strata <- rep(0, length(s$time))
  } else {
    s$strata <- as.numeric(sub("strata=", "", s$strata))
  }
  data_summary <- data.frame(s[c("time", "surv", "lower", "upper", "strata")])
  data_step <- ggKM.step(data_summary, fit)
  data_step$strata <- as.factor(data_step$strata)
  x_lim <- c(0, max(time))
  g_KM <- ggplot2::ggplot(data_step, ggplot2::aes(time, surv, color = strata)) +
    ggplot2::coord_cartesian(ylim = c(0, 1), xlim = x_lim, clip = "off") +
    ggplot2::geom_hline(yintercept = grid.s, color = grid.color, linewidth = grid.width) +
    ggplot2::geom_vline(xintercept = grid.t, color = grid.color, linewidth = grid.width) +
    ggplot2::guides(color = ggplot2::guide_legend(label.position = legend.label.position),
                    fill  = ggplot2::guide_legend(label.position = legend.label.position)) +
    ggplot2::labs(x = title.t, y = title.s, color = NULL) +
    ggplot2::scale_color_manual(name = "", values = colors, labels = legend.labels) +
    ggplot2::scale_fill_manual(name = "", values = colors, labels = legend.labels) +
    ggplot2::scale_x_continuous(breaks = breaks.t, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = breaks.s, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = textsize.axis),
                   axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 1)),
                   axis.title.x = ggplot2::element_text(size = textsize.axis),
                   legend.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.box.background = ggplot2::element_rect(fill = NA, color = NA),
                   legend.position = legend.position,
                   legend.text = ggplot2::element_text(size = textsize.legend),
                   legend.text.align = legend.text.align)
  if (CI) {
    g_KM <- g_KM +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = strata),
                           alpha = CI.alpha, colour = NA)
  }
  data_censor <- ggKM.censor(data_summary, data_input)
  data_censor$strata <- as.factor(data_censor$strata)
  g_KM <- g_KM + ggplot2::geom_segment(
    data = data_censor,
    ggplot2::aes(x = time, xend = time,
                 y = surv - line.height/2, yend = surv + line.height/2,
                 color = strata),
    linewidth = line.width
  )
  g_KM <- g_KM + ggplot2::geom_line(linewidth = line.width)
  if (risk.table) {
    s <- summary(fit, times = breaks.t)
    data_risk <- data.frame(
      time = s$time,
      strata = as.factor(as.numeric(sub("strata=", "", s$strata))),
      n_risk = s$n.risk
    )
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
