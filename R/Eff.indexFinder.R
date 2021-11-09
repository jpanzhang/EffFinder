#' Eff.indexFinder, to detect the Normalization efficiency of any normalization factor (single reference gene or possible multi-RG combination)
#' according to its serial number from Index file
#'
#' @param NF_index The serial number of any NF, you can find from the Index file
#'
#' @return NEdata_order The normalization efficiency (NE) values of normalization factors you want to detect
#' @export
#'
#' @examples Eff.indexFinder(NF_index=c(1:12))
#' @examples Eff.indexFinder(NF_index=12:20)
#'
Eff.indexFinder <- function(NF_index)
{
  library(pROC)
  library(grDevices)
  par(mfrow = c(1, 1),
      mar = c(4, 4, 1, 1),
      oma = c(0.5, 0.5, 0.5, 0.5)  )
  ####### define the internal Function ########
  Eff.ne.1 <- function(RG_order) {
    pvdata2 = data.pv
    pvall_singleRGs = t(pvdata2[RG_order, ])               #loading p-value of each single RG
    tp_real = sum(gold_standard == 1)                      #positive frequency of gold standard
    tn_real = sum(gold_standard == 2)                      #negative frequency of gold standard
    kk = length(RG_order)
    ramp <- colorRamp(c("red", "orange", "darkcyan", "blue"))
    colorpool <- rgb(ramp(seq(0, 1, length = kk)), max = 255)
    for (ii in 1:kk)
    {
      if (ii == 1) {
        rocobj <- plot.roc(
          gold_standard,
          pvall_singleRGs[, ii],
          percent = TRUE,
          col = colorpool[ii],
          grid = TRUE,
          smooth = TRUE,
          xlab = "FPR (%)",
          ylab = "TPR (%)",
          legacy.axes = TRUE,
          main = paste0("ROC space (Gold standard at p=", pv.th[1], ")")
        )
        NE = round((rocobj$auc / 100 - 0.5) / 0.5 * 100, digit <- 4)
        pv_x =  pvall_singleRGs[, ii]
        (pv_x[pv_x > pv.th] <-          2) &
          (pv_x[pv_x < pv.th] <-
             1)            #2 means negative (p>x), #1 means positive (p<x)
        pv_x_GSvsRG = apply(rbind(gold_standard, pv_x), 2, mean)
        tp = sum(pv_x_GSvsRG == 1)
        tn = sum(pv_x_GSvsRG == 2)
        fn = tp_real - tp
        fp = tn_real - tn
        TPR = round(tp / tp_real * 100 , digit <- 4)
        TNR = round(tn / tn_real * 100 , digit <- 4)
        FPR = round(fp / tn_real * 100 , digit <- 4)
        FNR = round(fn / tp_real * 100 , digit <- 4)
        Precision = round(tp / (tp + fp) * 100 , digit <- 4)
        Recall = round(tp / tp_real * 100 , digit <- 4)
        NEdata_order = data.frame(
          RG_order[ii],
          paste0("NF", RG_order[ii]),
          index[,2][RG_order[ii]],
          NE,
          TPR,
          TNR,
          FPR,
          FNR,
          Precision,
          Recall
        )
        text(
          20,65,
          labels = paste0("——", index[,2][RG_order[ii]], " (NE) = ", sprintf("%0.2f", NE), "%"),
          col = colorpool[ii],
          cex = 1.2
        )
      }
      else {
        rocobj <-          lines.roc(
          gold_standard,
          pvall_singleRGs[, ii],
          percent = TRUE,
          col = colorpool[ii],
          grid = TRUE,
          smooth = TRUE
        )
        NE = round((rocobj$auc / 100 - 0.5) / 0.5 * 100, digit <- 4)
        pv_x =  pvall_singleRGs[, ii]
        (pv_x[pv_x > pv.th] <-  2) &
          (pv_x[pv_x < pv.th] <-  1)            #2 means negative (p>x), #1 means positive (p<x)
        pv_x_GSvsRG = apply(rbind(gold_standard, pv_x), 2, mean)
        tp = sum(pv_x_GSvsRG == 1)
        tn = sum(pv_x_GSvsRG == 2)
        fn = tp_real - tp
        fp = tn_real - tn
        TPR = round(tp / tp_real * 100 , digit <- 4)
        TNR = round(tn / tn_real * 100 , digit <- 4)
        FPR = round(fp / tn_real * 100 , digit <- 4)
        FNR = round(fn / tp_real * 100 , digit <- 4)
        Precision = round(tp / (tp + fp) * 100 , digit <- 4)
        Recall = round(tp / tp_real * 100 , digit <- 4)
        NEdata_order2 = data.frame(
          RG_order[ii],
          paste0("NF", RG_order[ii]),
          index[,2][RG_order[ii]],
          NE,
          TPR,
          TNR,
          FPR,
          FNR,
          Precision,
          Recall
        )
        NEdata_order = rbind(NEdata_order, NEdata_order2)
        text(
          20,
          (70 - ii * 5),
          labels = paste0("——", index[,2][RG_order[ii]], " (NE) = ", sprintf("%0.2f", NE), "%"),
          col = colorpool[ii],
          cex = 1.2
        )
      }
    }
    colnames(NEdata_order) = c(
      "NF_index",
      "NF_index",
      "NF_name",
      "NE(%)",
      "TPR(%)",
      "TNR(%)",
      "FPR(%)",
      "FNR(%)",
      "Precision(%)",
      "Recall(%)"
    )
    NEdata_order <<- NEdata_order
  }
  #####main codes ##############
  RG_order1 = NF_index
  Eff.ne.1(RG_order1)
  RG_order2 = as.numeric(NEdata_order[order(NEdata_order$NE, decreasing = T),][, 1])
  Eff.ne.1(RG_order2)
  NEdata_order$Ranking_order = c(1:length(RG_order2))
  NEdata_order = NEdata_order[,-1]
  NEdata_order <<- NEdata_order
  print(NEdata_order)
}
