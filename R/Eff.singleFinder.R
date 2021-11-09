#' Eff.singleFinder, to detect the Normalization efficiencies of single reference genes (RGs)
#'
#' @family
#' @param rawdata The qRT-PCR data matrix. The column-names and row-names are the sample names and gene symbols, respectively.
#' The row-names preceded by reference genes (RGs) and followed by target genes (TGs)
#' @param num_RG The number of reference genes (RGs) in your raw data matrix
#' @param num_TG The number of target genes (TGs) in your raw data matrix
#' @param num_group The number of experimental groups in raw data
#' @param num_sample The number of biological samples in your raw data matrix, is equal to Num_rep multiplied by num_group
#' @param bio_rep The number of biological replicate in all experimental groups
#' @param name_RG The gene names of your all reference genes (RGs)
#' @param name_TG The gene names of your all target genes (TGs)
#' @param Dis.cq The Dis.cq is the setting gap of Cq value, between the maximum value and second maximum,
#' or between the minimum and second minimum in the same experimental group. It is generally set to 2
#' @param gsg.Ratio The ratio of the gold standard group filtered from all NFs, is generally set to 0.01
#' @param pv.th The threshold for identifying positive and negative results between any two paired comparison groups, it is generally set to 0.05
#' @param Primer.E The amplification efficiency (E) of all primers, should be a vector here.
#' If you have not tested the primers' amplification efficiencies, please let Primer.E=1
#'
#' @return
#' @export
#'
#' @examples Eff.singleFinder(rawdata, num_RG=12, num_TG=20, num_sample=54, num_group=18, bio_rep=3,
#'                   name_RG, name_TG, Dis.cq=2, gsg.Ratio=0.01, pv.th=0.05, Primer.E)
#'
Eff.singleFinder <- function(rawdata,
                             num_RG,
                             num_TG,
                             num_group,
                             num_sample,
                             bio_rep,
                             name_RG,
                             name_TG,
                             Dis.cq,
                             gsg.Ratio,
                             pv.th,
                             Primer.E)
{
  ###########################################################################
  ###########################  Step 1: Eff.dpp  #############################
  ###########################################################################
  Eff.dpp <- function(rawdata,
                      num_RG,
                      num_TG,
                      num_group,
                      bio_rep,
                      Dis.cq,
                      Primer.E)
  {
    rwd = rawdata
    message = paste(
      " =========================================================================================================",
      "\n",
      "The Eff.dpp:"
    )
    cat(message)
    ####################--Define the internal function 1--####################
    Eff.dpp.1 <- function(rwd)
    {
      if (length(Primer.E) == nrow(rawdata))
      {
        aa = (1 + Primer.E) ^ rwd
        rwd <<- log2(aa)  #Outputting global variable
        print("Eff.dpp.1 has been completed!")
      }
      else {
        stop(
          "The vector length of Primer.E is not equal to the row number of rawdata matrix, please check it!"
        )
      }
    }
    ####################--Define the internal function 2--####################
    Eff.dpp.2 <- function(rwd)
    {
      repeat {
        rwd = rwd
        mean_sample = apply(rwd, 2, mean, na.rm = T)
        mean_all = mean(mean_sample)
        cf = mean_sample - mean_all
        if (as.numeric(sprintf("%0.5f", min(cf))) == 0 &
            as.numeric(sprintf("%0.5f", min(cf))) == 0)
        {
          print("Eff.dpp.2 has been completed!")
          break
        }
        else {
          rwd = (t(t(rwd) - cf))
        }
      }
      rwd <<- rwd  #Outputting global variable
    }
    ####################--Define the internal function 3--####################
    Eff.dpp.3 <- function(rwd)
    {
      for (ii in 1:num_group)
      {
        detect.col = c(((ii - 1) * bio_rep + 1):(ii * bio_rep))
        mean.col = apply(rwd[, detect.col], 1, mean, na.rm = T)   #Calculate the mean of same gene in one group
        for (jj in 1:(num_RG + num_TG))
        {
          rwd[jj, detect.col][is.na(rwd[jj, detect.col])] = mean.col[jj]  #Replace null value with corresponding mean value
        }
        rwd <<- rwd   #Outputting global variable
      }
      print("Eff.dpp.3 has been completed!")
    }
    ####################--Define the internal function 4--####################
    Eff.dpp.4 <- function(rwd)
    {
      Eff.thva <- function(thva)
      {
        repeat {
          thva_order = order(thva, decreasing = T)
          thva_ordva = thva[thva_order]
          ll = length(thva)
          A = thva_ordva[1]
          B = thva_ordva[2]
          C = thva_ordva[ll - 1]
          D = thva_ordva[ll]
          if ((A - B >= Dis.cq &
               A > (mean(thva_ordva[c(2:ll)]) + sd(thva_ordva[c(2:ll)]) * 3)) |
              (C - D >= Dis.cq &
               D < (mean(thva_ordva[c(1:(ll - 1))]) - sd(thva_ordva[c(1:(ll - 1))]) *
                    3)))
          {
            if (A - B >= Dis.cq &
                A > (mean(thva_ordva[c(2:ll)]) + sd(thva_ordva[c(2:ll)]) * 3))
            {
              A <- mean(thva_ordva[c(2:ll)])
            }
            if (C - D >= Dis.cq &
                D < (mean(thva_ordva[c(1:(ll - 1))]) - sd(thva_ordva[c(1:(ll - 1))]) *
                     3))
            {
              D <- mean(thva_ordva[c(1:(ll - 1))])
            }
            thva[thva_order[1]] <- A
            thva[thva_order[ll]] <- D
          }
          else {
            thva <<- thva   #Outputting global variable
            break
          }
        }
      }
      for (ii in 1:num_group)
      {
        detect.col = c(((ii - 1) * bio_rep + 1):(ii * bio_rep))
        mean.col = apply(rwd[, detect.col], 1, mean, na.rm = T)   #Calculate the mean of same gene in one group
        for (jj in 1:(num_RG + num_TG))
        {
          thva = rwd[jj, detect.col]
          Eff.thva(thva)
          rwd[jj, detect.col] = thva  #Replace null value with corresponding mean value
        }
        rwd <<- rwd   #Outputting global variable
      }
      print("Eff.dpp.4 has been completed!")
    }
    ####################--Eff.dpp (main command)--####################
    Eff.dpp.1(rwd)
    Eff.dpp.2(rwd)
    Eff.dpp.3(rwd)
    Eff.dpp.2(rwd)
    it = 0
    repeat {
      Eff.dpp.4(rwd)
      mean_sample = apply(rwd, 2, mean, na.rm = T)
      mean_all = mean(mean_sample)
      cf = mean_sample - mean_all
      it = it + 1   #number of iteration
      if (as.numeric(sprintf("%0.5f", min(cf))) == 0 &
          as.numeric(sprintf("%0.5f", min(cf))) == 0)
      {
        message = paste0(
          " This step (Eff.dpp for qPCR data pre-processing) has all been completed!" ,
          "\n",
          " The number of iterations (Eff.dpp.2 to Eff.dpp.4) is ",
          it,
          "!",
          "\n"
        )
        cat(message)
        break
      }
      else {
        rwd = (t(t(rwd) - cf))
        print("Eff.dpp.2 has been completed!")
      }
    }
    rownames(rwd) = c(name_RG, name_TG)
    colnames(rwd) = c(paste0("sample", 1:(num_group * bio_rep)))
    cleandata <<-
      rwd                                                    #Outputting global variable
    write.csv(cleandata, "01-cleandata.csv")
    cleandata.RG <<-
      cleandata[(1:num_RG),]                           #Outputting global variable
    cleandata.TG <<-
      cleandata[(num_RG + 1):(num_TG + num_RG),]       #Outputting global variable
  }

  ###########################################################################
  ###########################  Step 2: Eff.gsg  #############################
  ###########################################################################
  Eff.gsg <-  function(cleandata.RG,
                       num_RG,
                       num_sample,
                       name_RG,
                       gsg.Ratio)
  {
    nc = ncol(cleandata.RG)
    nr = nrow(cleandata.RG)
    num_allRGcomb = sum(choose(num_RG, (1:num_RG)))      #number of RG combinations
    if (nc == num_sample)
    {
      if (nr == num_RG && nr == length(name_RG))
      {
        sdmean = data.frame()
        for (i in 1:nr) {
          if (i == 1) {
            b = t(combn(1:nr, i))
            for (j in 1:nrow(b)) {
              c = cbind(b, b)
              t = apply(cleandata.RG[c[j, ], ], 2, mean)                        #calculating the mean value of selected RGs within one sample
              t[nc + 1] = i
              t[nc + 2] = j
              t[nc + 3] = j
              t = as.data.frame(t(t))
              t[nc + 4] = paste(name_RG[b[j, ]], collapse = "+")
              sdmean = rbind(sdmean, t)
            }
          }
          else {
            b = t(combn(1:nr, i))
            for (j in 1:nrow(b)) {
              t = apply(cleandata.RG[b[j, ], ], 2, mean)                       #calculating the mean-value of selected RGs within one sample
              t[nc + 1] = i
              t[nc + 2] = j
              t[nc + 3] = sum(choose(num_RG, 1:(i - 1)), j)
              t = as.data.frame(t(t))
              t[nc + 4] = paste(name_RG[b[j, ]], collapse = "+")
              sdmean = rbind(sdmean, t)
            }
          }
        }
        sdmean$sd = apply(sdmean[, 1:nc], 1, sd)      #calculating the SD (standard deviation)
        colnames(sdmean) = c(
          paste0("Mean_sample", 1:nc),
          "Number of RG(s)",
          "Order(in same number of NFs)",
          "Order(all NFs)",
          "NF_name",
          paste("SD of all", num_sample, "mean values")
        )          #add the column names
        rownames(sdmean) = paste0("NF", 1:sum(choose(num_RG, (1:num_RG))))       #add the rownames
        data.NFs <<- sdmean         #Outputting global variable
        NF_index = paste0("NF", 1:sum(choose(num_RG, (1:num_RG))))                                    #creating the index file
        index <<-
          cbind(NF_index, NF_name = data.NFs[, (num_sample + 4)])                           #Outputting global variable
        write.csv(
          data.NFs,
          file = paste0(
            "02-The geometric value of all ",
            sum(choose(num_RG, 1:num_RG)),
            " Normalization factors (NFs) in all ",
            num_sample,
            " samples and their corresponding parameters.csv"
          ),
          row.names = TRUE
        )
        write.csv(index,
                  file = "Index file.csv",
                  row.names = TRUE)
        #extract the top gsg.Ratio as the gold standard group
        data.NFs_order = sdmean[order(sdmean[, (num_sample + 5)], decreasing = F), ]           #ordered by SD value
        data.NFs.gsgNo <<-
          data.NFs_order[1:(ceiling(num_allRGcomb * gsg.Ratio)), (num_sample + 3)]              #extract the row number in the top 1% SD
        data.NFs.gsg <<-
          data.NFs[data.NFs.gsgNo, ]                                            #Outputting global variable
        write.csv(
          data.NFs.gsg,
          file = paste0(
            "02-The gold standard group in the top ",
            gsg.Ratio * 100,
            "% SD values.csv"
          ),
          row.names = TRUE
        )
      }
      else {
        print(
          "The size of cleandata.RG matrix should match num_RG or length(name_RG), please check your cleandata.RG or your settings (num_RG, name_RG)!"
        )
      }
    }
    else {
      print(
        "The size of cleandata.RG matrix should match num_sample, please check your cleandata.RG or your settings (num_sample)!"
      )
    }
    message = paste(
      " =========================================================================================================",
      "\n",
      "This step (Eff.gsg for discovering the gold standard group) has all been completed!",
      "\n"
    )
    cat(message)
  }

  ###########################################################################
  ###########################  Step 3: Eff.rep  #############################
  ###########################################################################
  Eff.rep <- function(cleandata.TG,
                      data.NFs,
                      num_sample,
                      name_TG)
  {
    nf <- as.matrix(data.NFs[, 1:(num_sample)])
    tg <- as.matrix(cleandata.TG)
    nc1 = ncol(tg)
    nr1 = nrow(tg)
    nc2 = ncol(nf)
    nr2 = nrow(nf)
    if (nc1 == nc2)
    {
      if (nr1 == length(name_TG))
      {
        delta_ct = as.data.frame(array(, dim = c(nr2, 0)))                        #calculating the delta ct (cq) values
        for (ii in 1:nr1)   {
          delta_ct_single = t(tg[ii,] - t(nf[, (1:nc1)]))
          colnames(delta_ct_single) = c(paste0("Sample", 1:nc1, "_TG", ii, "(", name_TG[ii], ")"))        #add the column names
          delta_ct = cbind(delta_ct, delta_ct_single) }
        data.rep = 2 ^ (-delta_ct)                                                        #relative expression levels=2^(-??ct(TG-RG))
        rownames(data.rep) = paste0("NF", 1:sum(choose(num_RG, (1:num_RG))))                  #add the rownames
        data.rep <<- data.rep
        write.csv(
          data.rep,
          file = paste0(
            "03-The expression levels of the ",
            num_TG,
            " target genes were normalized by every one of the ",
            sum(choose(num_RG, 1:num_RG)),
            " Normalization factors.csv"
          ),
          row.names = T
        )
      }
      else {
        print(
          "Error: the length of name_TG did not match columns of cleandata.TG, please check your cleandata.TG or your settings (name_TG)!"
        )
      }
    }
    else {
      print(
        "Error: the number of columns of data.NFs did not match cleandata.TG, check your data!"
      )
    }
    message = paste(
      " =========================================================================================================",
      "\n",
      "This step (Eff.rep for obtaining the TG relative expression levels) has all been completed!",
      "\n"
    )
    cat(message)
  }

  ###########################################################################
  ###########################  Step 4: Eff.pv   #############################
  ###########################################################################
  Eff.pv <- function(data.rep, num_TG, num_group, bio_rep)
  {
    start_time <- Sys.time()
    library(tcltk)                                       #library a R package for showing progress status
    pb <- tkProgressBar("progress", "%completed", 0, 100)
    num_sample = num_group * bio_rep
    num_comp = choose(num_sample / bio_rep, 2)
    data.rep2 = data.rep
    nc2 = ncol(data.rep2) / num_TG
    nr2 = nrow(data.rep2)
    pairgp = t(combn(1:(nc2 / bio_rep), 2))             #paired groups
    nr_pairgp = nrow(pairgp)
    data.pv = data.frame()                       #create an empty matrix for all p-values
    data.fc = data.frame()                       #create an empty matrix for all p-values
    #record the program start time
    cnn = num_TG * choose(num_group, 2)
    rnn = sum(choose(num_RG, 1:num_RG))
    crnn = cnn * rnn
    message = paste(
      " =========================================================================================================",
      "\n",
      "This step (Eff.pv for calculating the p-value data matrix) is very time consuming!",
      "\n",
      paste0("          Your p-value data matrix size is ", cnn, " GCCs ?? ", rnn, " NFs,"),
      "\n",
      paste0("          Student's t-test needs to be performed ", crnn, " times!"),
      "\n"
    )
    cat(message)
    for (ii in 1:nr2)
    {
      info <-
        sprintf("Eff.pv program, %d%% completed", round(ii * 100 / nr2))                   #progress bar information
      setTkProgressBar(pb, ii * 100 / nr2, sprintf("progress (%s)", info), info)           #progress bar information
      p_value_row <-
        as.data.frame(array(, dim = c(1, 0)))                    #create an empty vector for every row
      fc_row <- as.data.frame(array(, dim = c(1, 0)))
      for (geneNo in 1:num_TG)
      {
        p_value_gene <- c()
        fc_gene <- c()
        for (mm in 1:nr_pairgp)
        {
          pos_group1_start = (geneNo - 1) * num_sample + (pairgp[mm, 1] - 1) * bio_rep + 1                                            #the start position of group1
          pos_group2_start = (geneNo - 1) * num_sample + (pairgp[mm, 2] - 1) *
            bio_rep + 1                                            #the start position of group2
          value_group1 = as.numeric(data.rep2[ii, pos_group1_start:(pos_group1_start + (bio_rep - 1))])        #loading relative expression values of group1
          value_group2 = as.numeric(data.rep2[ii, pos_group2_start:(pos_group2_start + (bio_rep - 1))])        #loading relative expression values of group2
          p_value_single <-
            t.test(value_group1, value_group2, var.equal = T)$p.value                          #t.test for every paired group
          p_value_gene <- c(p_value_gene, p_value_single)
          fc_single = mean(value_group1) / mean(value_group2)
          fc_gene <- c(fc_gene, fc_single)
        }
        p_value_row <- cbind(p_value_row, t(p_value_gene))
        fc_row <- cbind(fc_row, t(fc_gene))
      }
      data.pv = rbind(data.pv, t(t(p_value_row)))
      data.fc = rbind(data.fc, t(t(fc_row)))
    }
    rownames(data.pv) = paste0("NF", 1:sum(choose(num_RG, (1:num_RG))))       #add the rownames
    rownames(data.fc) = paste0("NF", 1:sum(choose(num_RG, (1:num_RG))))       #add the rownames
    colnames(data.pv) = paste(rep(c(paste(
      paste0("groups ", data.frame(pairgp)$X1),
      "vs",
      paste0(data.frame(pairgp)$X2),
      sep = ""
    )), times = num_TG), rep(paste0("in TG", 1:num_TG), each = num_comp))
    data.fc = -log2(data.fc)                                                  #calculate the -log2(FC) values
    colnames(data.fc) = paste(rep(c(paste(
      paste0("groups ", data.frame(pairgp)$X1),
      "vs",
      paste0(data.frame(pairgp)$X2),
      sep = ""
    )), times = num_TG), rep(paste0("in TG", 1:num_TG), each = num_comp))
    data.pv <<-
      data.pv                             #Outputting global variable
    write.csv(
      data.pv,
      file = paste0(
        "04-The p-value data matrix (",
        sum(choose(num_RG, 1:num_RG)),
        " NFs × ",
        choose(num_group, 2) * num_TG,
        " GCCs).csv"
      ),
      row.names = T
    )                     #write the p-value data matrix into a csv file
    write.csv(data.pv[data.NFs.gsgNo,],
              "04-The p-value data matrix (gold standard group).csv",
              row.names = T)
    #data.fc <<- data.fc                             #Outputting global variable
    #write.csv(data.fc, file = "04-The fold change data matrix.csv", row.names =T)      #write the fc (fold change) data matrix into a csv file
    #write.csv(data.fc[data.NFs.gsgNo, ],"04-The fold change data matrix of top 1% (gold standard).csv", row.names = T)
    close(pb)                                                                                      #close the progress bar
    end_time <- Sys.time()
    run_time <- difftime(end_time, start_time, units = "mins")
    print(paste(
      "The Eff.pv program running time:",
      round(as.numeric(run_time), 1),
      "mins"
    ))
    message = paste(
      " =========================================================================================================",
      "\n",
      "This step (Eff.pv for calculating the p-value data matrix) has all been completed!",
      "\n"
    )
    cat(message)
  }

  ###########################################################################
  ###########################  Step 5: Eff.ne   #############################
  ###########################################################################
  Eff.ne <- function(data.pv, num_RG, name_RG, pv.th)
  {
    ##########################code and calculate to obtain the gold standard############################
    pvdata2 = data.pv
    pvall_singleRGs = t(pvdata2[1:num_RG, ])                          #loading p-value of each single RG
    pv.th <<- pv.th
    (pvdata2[pvdata2 > pv.th[1]] <-
        2) &
      (pvdata2[pvdata2 < pv.th[1]] <-
         1)                                              #2 means negative (p>x), #1 means positive (p<x)
    pv_x_GS = apply(pvdata2[data.NFs.gsgNo, ], 2, mean)   #integrate the mean performance gold standard group at pv.th
    (pv_x_GS[pv_x_GS > 1.5] <-
        2) &
      (pv_x_GS[pv_x_GS <= 1.5] <-
         1)                                             #the mean of gold standard group, 1 (positive) or 2 (negative)
    gold_standard <<-
      pv_x_GS                                            #It is the gold standard
    tp_real = sum(gold_standard == 1)                      #positive frequency of gold standard
    tn_real = sum(gold_standard == 2)                      #negative frequency of gold standard
    #######ROC curve plotting, and the AUC values calculating ########
    library(pROC)
    library(grDevices)
    par(
      mfrow = c(1, 1),
      mar = c(4, 4, 1, 1),
      oma = c(0.5, 0.5, 0.5, 0.5)
    )
    Eff.ne.1 <- function(RG_order) {
      kk = length(RG_order)
      ramp <- colorRamp(c("red", "orange", "darkcyan", "blue"))
      colorpool <- rgb(ramp(seq(0, 1, length = kk)), max = 255)
      for (ii in 1:kk)
      {
        if (ii == 1) {
          rocobj <- plot.roc(
            gold_standard,
            pvall_singleRGs[, RG_order[ii]],
            percent = TRUE,
            col = colorpool[ii],
            grid = TRUE,
            smooth = TRUE,
            xlab = "FPR (%)",
            ylab = "TPR (%)",
            legacy.axes = TRUE,
            main = "Comparison in ROC space"
          )
          NE = round((rocobj$auc / 100 - 0.5) / 0.5 * 100, digit <-
                       4)
          pv_x =  pvall_singleRGs[, RG_order[ii]]
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
            which(RG_order == ii),
            paste0("NF", RG_order[ii]),
            name_RG[RG_order[ii]],
            NE,
            TPR,
            TNR,
            FPR,
            FNR,
            Precision,
            Recall
          )
          text(
            20,
            (70 - ii * 5),
            labels = paste0("——", name_RG[RG_order[ii]], " (NE) = ", sprintf("%0.2f", NE), "%"),
            col = colorpool[ii],
            cex = 1.2
          )
        }
        else {
          rocobj <-          lines.roc(
            gold_standard,
            pvall_singleRGs[, RG_order[ii]],
            percent = TRUE,
            col = colorpool[ii],
            grid = TRUE,
            smooth = TRUE
          )
          NE = round((rocobj$auc / 100 - 0.5) / 0.5 * 100, digit <-
                       4)
          pv_x =  pvall_singleRGs[, RG_order[ii]]
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
          NEdata_order2  = c(
            which(RG_order == ii),
            paste0("NF", RG_order[ii]),
            name_RG[RG_order[ii]],
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
            labels = paste0("——", name_RG[RG_order[ii]], " (NE) = ", sprintf("%0.2f", NE), "%"),
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
    RG_order1 = c(1:num_RG)
    Eff.ne.1(RG_order1)
    RG_order2 = as.numeric(NEdata_order[order(NEdata_order$NE, decreasing = T),][, 1])
    Eff.ne.1(RG_order2)
    NEdata_order$Ranking_order = c(1:length(RG_order2))
    NEdata_order = NEdata_order[,-1]
    NEdata_order <<- NEdata_order
    print(NEdata_order)
    message = paste(
      " =========================================================================================================",
      "\n",
      "This step (Eff.ne for detecting Normalization efficiencies of all single RGs) has all been completed!",
      "\n",
      "Now, you can use Eff.indexFinder function to detect the Normalization efficiency of any multi-RG combination.",
      "\n",
      "The index file is saved in your working directory!",
      "\n",
      "==========================================================================================================",
      "\n",
      "For more technical details about this algorithm, please see our publication:",
      "\n",
      "Jipan Zhang, Yongju Zhao: Normalization efficiency changes in reference genes (RG): evidence calling for use of multi-RG combination in qPCR experiments",
      "\n"
    )
    cat(message)
  }

  ####################################################################################
  ##################################   Main codes  ###################################
  ####################################################################################
  Eff.dpp(rawdata, num_RG, num_TG, num_group, bio_rep, Dis.cq, Primer.E)
  Eff.gsg(cleandata.RG, num_RG, num_sample, name_RG, gsg.Ratio)
  Eff.rep(cleandata.TG, data.NFs, num_sample, name_TG)
  Eff.pv(data.rep, num_TG, num_group, bio_rep)
  Eff.ne(data.pv, num_RG, name_RG, pv.th)
}
