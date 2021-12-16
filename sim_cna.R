# ---------------------- CNA simulation ----------------------------------------#

# ---------------------- parameters
# 1. Sample size (n)
# 2. Number of chromosomal sites or probes (m)
# 3. Number of recurrent CNAs (num_amp,num_del)
# 4. CN for each recurrent CNAs
# 5. Lengths of recurrent CNAs (len_amp,len_del)
# 6. Frequency of recurrent CNA gains and losses of the samples (freq_amp,freq_del)
# 7. Signal to noise ratio (SNR)
# 8. dir:文件保存路径
# Based on the above parameters, we simulate log2 intensity ratios
"simcna" <- function(r = 100, n = 100, m = 2000, interval = c(100, 149, 500, 529, 900 - 919), dir) {
     for (loop in 1:r) {
          # Initialization
          simu_data <- matrix(0, n, m)
          simu_data[, ] <- 2

          # Insert a number of recurrent CNA gains and losses at fixed positions across genome

          ap_s1 <- interval[1]
          ap_e1 <- interval[2]
          ap_s2 <- interval[3]
          ap_e2 <- interval[4]
          ap_s3 <- interval[5]
          ap_e3 <- interval[6]
          # C()向量
          freq_amp <- c(0.1, 0.15, 0.20)


          minsd <- 0.2
          masd <- 0.6

          # ----------------------------------------------- amplication ground truth
          num <- freq_amp[1] * n
          for (i in 1:num) {
               # runif()生成均匀分布随机数的函数   (n,min=0,max=1) n表示生成的随机数数量,min表示均匀分布的下限,max表示均匀分布的上限
               # 内置函数 round ()，它四舍五入到给定的位数，如果没有提供四舍五入的位数，它会将数字四舍五入到最接近的整数
               rs <- round(runif(1, min = 1, max = n))
               simu_data[rs, ap_s1:ap_e1] <- 4
          }

          num <- freq_amp[2] * n
          for (i in 1:num) {
               rs <- as.integer(runif(1, min = 1, max = n))
               simu_data[rs, ap_s2:ap_e2] <- 5
          }

          num <- freq_amp[3] * n
          for (i in 1:num) {
               rs <- as.integer(runif(1, min = 1, max = n))
               simu_data[rs, ap_s3:ap_e3] <- 6
          }



          #-----------------------------------------------------------------------------

          #------ noise 1, random alteration
          for (i in 1:100) {
               rand_s <- as.integer(runif(1, min = 1, max = n - 1))
               rand_pr <- as.integer(runif(1, min = 1, max = m - 100))
               CN <- runif(1, min = 3, max = 4)
               len <- round(runif(1, min = 10, max = 100))
               simu_data[rand_s, rand_pr:(rand_pr + len)] <- CN
          }


          #------------------------------------------------- mixing of tumor cells and normal cells
          for (i in 1:n) {
               rand_p <- runif(1, min = 0.3, max = 0.7)
               simu_data[i, ] <- log2((simu_data[i, ] * rand_p + 2 * (1 - rand_p)) / 2)
          }

          # ------------------------------------------------ add noise to the data
          for (i in 1:n) {
               rand_sd <- runif(1, min = minsd, max = masd)
               Gauss_noise <- rnorm(m, mean = 0, sd = rand_sd)
               simu_data[i, ] <- simu_data[i, ] + Gauss_noise
          }


          # ------------------------------------------------ write out
          simu_data <- t(simu_data)

          #-------------------------------------------------- used for cmds.
          header <- matrix(0, 1, n)
          header[1, ] <- 1:n
          simu_data_two <- (2^simu_data) * 2
          simu_data_two <- rbind(header, simu_data_two)
          header <- matrix(0, m + 1, 2)
          header[1, 1] <- "chromosome"
          header[1, 2] <- "position"
          header[2:(m + 1), 1] <- 1
          header[2:(m + 1), 2] <- 1:m
          simu_data_two <- cbind(header, simu_data_two)

          out.file_cmds <- paste(dir, "simu_data_", loop)
          write.table(simu_data_two, file = out.file_cmds, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
     }
}