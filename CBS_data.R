# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#       install.packages("BiocManager")
#   }
# BiocManager::install("DNAcopy")

library(DNAcopy) # gistic,SCEICG,stac
# denoised文件表示
"CBS_data" <- function(n = 10, dir) {
    data.dir <- dir
    ind <- c(1:n)
    for (i in ind) {
        # cat("The ", i, "th data", "\n")
        input.data.dir <- paste(data.dir, "simu_data_", i)
        data <- read.table(input.data.dir)
        # 去掉第一行，去掉前两列
        data <- data[-1, -(1:2)]

        # 最后一列Segment_Mean值就是 log2(copy_number/ 2), 正常来说人是二倍体生物则此value值为0，如果拷贝数小于2（删除）则小于0，拷贝数大于2（扩增）则大于0.
        # 参考https://www.cnblogs.com/leezx/p/7132099.html
        data <- log2(data / 2)
        # denoised文件，存储行号
        head <- matrix(0, nrow(data), 3)
        head[, 1] <- 1
        head[, 2] <- 1:nrow(data)
        head[, 3] <- 1:nrow(data)
        # 复制1两千次
        # rep函数主要是重复输出
        chrom <- rep(1, nrow(data))
        # 1-2000
        maploc <- 1:nrow(data)
        # 生成1行6列0矩阵
        seg.file_g <- matrix(0, 1, 6)
        seg.file_g_one <- matrix(0, 1, 6)
        seg.file <- matrix(0, nrow(data), 1)

        #---------------------------------------------------------------------------stac
        stac_amp <- matrix(0, 1, nrow(data))
        stac_amp[1, ] <- 1:nrow(data)
        stac_amp_one <- matrix(0, 1, nrow(data))

        stac_del <- matrix(0, 1, nrow(data))
        stac_del[1, ] <- 1:nrow(data)
        stac_del_one <- matrix(0, 1, nrow(data))
        # 遍历每个文件每列
        for (j in 1:ncol(data)) {
            # cat("sampl No,", j, "\n")
            # 从data取出第j列，chrom:2000，maploc:1-2000
            # 通常，CNV指种系事件，是种群中的变体。 CNA通常指的是体细胞事件，就像在肿瘤中发现的那样
            # CNA(genomdat,chrom,maploc)创建一个“拷贝数阵列”数据对象，用于通过循环二进制分割 (CBS) 等程序进行 DNA 拷贝数分析
            # genomdat:来自阵列 CGH、ROMA 或其他拷贝数实验的数据向量或矩阵。 如果是矩阵，则行对应于标记，列对应于样本
            # chrom:标记来自的染色体（或其他组标识符）。 长度与 genomdat 的行数相同的向量。 如果希望染色体按自然顺序排列，这个变量应该是数字或
            #      有序类别
            # maploc:标记在基因组上的位置。 长度与genomdat 的行数相同的向量。 这必须是数字。
            # segment()该程序使用循环二进制分割 (CBS) 将 DNA 拷贝数数据分割为估计相等拷贝数的区域
            # 该函数实现了 Olshen 和 Venkatraman (2004) 的圆形二元分割 (CBS) 算法。
            # 给定一组连续或二进制的基因组数据，该算法根据最大 t 统计量递归地将染色体分成两个或三个子段。
            # 用于决定是否拆分的参考分布是通过排列来估计的。当相邻段的均值相距不够远时，可以选择消除分裂。
            # 请注意，在第一次拆分后，用于拆分的测试的 α 水平不是无条件的。我们建议使用其中一个撤消选项来删除由于局部趋势
            # segment返回数据，表头解释参考https://www.jianshu.com/p/4312a453b4a4
            #        ID chrom loc.start loc.end num.mark seg.mean
            # 1 Sample.1     1         1    2000     2000  -0.0187
            seg <- segment(CNA(data[, j], chrom, maploc))

            n <- 0
            # length(seg$output$loc.start):返回seg的行数
            for (k in 1:length(seg$output$loc.start)) {
                seg.file_g_one[1, 1] <- j
                seg.file_g_one[1, 2] <- 1
                seg.file_g_one[1, 3] <- seg$output$loc.start[k] * 10000
                seg.file_g_one[1, 4] <- seg$output$loc.end[k] * 10000
                seg.file_g_one[1, 5] <- seg$output$num.mark[k]
                seg.file_g_one[1, 6] <- seg$output$seg.mean[k]
                seg.file_g <- rbind(seg.file_g, seg.file_g_one)
                seg.file_g_one <- matrix(0, 1, 6)

                len <- seg$output$num.mark[k]
                for (l in 1:len) {
                    n <- n + 1
                    # denoised文件第一列
                    seg.file[n, 1] <- seg$output$seg.mean[k]
                    # 生成另一组数据，是将序列区间聚焦到了基因上并评估基因是否是扩增或者缺失的
                    #------------------------------------------stac
                    if (as.numeric(seg$output$seg.mean[k]) > 0.1) {
                        stac_amp_one[1, n] <- 1
                    } else {
                        stac_amp_one[1, n] <- 0
                    }

                    if (as.numeric(seg$output$seg.mean[k]) < (-0.1)) {
                        stac_del_one[1, n] <- 1
                    } else {
                        stac_del_one[1, n] <- 0
                    }
                }
            }
            # cbind,rbind合并矩阵
            head <- cbind(head, seg.file)
            stac_amp <- rbind(stac_amp, stac_amp_one)
            stac_del <- rbind(stac_del, stac_del_one)
        }
        seg.file_g <- seg.file_g[-1, ]
        out.file <- paste(data.dir, "seg", i, sep = "")
        write.table(seg.file_g, file = out.file, row.names = F, col.names = F, quote = F, sep = "\t")

        out.file <- paste(data.dir, "data_denoised_", i, sep = "")
        write.table(head, file = out.file, row.names = F, col.names = F, quote = F, sep = "\t")

        #-----------------------------------------------------------------------------------stac

        firstLine <- stac_amp[1, ]
        firstLine[1] <- paste("\t", firstLine[1], sep = "")
        firstLine <- as.matrix(firstLine)
        firstLine <- t(firstLine)
        stac_amp <- stac_amp[-1, ]
        stac_del <- stac_del[-1, ]


        h <- matrix(0, nrow(stac_amp), 1)
        h[, 1] <- 1:nrow(stac_amp)
        stac_amp <- cbind(h, stac_amp)
        stac_del <- cbind(h, stac_del)


        out.file <- paste(data.dir, "\\", "data_stac", "\\", "data_stac_amp", i, sep = "")
        # write.table(firstLine, file = out.file, row.names = F, col.names = F, quote = F, sep = "\t")
        write.table(stac_amp, file = out.file, row.names = F, col.names = F, quote = F, sep = "\t")

        # out.file <- paste(data.dir, "\\", "data_stac", "\\", "data_stac_del", i, sep = "")
        # write.table(firstLine, file = out.file, row.names = F, col.names = F, quote = F, sep = "\t")
        # write.table(stac_del, file = out.file, row.names = F, col.names = F, quote = F, sep = "\t")
    }
}