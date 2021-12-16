# @uppertime
# 需先执行sim_cna.R,CBS_data.R

# dir:文件绝对路径
# num:文件个数
# return:文件合并后的data

"read_data" <- function(dir, num) {
    data_all <- matrix(nrow = 0, ncol = 2000)
    for (i in 1:num) {
        cat("The ", i, "th data", "\n")
        input.dir <- paste(dir, "data_stac_amp", i, sep = "")
        data <- read.table(input.dir)
        data <- data[, -1]
        data_all <- rbind(data_all, data)
    }
    # out.file <- paste(dir, "test", sep = "")
    # write.table(data_all, file = out.file, row.names = F, col.names = F, quote = F, sep = "\t")
    return(data_all)
}
# 零假设:没有CNV突变
# data:样本
# permutation_n:置换次数
# alpha:显著性水平
# return：CNV的列
"evaluate" <- function(data, permutation_n = 100, alpha = 0.01) {
    # 求data每列均值,构建统计量,保存初始值
    pre_fi <- colMeans(data, na.rm = TRUE)
    # 构建统计每列置换样本大于原始样本统计量的零矩阵
    t_matrix <- matrix(0, nrow = 1, ncol = ncol(data))
    # 构造置换样本,行置换
    for (i in 1:permutation_n) {
        cat("The ", i, "th permutation", "\n")
        # 创建空矩阵，用以保存原始矩阵每一列都置换后的矩阵
        data_p <- matrix(nrow = 0, ncol = ncol(data))
        # 对每一行做置换,并求均值
        for (j in 1:nrow(data)) {
            # 置换每一行
            # cat("The ", j, "th row", "\n")
            fij_p <- sample(data[j, ])
            # 将置换后的每一行合并到置换后的矩阵
            data_p <- rbind(data_p, fij_p)
        }
        # 求取每列的均值
        # cat("The ", i, "th mean", "\n")
        fi_p <- colMeans(data_p, na.rm = TRUE)
        # 统计大于等于初始值的个数
        for (m in 1:length(fi_p)) {
            if (fi_p[m] >= pre_fi[m]) {
                t_matrix[1, m] <- t_matrix[1, m] + 1
            }
        }
    }
    # 计算每一列p值
    pfi <- t_matrix / permutation_n
    # 比较每一列是否接受原假设,并返回变异位点
    return(which(pfi <= alpha, arr.ind = TRUE)[, 2])
}
# data:原始数据
# pre_data:预测数据
"sensitivity" <- function(CNV_all, pre_CNVs) {
    tp <- length(intersect(CNV_all, pre_CNVs))
    sene <- tp / length(CNV_all)
    cat("sensitivity:", sene, "\n")
    return(sene)
}
# data:原始数据
# pre_data:预测数据
"precision" <- function(CNV_all, pre_CNVs) {
    prec <- length(intersect(CNV_all, pre_CNVs)) / length(pre_CNVs)
    cat("precision:", prec, "\n")
    return(prec)
}
# dir:根目录
# num:num*100条数据
# permutation_n：置换次数
# alpha：显著性水平
"main" <- function(dir = "D:\\data_experiment\\", num = 10, permutation_n = 1000, alpha = 0.01) {
    cat("处理规模：", (num * 100), "\n")
    cat("置换次数：", permutation_n, "\n")
    CNVs <- c(100, 149, 500, 529, 900, 919)
    CNV_all <- c()

    for (i in 1:length(CNVs)) {
        if (i %% 2 != 0) {
            CNV_all <- c(CNV_all, c(CNVs[i]:CNVs[i + 1]))
        }
    }
    cat("原始CNV片段位置：", CNV_all, "\n")
    cat("生成数据集中。。。。。\n")
    # 生成n条染色体，每条染色体m个位置，CNV插入位置为interval
    simcna(r = num, n = 100, m = 2000, interval = CNVs, dir = dir)

    cat("CBS算法处理数据集中。。。。。\n")
    CBS_data(n = num, dir = dir)

    # 读取num个文件，每个文件100条数据
    cat("读取文件中。。。。。\n")
    data <- read_data(paste(dir, "data_stac\\", sep = ""), num = num)
    cat("置换检验中。。。。。\n")
    pre_CNVs <- evaluate(data, permutation_n, alpha)
    cat("处理规模：", (num * 100), "\n")
    cat("置换次数：", permutation_n, "\n")
    cat("原始CNV片段位置：", CNV_all, "\n")
    cat("预测CNV片段位置：", pre_CNVs, "\n")
    sensitivity(CNV_all, pre_CNVs)
    precision(CNV_all, pre_CNVs)
    return()
}