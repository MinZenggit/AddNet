# cat(1, "\n")
# source("test_simu_1.R")
# # 35.01474 7.80938 2.190632 0.4317256 0.1359093
# cat(2, "\n")
# source("test_simu_3.R")
# # 19.66861 4.795508 1.115385 0.2134667 0.07247049 
# cat(3, "\n")
# source("test_simu_4.R")  
# # 19.35373 4.333466 1.061593 0.19179 0.06726152
# cat(4, "\n")
# source("test_simu_5.R")
# # 19.38905 4.4096 1.067333 0.2006518 0.06547802


# 创建数据框
case1 <- c(35.01474, 7.80938, 2.190632, 0.4317256, 0.1359093)
case2 <- c(19.66861, 4.795508, 1.115385, 0.2134667, 0.07247049)
case3 <- c(19.35373, 4.333466, 1.061593, 0.19179, 0.06726152)
case4 <- c(19.38905, 4.4096, 1.067333, 0.2006518, 0.06547802)

# 创建数据框（长格式）
df <- data.frame(
  x = rep(c(2000, 1000, 500, 200, 100), 4),
  time = c(case1, case2, case3, case4),
  Case = factor(rep(1:4, each = 5))
)

# 加载ggplot2包
library(ggplot2)

# 创建散点图并添加二次拟合曲线
shapes <- c(16, 17, 18, 15)   # 不同的点形状：圆形，三角形，菱形，方形
linetypes <- c("dashed", "dashed", "dotted", "dotdash")  # 不同的线型

# 创建散点图并添加二次拟合曲线
p <- ggplot(df, aes(x = x, y = time, color = Case, shape = Case, linetype = Case)) +
  geom_point(alpha = 0.9, size = 3) +  # 散点
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = 0.6,
              se = FALSE, size = 1) +  # 二次拟合，使用不同的线型
  scale_color_manual(values = c("#377EB8", "#4DAF4A", "#E41A1C", "#984EA3")) +
  scale_shape_manual(values = shapes) +
  scale_linetype_manual(values = linetypes) +
  labs(x = "Number of Nodes", y = "Computation Time (s)", 
       title = NULL) +
  theme_minimal() +
  theme(legend.position = "right")

# 显示图形
print(p)
# 保存图形（可选）
ggsave("result/test_simu_all_cases/simu_time.pdf", plot = p, width = 5, height = 4)
