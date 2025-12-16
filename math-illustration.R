library(igraph)
library(dplyr)
library(ggplot2)
library(ggsci)

set.seed(123)
rm(list = ls())

# -----------------------------
# 数据准备
# -----------------------------
math <- read.table("~/Downloads/sx-mathoverflow.txt") %>% as.data.frame()
math <- math[order(math[, 3]), ]
names(math) <- c("s", "r", "t")

all_nodes <- sort(unique(c(math$s, math$r)))
node_map <- setNames(1:(length(all_nodes)), all_nodes)
math_cleaned <- data.frame(
  s = unname(node_map[as.character(math$s)]),
  r = unname(node_map[as.character(math$r)]),
  t = math$t
)

n_lim = 24818
df <- subset(math_cleaned, s <= n_lim & r <= n_lim)

g <- graph_from_data_frame(df[, c("s", "r")], directed = TRUE)

# -----------------------------
# 参数
# -----------------------------
n_sim <- 20
step <- 500
track_nodes <- c(5, 15, 105, 1005, 10005)

# -----------------------------
# 完整图度数
# -----------------------------
deg_full <- degree(g)

# -----------------------------
# 主循环：采样与记录
# -----------------------------
all_records <- list()

for (sim in 1:n_sim) {
  records <- data.frame()
  for (v in track_nodes) {
    v_chr <- as.character(v)
    if (!(v_chr %in% as_ids(V(g)))) next
    
    other_nodes <- setdiff(as_ids(V(g)), v_chr)
    steps <- unique(c(seq(100, vcount(g), by = step), vcount(g)))
    
    for (k in steps) {
      chosen <- sample(other_nodes, size = min(k-1,length(other_nodes)), replace = FALSE)
      chosen <- c(v_chr, chosen)
      subg <- induced_subgraph(g, vids = chosen)
      
      deg <- degree(subg, v_chr)
      deg_norm <- ifelse(deg_full[v_chr] > 0, deg / deg_full[v_chr], NA)
      
      records <- rbind(records, data.frame(
        sim        = sim,
        node       = v,
        n          = k,
        degree     = deg,
        degree_norm= deg_norm
      ))
    }
  }
  message("Simulation ", sim, " finished")
  all_records[[sim]] <- records
}

records <- bind_rows(all_records)

library(ggrepel)

# 平均度数曲线 + log-log
avg_records_deg <- records %>%
  group_by(node, n) %>%
  summarise(mean_degree = mean(degree, na.rm=TRUE), .groups="drop")
c_estimates <- records %>%
  group_by(node, sim) %>%
  filter(degree > 0) %>%
  summarise(
    fit = list(lm(log(degree) ~ log(n))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(c = coef(fit)[2]) %>%
  ungroup()

c_summary <- c_estimates %>%
  group_by(node) %>%
  summarise(mean_c = mean(c), sd_c = sd(c))
# 用于标注 c 的数据：
# 取每个节点在最大 n 处的位置作为标注点
# 调整label_positions时直接把坐标右移
label_positions <- avg_records_deg %>%
  group_by(node) %>%
  filter(n == max(n)) %>%
  left_join(c_summary, by="node") %>%
  mutate(label = paste0("c==", round(mean_c, 2)),   # 字符串
         n = 0.9*n)                                     # 右移一些

p_loglog <- ggplot(avg_records_deg %>% filter(mean_degree > 0),
                   aes(x=n, y=mean_degree, color=factor(node))) +
  geom_line(size=0.5, c=0.8) +
  geom_point(size=0.8, c=0.8) +
  geom_text_repel(
    data = label_positions,
    aes(x = n, y = mean_degree, label = label, color = factor(node)),
    parse = TRUE,  
    nudge_y = 0.06,              # 让字符串转成表达式
    segment.color = "grey50",
    size = 3.5,
    show.legend = FALSE
  )+
  scale_x_log10(limits=c(1000, max(avg_records_deg$n)*1.3)) +  # 给右边留余地
  scale_y_log10() +
  scale_color_nejm(labels=c("1","2","3","4","5")) +
  labs(x="Network Size (log scale)", 
       y="Out degree (log scale)", 
       color="Nodes") +
  theme_minimal(base_size=10) +
  theme(plot.title = element_blank())+
  theme(legend.position = "right")

pdf(file = "math/p_loglog.pdf", width = 4, height = 3.5)
plot(p_loglog)
dev.off()

# =============================
# Rank–Degree log-log 图 + 横线 + LaTeX 标注
# =============================

deg <- degree(g, mode = "all")
n_total <- vcount(g)

deg_sorted <- sort(deg, decreasing = TRUE)
rank <- seq_along(deg_sorted)

df_rank <- data.frame(
  rank = rank,
  degree = deg_sorted,
  node = as.numeric(names(deg_sorted))
)

# 标记 track_nodes
df_rank <- df_rank %>%
  mutate(track_flag = ifelse(node %in% track_nodes, as.character(node), "Other"))

# NEJM 调色盘 (和 p_loglog 保持一致)
nejm_colors <- pal_nejm()(5)
names(nejm_colors) <- as.character(track_nodes)

# === 横线数据： y = n^{c} ===
hlines <- c_summary %>%
  filter(node %in% track_nodes) %>%
  mutate(
    y = n_total ^ mean_c,
    node = as.character(node),
    # latex 标签
    label = paste0("n^{", round(mean_c, 2), "}")
  )

# 作图
p_rank <- ggplot(df_rank, aes(x = rank, y = degree)) +
  geom_point(data = subset(df_rank, track_flag == "Other"),
             aes(color = track_flag),
             size = 0.8, c = 0.5, show.legend = FALSE) +
  geom_point(data = subset(df_rank, track_flag != "Other"),
             aes(color = track_flag),
             size = 2) +
  # 横线
  geom_hline(data = hlines,
             aes(yintercept = y, color = node),
             linetype = "dashed", show.legend = FALSE) +
  # 横线标注 (LaTeX 格式 + 不影响 legend)
  geom_text(
    data = hlines,
    aes(x = 10, y = y * 1.3, label = label, color = node),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 1.5, size = 3.5,
    show.legend = FALSE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = c(nejm_colors, Other = "grey70"),
    name   = "Nodes",
    breaks = as.character(track_nodes),
    labels = c("1","2","3","4","5")
  ) +
  labs(x = "Node rank (log scale)", y = "Out degree (log scale)") +
  theme_minimal(base_size = 10)+
  theme(legend.position = "right")

library(ggpubr)
pdf(file = "math/p_loglog_rank.pdf", width = 7, height = 3.5)
ggarrange(p_loglog,p_rank, common.legend = T, legend = "right")
dev.off()
pdf(file = "math/p_rank.pdf", width = 4, height = 3.5)
plot(p_rank)
dev.off()
