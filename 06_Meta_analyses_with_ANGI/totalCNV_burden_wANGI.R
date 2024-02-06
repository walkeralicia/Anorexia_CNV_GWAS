

############################### Combining results with ANGI study #########################################################

#remotes::install_version("nloptr", version = "1.2.2.3")
ANGI <- read.csv("/afm01/UQ/Q4399/Anorexia/ANGI/results/angi_length_type.csv", header=TRUE)
UKB <- read.csv("/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/overall_burden/UKBB_length_type.csv", header=TRUE)
ANGI$Study <- "ANGI"
UKB$Study <- "UKB"

library(ggplot2)
library(ggpubr)
both <- rbind(ANGI, UKB)
both$index <- 1:dim(both)[1]
both$lower <- as.numeric(both$lower)
both$upper <- as.numeric(both$upper)
both$OR <- as.numeric(both$OR)
both$Pvalue <- as.numeric(both$Pvalue)
both$Pvalue <- signif(ifelse(both$OR>1, both$Pvalue/2, 1-(both$Pvalue/2)),3)
both$Type <- ifelse(both$Type=="Deletions and Duplications", "Both CNV Types", 
                    ifelse(both$Type=="Deletions", "Deletion-only", "Duplication-only"))
both$Type <- factor(both$Type, levels=c("Both CNV Types", "Deletion-only", "Duplication-only"))
write.table(both, "/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/overall_burden/both_length_type.txt", row.names=FALSE, col.names=TRUE, sep ="\t", quote=TRUE)


p1 <- ggplot(data=both, aes(x=OR, y=Length, color=Study)) +
  geom_errorbar(aes(xmin=lower, xmax = upper), width=0.4, linewidth = 6, position = position_dodge(width = 0.6)) +
  geom_point(size=12, position = position_dodge(width = 0.6))  +
  labs(x='OR', y = 'CNV Length', tag = "A") +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=2, linewidth=1.5) +
  xlim(0.5, 1.7) +
  scale_color_manual(values = c("ANGI" = "#AF58BA", "UKB" = "#F28522")) + 
  scale_y_discrete(limits = c("All", "20-100kB", "100-200kB", "200-500kB", ">500kB")) +
  theme_minimal(base_size=90) +
  theme(legend.position="none") + 
  scale_x_continuous(breaks = c(0.75, 1, 1.25)) +
  facet_grid(~Type)

png(filename = "/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/overall_burden/both_length_type.png",width = 50, height = 15, units='in', bg="white", res=150, type=c("cairo"))
print(p1)
dev.off()



### combine rare cnvs counts across ukb and angi
library(gridExtra)

ANGI <- read.csv("/QRISdata/Q4399/Anorexia/ANGI/results/ANGI_rare_count.csv", header=TRUE)
UKB <- read.csv("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden/rare/UKB_rare_counts.csv", header=TRUE)

UKB$lower <- as.numeric(UKB$lower)
UKB$upper <- as.numeric(UKB$upper)
UKB$OR <- as.numeric(UKB$OR)
UKB$Pvalue <- as.numeric(UKB$Pvalue)
UKB$Pvalue <- signif(ifelse(UKB$OR>1, UKB$Pvalue/2, 1-(UKB$Pvalue/2)),3)


p2 <- ggplot(data=UKB, aes(x=OR, y=CNV_Count)) +
  geom_errorbar(aes(xmin=lower, xmax = upper), width=0.4, linewidth = 6, color="#F28522") +
  geom_point(size=12, color="#F28522")  +
  labs(y='UKB rCNV Count', x = 'OR', tag = " ") +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=2, linewidth=1.5) +
  scale_y_discrete(limits = c("Singleton", "2-5", "6-10", "11-100", "101-1000","1001-1500", "1501-3872")) +
  theme_minimal(base_size=90) +
  xlim(c(0.7, 1.3))

ANGI$lower <- as.numeric(ANGI$lower)
ANGI$upper <- as.numeric(ANGI$upper)
ANGI$OR <- as.numeric(ANGI$OR)
ANGI$Pvalue <- as.numeric(ANGI$Pvalue)
ANGI$Pvalue <- signif(ifelse(ANGI$OR>1, ANGI$Pvalue/2, 1-(ANGI$Pvalue/2)),3)

p3 <- ggplot(data=ANGI, aes(x=OR, y=CNV_Count)) +
  geom_errorbar(aes(xmin=lower, xmax = upper), width=0.4, linewidth = 6, color="#AF58BA") +
  geom_point(size=12, color="#AF58BA")  +
  labs(y='ANGI rCNV Count', x = 'OR', tag = "B") +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=2, linewidth=1.5) +
  scale_y_discrete(limits = c("Singleton", "2-5", "6-10", "11-20", "21-40","41-80", "81-125")) +
  theme_minimal(base_size=90) +
  xlim(c(0.7, 1.3))

all_plots <- ggarrange(p1, 
                       ggarrange(p3, p2, nrow = 2), 
                       ncol = 2, widths=c(2,1))



png(filename = "/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/overall_burden/both_rare_counts.png",width = 60, height = 30, units='in', bg="white", res=150, type=c("cairo"))
print(all_plots)
dev.off()


### combining overall burden results

R
library(dplyr)
ANGI <- read.csv("/QRISdata/Q4399/Anorexia/ANGI/results/burden_results.csv", header=TRUE)
UKB <- read.csv("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden/burden_results.csv", header=TRUE)
ANGI$Study <- "ANGI"
UKB$Study <- "UKB"
colnames(ANGI) <- c("Type", "Test", "OR", "lower", "upper", "pval", "Avg_case", "Avg_cont", "Study")

both <- rbind(ANGI, UKB)
both$CI <- paste0("(", both$lower, "-", both$upper, ")")
both <- both %>% select(Type, Test, Study, OR, CI, pval, Avg_case, Avg_cont)
both$pval <- signif(ifelse(both$OR>1, both$pval/2, 1-(both$pval/2)),3)
both <- both[order(both$Type, both$Test),]
both$Type <- ifelse(both$Type=="ALL", "Deletions and Duplications", 
                    ifelse(both$Type=="DEL", "Deletions", "Duplications"))
both$Test <- ifelse(both$Test=="zooConsAm", "PhyloP", ifelse(both$Test=="zooConsPr", "PhastConst", both$Test))
write.table(both, "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/overall_burden/both_overall_burden.txt", row.names=FALSE, col.names=TRUE, sep ="\t", quote=TRUE)



