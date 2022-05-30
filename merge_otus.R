library(dplyr)
library(tidyr)
library(stringr)

otu_taxonomy <- read.csv("otu_taxonomy.txt", sep = '\t', header = F)
colnames(otu_taxonomy) <- c("X.OTU.ID", "V2", "V3", "V4")
otu_taxonomy <- otu_taxonomy %>% separate(V4, c("d", "p", "c", "o", "f", "g"), ",") 
otu_taxonomy <- data.frame(otu_taxonomy$X.OTU.ID, otu_taxonomy$c)
colnames(otu_taxonomy) <- c("X.OTU.ID", "taxa")

otus <- read.csv("egor_metagenome_otu.txt", sep='\t')


df <- round(as.data.frame(lapply(otus[,-c(1)],function(x)x/sum(otus[,-c(1)]))),7)

#test
test_otus <- otus[c(1:10), c(1:10)]
test_df <- round(test_otus[,-c(1)]/colSums(test_otus[,-c(1)]), 2)
test_df_2 <- round(as.data.frame(lapply(test_otus[,-c(1)],function(x)x/sum(test_otus[,-c(1)]))),2)

df_1 <- round(otus[,-c(1)]/colSums(otus[,-c(1)]), 7)
last_df <- cbind(otus[,c(1)], df_1)


colnames(last_df)[1] <- "X.OTU.ID"

merged <- merge(otu_taxonomy, last_df, by="X.OTU.ID")
write.table(merged, "merged_otus_per.tsv", sep='\t', row.names = F)



data_for_plot <- read.csv('merged_otus_per.tsv', sep='\t')



data_for_plot_1 <- read.csv('merged_otus_1.tsv', sep='\t')


test_data <- read.csv('table.csv', sep = ',', header = T)

library(ggplot2)

library(dplyr)
library(tidyr)
rownames(merged) <- merged$X.OTU.ID
data_long <- gather(test_data, taxa, sample, factor_key=FALSE)

library(data.table)
long <- melt(setDT(test_data), id.vars = c("X.OTU.ID","taxa"), variable.name = "SRR_ID")
metadata <- read.csv('metadata.csv', header = TRUE)

long_add <- merge(metadata, long, by = 'SRR_ID')
ggplot(data=long_add, aes(x=SRR_ID, y=value, fill=taxa)) +
  geom_bar(stat="identity") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(.~Healfy_status, scales="free_x", labeller = labeller(Healfy_status = c("TRUE" = "Healthy","FALSE" = "Sick")))

ggsave(file="healp.svg", plot=healp)

ggplot(wl_test, aes(x = Healfy_status, y = `diversity(diversity_data, index = "shannon", MARGIN = 2)`)) + 
  geom_boxplot() + 
  geom_jitter()



# plot places
long_add_place <- long_add[long_add$Healfy_status == FALSE, ]
s
ggplot(data=long_add_place, aes(x=SRR_ID, y=value, fill=taxa)) +
  geom_bar(stat="identity") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(.~Type.of.inflammation, scales="free_x", 
             labeller = labeller(Type.of.inflammation = c("left-sided colitis" = "left-sided","proctitis" = "proctitis", "total" = "total")))



diversity_data <- subset(otu_data, select = -X.OTU.ID)
div_by_healf <- as.data.frame(diversity(diversity_data, index = "shannon",  MARGIN = 2))
div_by_healf$SRR_ID <- rownames(div_by_healf)

wl_test <- merge(div_by_healf, metadata, by = 'SRR_ID')

wl_test_right <- wl_test[wl_test$Healfy_status == FALSE, ]
ggplot(wl_test_right, aes(x = Type.of.inflammation, y = `diversity(diversity_data, index = "shannon", MARGIN = 2)`)) + 
  geom_boxplot() + 
  geom_jitter()


left_sided <- wl_test_right[wl_test_right$Type.of.inflammation != "proctitis" & wl_test_right$Type.of.inflammation != "total",]
proctitis <- wl_test_right[wl_test_right$Type.of.inflammation == "proctitis",]
total_data <- wl_test_right[wl_test_right$Type.of.inflammation == "total",]

wilcox.test(left_sided$`diversity(diversity_data, index = "shannon", MARGIN = 2)`,proctitis$`diversity(diversity_data, index = "shannon", MARGIN = 2)`)
# W = 12, p-value = 0.414
wilcox.test(left_sided$`diversity(diversity_data, index = "shannon", MARGIN = 2)`,total_data$`diversity(diversity_data, index = "shannon", MARGIN = 2)`)
# W = 52, p-value = 0.917
wilcox.test(proctitis$`diversity(diversity_data, index = "shannon", MARGIN = 2)`,total_data$`diversity(diversity_data, index = "shannon", MARGIN = 2)`)
#W = 31, p-value = 0.4462
#plot_AR

data_AR <- read.csv('/media/lavrentydanilov/Linux_Home/lavrentydanilov/Documents/Documents/Scientific_work/ITMO/Scince_projects/metagenomes/Analysis/AP_plot.csv', header = TRUE)

ggplot(data_AR, aes(Antibiotics, Numer.of.resistant.isolates)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  coord_flip()+
  labs(x = 'Antibiotics', y = 'Number of resistant isolates')+
  theme_bw()


layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(MyResult.splsda.fixed, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 10)
plotLoadings(MyResult.splsda.fixed, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 10,  legend = FALSE, col.ties="black")
plotIndiv(MyResult.splsda.fixed, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, pch = 19, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "sPLS-DA ordination of samples", size.title = 1.5)
legend("bottomright", legend = levels(Y_1), cex = 1.5, fill = color.mixo(1:2), bty = "n")


layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(MyResult.splsda.fixed, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 10)
plotLoadings(MyResult.splsda.fixed, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 10,  legend = FALSE, col.ties="black")
plotIndiv(MyResult.splsda.fixed, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, pch = 19, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "sPLS-DA ordination of samples", size.title = 1.5)
legend("bottomright", legend = levels(Y_1), cex = 0.8, fill = color.mixo(1:3), bty = "n")