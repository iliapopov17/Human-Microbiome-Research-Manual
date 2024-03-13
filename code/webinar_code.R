################################
# Устанавливаем нужные пакеты
###############################

# install.packages(c("data.table",
#                    "openxlsx",
#                    "ggplot2",
#                    "vegan",
#                    "zCompositions",
#                    "devtools",
#                    "propr",
#                    "partitions",
#                    "data.tree",
#                    "e1071",
#                    "GUniFrac"
# ))
# 
# install_github("tpq/balance")
# install_bitbucket("knomics/nearestbalance")
# install_github("jbisanz/MicrobeR")



#####################
# Загружаем данные
#####################
library(data.table)
library(openxlsx)
# загружаем метаданные и упорядочиваем их по номеру участника
metadata = fread("data_Crohns_disease/metadata.csv")
metadata[, .N, by = diagnosis_full]

# Загружаем таблицу каунтов
counts <- read.csv("data_Crohns_disease/counts.csv",row.names = 1)
counts <- counts[metadata$sample, colSums(counts) >0]
dim(counts) # сколько образцов и микробов?
range(rowSums(counts)) # какое покрытие?
raw_abund <- counts/rowSums(counts) # состав образцов
library(MicrobeR)
metadata$SampleID <- metadata$sample
Microbiome.Barplot(t(raw_abund), metadata, CATEGORY = "diagnosis_full")

# достаточное покрытие?
metadata$coverage <- rowSums(counts)
library(ggplot2)
ggplot(metadata) + 
  geom_histogram(aes(coverage)) + 
  theme_minimal() + 
  xlab("N samples")

#####################################################
# Фильтрация от редких и малопредставленных таксонов
#######################################################
filt_counts <- counts[, colSums(counts>0)>0.3*nrow(counts)]
metadata[, filt_coverage := rowSums(filt_counts)]
metadata[, proportion_of_prevalent_taxa := 100*filt_coverage/coverage]

# покрытие образцов после фильтрации
ggplot(metadata)+
  geom_histogram(aes(filt_coverage)) + 
  theme_minimal()+
  xlab("покрытие после фильтрации") + 
  ylab("N образцов")
# какая доля микробов осталась в анализе
ggplot(metadata)+
  geom_histogram(aes(proportion_of_prevalent_taxa)) + 
  theme_minimal()+ 
  xlab("доля микробов, оставшихся в анализе") + 
  ylab("N образцов")

# считаем относительные  представленности
library(zCompositions)
library(NearestBalance)
abundance <- cmultRepl(filt_counts)
heatmap_with_split(abundance,
                   metadata, ~ diagnosis_full,
                   show_samp_names = F) + 
  theme(axis.text.y = element_text(size =5))

###############################
# Считаем альфа-разнообразие
#############################
# Считаем альфа-разнообразие
library(GUniFrac)
library(vegan)
alpha_div <- rowMeans(sapply(1:5, function(i){
  counts_rar_i = Rarefy(counts, 19000)$otu.tab.rff
  alpha_div_i = vegan::diversity(counts_rar_i)
}))
metadata$Shannon.index <- alpha_div[metadata$sample]

# Рисуем альфа-разнообразие
ggplot(metadata) + 
  geom_boxplot(aes(diagnosis_full, Shannon.index, fill = diagnosis_full)) + 
  theme_minimal() +
  theme(legend.position = 'none') + 
  xlab("") 
# оно поменялось?
wilcox.test(Shannon.index ~ diagnosis_full, metadata)$p.value

########################################
# бета-разнообразие Эйтчисона: есть ли различие по пропорциям?
########################################
clr <- log(abundance) - rowMeans(log(abundance))
beta_div <- dist(clr)

library(ape)
pcoa_res <- pcoa(beta_div)$vectors
var <- apply(pcoa_res, 2, var)
var_rel <- round(var*100/sum(var), 1)
ggplot(cbind(metadata,pcoa_res)) + 
  geom_point(aes(Axis.1, Axis.2, col=diagnosis_full)) +
  coord_fixed() + 
  theme_minimal() + 
  labs(col="") + 
  xlab(paste0("Axis.1 (",var_rel[1], "%)")) + 
  ylab(paste0("Axis.2 (",var_rel[2], "%)"))

# отличие статистически значимо?
adonis2(beta_div ~ metadata$diagnosis_full)

#####################################
# В чем именно состоит различие
#######################################
nb <- nb_lm(abundance = abundance,
            metadata = metadata,
            pred = "diagnosis_full")
heatmap_with_split(abundance = abundance[, unlist(nb$nb$b1)],
                   metadata = metadata,
                   formula = ~ diagnosis_full,
                   num_name = "ассоциированны со здоровьем",
                   den_name = "ассоциированны со болезнью",
                   show_samp_names = F,
                   balance = nb$nb$b1)
# рисуем пояснительную картиночку - различие между средней микробиотой здоровых и больных
psi <- make_psi_from_sbp(nb$coord$sbp)
mean_diff_clr <- drop(nb$lm_res$coefficients[2,] %*% psi)
bal_unit <- balance_to_clr(nb$nb$b1, colnames(abundance))
bal_diff_clr <- drop(bal_unit %*% mean_diff_clr) * bal_unit
tab <- data.table(taxon = names(mean_diff_clr),
                  clr_diff = mean_diff_clr,
                  bal = bal_diff_clr)
setorderv(tab, "clr_diff")
tab$taxon <- factor(tab$taxon, levels = tab$taxon)

# среднее различие между больными и здоровыми
ggplot(tab) + 
  geom_col(aes(clr_diff, taxon), fill = "darkblue") + 
  theme_minimal() + xlab("CLR(v)") + ylab("") + 
  theme(axis.text.y = element_text(size =5))

# приближенное, упрощенное различие
ggplot(tab) + 
  geom_col(aes(bal, taxon), fill = "darkblue") + 
  theme_minimal() + xlab("CLR(b)") + ylab("") +
  theme(axis.text.y = element_text(size =5))

# значение баланса в каждом образце
library(selbal)
metadata$balance <- balance.fromSBP(abundance, nb$nb$sbp)
ggplot(metadata) + 
  geom_boxplot(aes(diagnosis_full, balance, fill = diagnosis_full)) + 
  theme_minimal() +
  theme(legend.position = 'none') + 
  xlab("") 
# оно поменялось?
wilcox.test(balance ~ diagnosis_full, metadata)$p.value

