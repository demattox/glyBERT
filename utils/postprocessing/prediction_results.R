
library(reshape2)
library(ggplot2)
library(maptools)

library(ggpubr)
library(car)

######################################

coff_etal_aucs = as.data.frame(matrix(0, nrow = 20, ncol = 3))
colnames(coff_etal_aucs) = c('train','validate', 'shared')
coff_etal_aucs$shared = T

row.names(coff_etal_aucs) = c('ABA',
                              'ConA',
                              'DBA',
                              'DC-SIGN',
                              'GSL.I.B4',
                              'H1N1',
                              'H3N8',
                              'jacalin',
                              'LCA',
                              'MAL.I',
                              'MAL.II',
                              'PHA.E',
                              'PHA.L',
                              'PNA',
                              'PSA',
                              'RCA.I',
                              'SBA',
                              'SNA.I',
                              'UEA.I',
                              'WGA')

# coff_etal_aucs['DC-SIGN',3] = F
# coff_etal_aucs['H1N1',3] = F
# coff_etal_aucs['H3N8',3] = F

coff_etal_aucs$validate = c(0.934,
                         0.971,
                         0.839,
                         0.841,
                         0.867,
                         0.913,
                         0.959,
                         0.882,
                         0.964,
                         0.833,
                         0.718,
                         0.959,
                         0.914,
                         0.914,
                         0.890,
                         0.953,
                         0.875,
                         0.950,
                         0.861,
                         0.882)

coff_etal_aucs$train = c(0.947,
                         0.982,
                         0.897,
                         0.955,
                         0.953,
                         0.973,
                         0.958,
                         0.896,
                         0.976,
                         0.848,
                         0.814,
                         0.975,
                         0.967,
                         0.943,
                         0.929,
                         0.958,
                         0.938,
                         0.979,
                         0.895,
                         0.901)



mean(coff_etal_aucs$validate)
mean(coff_etal_aucs$validate[coff_etal_aucs$shared])

# write.csv(coff_etal_aucs, file = '~/Desktop/Coff_etal_AUCs.csv', quote = F)


bowen_aucs = read.delim('~/cbk/sugarTrees/glyPredict/cv_updated_new.csv', sep = ',', stringsAsFactors = F)
row.names(bowen_aucs) =  bowen_aucs$X
bowen_aucs$X <- NULL

# prev = read.delim('~/Downloads/cv5.csv', sep = ',', stringsAsFactors = F)
# prev = read.delim('~/cbk/sugarTrees/glyPredict/cv5_updated.csv', sep = ',', stringsAsFactors = F)

# row.names(prev) =  prev$X
# prev$X <- NULL
# 
# for(i in 1:nrow(prev)){
#   lec = row.names(prev)[i]
#   if (lec %in% row.names(bowen_aucs)){
#     cat(lec,'\n\tOld: ', prev$X0[i], '\n\tNew: ', bowen_aucs$X0[row.names(bowen_aucs) == lec], '\n')
#   }
# }
# prev = prev[row.names(bowen_aucs),]
# prev = prev[!is.na(prev$X0),]
# tag = row.names(bowen_aucs) %in% row.names(prev)
# all(row.names(bowen_aucs)[tag] == row.names(prev))
# plot(bowen_aucs$X0[tag], prev$X0,
#      xlab = 'New AUC vals', ylab = 'Prev AUC vals',
#      xlim = c(0.80,1), ylim = c(0.80,1),
#      pch = 19, col = alpha('black', 0.5))
# abline(a = 0, b=1)
# text(x = bowen_aucs$X0[tag], y = prev$X0, labels = row.names(prev))
# # up-down error bars
# arrows(y0 = prev$X0 - prev$X1, y1 = prev$X0 + prev$X1, x0 = bowen_aucs$X0[tag], x1 =  bowen_aucs$X0[tag],
#        code=3, angle=90, length=0.05, col="black", lwd=1)
# #left right error bars
# arrows(x0 = bowen_aucs$X0[tag] - bowen_aucs$X1[tag], x1 =  bowen_aucs$X0[tag] + bowen_aucs$X1[tag], y0 = prev$X0, y1 = prev$X0,
#        code=3, angle=90, length=0.025, col="black", lwd=.5)


row.names(bowen_aucs)[! row.names(bowen_aucs) %in% row.names(coff_etal_aucs)]

row.names(bowen_aucs)[row.names(bowen_aucs) == 'GSL'] = 'GSL.I.B4'

row.names(bowen_aucs)[row.names(bowen_aucs) == 'MAL-I'] = 'MAL.I'
row.names(bowen_aucs)[row.names(bowen_aucs) == 'MAL_II'] = 'MAL.II'

row.names(bowen_aucs)[row.names(bowen_aucs) == 'PHA-E'] = 'PHA.E'
row.names(bowen_aucs)[row.names(bowen_aucs) == 'PHA-L'] = 'PHA.L'

row.names(bowen_aucs)[row.names(bowen_aucs) == 'UEAI'] = 'UEA.I'

row.names(bowen_aucs)[row.names(bowen_aucs) == 'SNA'] = 'SNA.I'

row.names(bowen_aucs)[row.names(bowen_aucs) == 'RCAI'] = 'RCA.I'

# row.names(bowen_aucs)[row.names(bowen_aucs) == 'jacalin'] = 'jacalin'

row.names(bowen_aucs)[row.names(bowen_aucs) == 'HA'] = 'H1N1'

row.names(bowen_aucs)[row.names(bowen_aucs) == 'Human'] = 'DC-SIGN'

all(row.names(bowen_aucs) %in% row.names(coff_etal_aucs))

coff_etal_aucs = coff_etal_aucs[row.names(bowen_aucs),]

row.names(coff_etal_aucs) == row.names(bowen_aucs)

coff_etal_aucs$val_sd = c(0.061, # gslib4
                          0.028, #H3N8
                          0.049, # UEA I
                          0.032, #lca
                          0.060, # sna
                          0.048, #PNA
                          0.026, # rca.i
                          0.061, #sba
                          0.018, #PHA.E
                          0.126, #PHAL.L
                          0.021, # wga
                          0.055, #jacalin
                          0.031, # conA
                          0.078, #MAL-II
                          0.035, # MAL-I
                          0.034, #aba
                          0.105, #H1N1
                          0.069, #dba
                          0.053, #psa
                          0.062) # DC-sign

colnames(bowen_aucs) = c('auc', 'sd')

bowen_aucs$coff_Val = coff_etal_aucs$validate
bowen_aucs$coff_sd = coff_etal_aucs$val_sd

mean(bowen_aucs$auc)
median(bowen_aucs$auc)



# bowen_aucs$diff = abs(bowen_aucs$auc - bowen_aucs$coff_Val)
# min(bowen_aucs[,c(1,3)])
# 
# 
# color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') { # https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
#   scale = (length(lut)-1)/(max-min)
#   
#   dev.new(width=1.75, height=5)
#   plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
#   axis(2, ticks, las=1)
#   for (i in 1:(length(lut)-1)) {
#     y = (i-1)/scale + min
#     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
#   }
# }
# 
# cols = colorspace::heat_hcl(250, h=c(0,-100), l=c(75,40), c=c(40,80), power = 1)
# summary(bowen_aucs$diff)
# breaks = seq(0,0.13, length.out = 250)
# 
# bowen_aucs$colDiffs = ''
# for (i in 1:nrow(bowen_aucs)){
#   bowen_aucs$colDiffs[i] = cols[breaks > bowen_aucs$diff[i]][1]
# }
# 
# 
# color.bar(cols, min = 0, max = 0.13)
# 
# 
# plot(bowen_aucs$auc, bowen_aucs$coff_Val,
#      xlim = c(0.71,1.01), ylim = c(0.71,1.01),
#      pch = 19, col = bowen_aucs$colDiffs,
#      xlab = "",
#      ylab = "", 
#      cex = 1.5, bty = 'n')
# # abline(a=0, b=1, lty =2)
# lines(x = c(0,1), y = c(0,1), lty = 2)
# lines(x = c(1,1), y = c(0,1), lwd = 0.75)
# lines(x = c(0,1), y = c(1,1), lwd = 0.75)
# title(main = 'Prediction of lectin binding from microarray data', cex.main = 2)
# title(ylab='CCARL Mean Validation AUC',
#       col.lab = 'black', cex.lab = 1.5, font = 2)
# title(xlab='glyBERT Mean Validation AUC',
#       col.lab = 'black', cex.lab = 1.5, font = 2)
# #Coff error bars
# # arrows(y0 = bowen_aucs$coff_Val - bowen_aucs$coff_sd, y1 = bowen_aucs$coff_Val + bowen_aucs$coff_sd, x0 = bowen_aucs$auc, x1 = bowen_aucs$auc,
# #        code=3, angle=90, length=0.05, col="dodgerblue", lwd=1)
# # #Bert error bars
# # arrows(x0 = bowen_aucs$auc - bowen_aucs$sd, x1 = bowen_aucs$auc + bowen_aucs$sd, y0 = bowen_aucs$coff_Val, y1 = bowen_aucs$coff_Val,
# #        code=3, angle=90, length=0.025, col="forestgreen", lwd=.5)
# # text(x = bowen_aucs$auc, y = bowen_aucs$coff_Val, labels = row.names(bowen_aucs), cex = 1.5, font = 2)
# pointLabel(x = bowen_aucs$auc, y = bowen_aucs$coff_Val, labels = row.names(bowen_aucs), cex = 1.5, font = 2,
#            allowSmallOverlap = F)
#         # ylim = c(0.8,1))



bowen_aucs$GlyNet = 0
bowen_aucs['MAL.II', 'GlyNet'] = 0.706
bowen_aucs['DC-SIGN', 'GlyNet'] = 0.805
bowen_aucs['UEA.I', 'GlyNet'] = 0.824
bowen_aucs['DBA', 'GlyNet'] = 0.838
bowen_aucs['jacalin', 'GlyNet'] = 0.848
bowen_aucs['MAL.I', 'GlyNet'] = 0.869
bowen_aucs['H1N1', 'GlyNet'] = 0.918
bowen_aucs['SBA', 'GlyNet'] = 0.923
bowen_aucs['PHA.L', 'GlyNet'] = 0.941
bowen_aucs['WGA', 'GlyNet'] = 0.941
bowen_aucs['ABA', 'GlyNet'] = 0.943
bowen_aucs['PNA', 'GlyNet'] = 0.952
bowen_aucs['H3N8', 'GlyNet'] = 0.954
bowen_aucs['PSA', 'GlyNet'] = 0.955
bowen_aucs['PHA.E', 'GlyNet'] = 0.956
bowen_aucs['GSL.I.B4', 'GlyNet'] = 0.963
bowen_aucs['RCA.I', 'GlyNet'] = 0.969
bowen_aucs['SNA.I', 'GlyNet'] = 0.969
bowen_aucs['LCA', 'GlyNet'] = 0.974
bowen_aucs['ConA', 'GlyNet'] = 0.996



tmp = bowen_aucs[,c(1,3,5)]
tmp = tmp[order(tmp$auc, decreasing = T),]
colnames(tmp) = c('glyBERT', 'CCARL', 'GlyNet')
tmp = tmp[,c(1,3,2)]


barplot(t(as.matrix(tmp)), beside = T, ylim = c(0.7, 1)) # Barplot option
plot(x = c(jitter(rep(1,20)), rep(2,20), rep(3,20)), y = c(tmp$glyBERT, tmp$GlyNet, tmp$CCARL)) # Lineplot option

cols = colorspace::rainbow_hcl(20)
tmpJitter = cbind(tmp, data.frame(x1 = jitter(rep(1,20), amount = 0.05), x2 = jitter(rep(2,20), amount = 0.05), x3 = jitter(rep(3,20), amount = 0.05)))

# plot(0,0,
#      xlim = c(0.8,3.1), ylim = c(0.7,1),
#      type = 'n', xaxt = 'n',
# 
# par(new = T)
boxplot(tmp, at = c(1,2,3), notch = F,
        xlim = c(0.5,3.5), ylim = c(0.7,1),
        outline = F,
        xlab = 'Method', 
        ylab = 'Mean AUC')
for(i in 1:20){
  par(new=T)
  plot(x = as.numeric(tmpJitter[i,4:6]), y = tmpJitter[i,1:3],
       xlim = c(0.5,3.5), ylim = c(0.7,1),
       axes = F, xlab = '', ylab = '',
       type = 'b', col = cols[i], pch = 19)
  pointLabel(x = as.numeric(tmpJitter[i,4:6]), y = tmpJitter[i,1:3], labels = row.names(tmpJitter)[i], cex = 0.7, font = 2, col = cols[i],
             allowSmallOverlap = F)
}
axis(side=1,at=c(1,2,3), labels = colnames(tmp))


aucMelt = as.data.frame(matrix(0, nrow = nrow(bowen_aucs)*3, ncol = 4))
colnames(aucMelt) = c('Lectin', 'Method', 'AUC', 'sd')

aucMelt$Lectin = rep(row.names(bowen_aucs),3)

aucMelt$Method[1:20] = 'glyBERT'
aucMelt$AUC[1:20] = bowen_aucs$auc
aucMelt$sd[1:20] = bowen_aucs$sd

aucMelt$Method[21:40] = 'GlyNet'
aucMelt$AUC[21:40] = bowen_aucs$GlyNet
aucMelt$sd[21:40] = NA

aucMelt$Method[41:60] = 'CCARL'
aucMelt$AUC[41:60] = bowen_aucs$coff_Val
aucMelt$sd[41:60] = bowen_aucs$coff_sd


mean(bowen_aucs$auc)
median(bowen_aucs$auc)

mean(bowen_aucs$GlyNet)
median(bowen_aucs$GlyNet)

mean(bowen_aucs$coff_Val)
median(bowen_aucs$coff_Val)

qqnorm(aucMelt$AUC)
qqline(aucMelt$AUC)
qqPlot(aucMelt$AUC) # Some evidence of non-normality

t.test(bowen_aucs$auc, bowen_aucs$GlyNet)
t.test(bowen_aucs$auc, bowen_aucs$coff_Val)
t.test(bowen_aucs$GlyNet, bowen_aucs$coff_Val)

oneWay = aov(data = aucMelt, formula = AUC ~ Method)
summary(oneWay)
tukey.oneWay = TukeyHSD(oneWay)
tukey.oneWay

plot(oneWay, 1)
leveneTest(AUC ~ Method, data = aucMelt) # No significant differences in variance

plot(oneWay, 2)
resis = residuals(oneWay)
shapiro.test(x = resis) # Shapiro-Wilk test

shapiro.test(aucMelt$AUC)

# Significant non-normality, use non-parametric method --> 

ft = friedman.test(AUC ~ Method | Lectin, data = aucMelt) # Matched sets of measurements, Friedman test (can be thought of as non-parametric analog of repeated measures ANOVA)
ft
# post-hoc with paired Wilcoxon signed-rank tests and Bonferonni adjustment
bc = wilcox.test(aucMelt$AUC[aucMelt$Method == 'glyBERT'], aucMelt$AUC[aucMelt$Method == 'CCARL'], paired = T) # post-hoc with paired Wilcoxon signed-rank test
bg = wilcox.test(aucMelt$AUC[aucMelt$Method == 'glyBERT'], aucMelt$AUC[aucMelt$Method == 'GlyNet'], paired = T) # post-hoc with paired Wilcoxon signed-rank test
gc = wilcox.test(aucMelt$AUC[aucMelt$Method == 'GlyNet'], aucMelt$AUC[aucMelt$Method == 'CCARL'], paired = T) # post-hoc with paired Wilcoxon signed-rank test

pwc.p = c(bc$p.value, bg$p.value, gc$p.value) 
names(pwc.p) = c('glyBERT-CCARL', 'glyBERT-GlyNet', 'GlyNet-CCARL')
pwc.p.adj = p.adjust(pwc.p, method = 'holm')

# with(aucMelt , boxplot( AUC ~ Method ))
# friedman.test.with.post.hoc(AUC ~ Method | Lectin, data = aucMelt)

boxplot(bowen_aucs$auc - bowen_aucs$coff_Val)
boxplot(bowen_aucs$auc - bowen_aucs$GlyNet)
boxplot(bowen_aucs$GlyNet - bowen_aucs$coff_Val)

aucMelt$Lectin = factor(aucMelt$Lectin, levels = row.names(bowen_aucs)[order(bowen_aucs$auc, decreasing = T)])

aucMelt$Method = factor(aucMelt$Method, levels = c('glyBERT', 'GlyNet', 'CCARL'))

# p1 = ggplot(aucMelt, aes(Lectin, y = AUC, fill = Method)) + geom_bar(position = position_dodge(0.5), width = .5, stat = 'identity') + theme_minimal()
# p1
# dev.off()
# p1 + geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), colour="black", width=.1, position=position_dodge(0.5)) +
#   coord_cartesian(ylim = c(0.65,1)) +
#   scale_fill_manual(values = c('dodgerblue1', 'gold2', 'plum3')) +
#   labs(x = 'Lectins', y = '5x CV mean AUC') + 
#   theme(text = element_text(size = 20), axis.text = element_text(size = 12)) 
# 
# ggsave(filename = 'lectinAUCs.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
#        width = 18,
#        height = 5.4)
# 
# 
# dev.off()
# p2 = ggplot(aucMelt, aes(y = Method, x = AUC, fill = Method)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width=0, height = 0.1, alpha=0.7, size = 2) +
#   theme_minimal() +
#   scale_fill_manual(values = alpha(c('dodgerblue1', 'gold2', 'plum3'), 0.6)) +
#   stat_summary(fun=mean, geom="point", shape='|', size=18, color=alpha("red", 0.5), show.legend = F) +
#   scale_x_reverse(limits = c(1, 0.7)) +  
#   scale_y_discrete(limits = c('CCARL', 'GlyNet', 'glyBERT')) +
#   labs(x = '5x CV mean AUCs from all lectins', y = '') + 
#   theme(text = element_text(size = 20)) 
# 
# p2
# ggsave(filename = 'lectinAUC_boxplots.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
#        width = 18,
#        height = 4)

# Reverse order?
aucMelt$Lectin = factor(aucMelt$Lectin, levels = row.names(bowen_aucs)[order(bowen_aucs$auc, decreasing = T)])

p3 = ggplot(aucMelt, aes(Lectin, y = AUC, fill = Method)) + geom_bar(position = position_dodge(0.65), width = .65, stat = 'identity',color = 'grey33', size=.75) + theme_minimal()
p3
dev.off()
p3 + geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), colour="black", width=.1, position=position_dodge(0.65)) +
  # coord_flip(ylim = c(0.65,1)) + 
  coord_cartesian(ylim = c(0.65,1)) +
  scale_fill_manual(values = c('#16a3f5', '#f5d416', '#f5165d')) +
  labs(x = 'Lectins', y = '5x CV mean AUC') + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'bold'),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey")) 

ggsave(filename = 'lectinAUCs.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 7.5*1.5,
       height = 6*1.5)
ggsave(filename = 'lectinAUCs_wide.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 9*1.5,
       height = 6*1.5)
ggsave(filename = 'lectinAUCs_xtraWide.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 11.5*1.5,
       height = 6*1.5)
ggsave(filename = 'lectinAUCs_short_xtraWide.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 11.5*1.5,
       height = 4*1.5)


# p-values
pwc.p.adj
stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,    ~p.ast,
  "GlyNet",     "CCARL", 1.707e-01, 'ns', # pwc.p.adj['GlyNet-CCARL']
  "glyBERT",     "GlyNet", 3.050e-03, '**', # pwc.p.adj['glyBERT-GlyNet']
  "glyBERT",     "CCARL", 1.717e-05, '***' # pwc.p.adj['glyBERT-CCARL']
)
stat.test

dev.off()
lecCols = colorspace::rainbow_hcl(20)
p4 = ggplot(aucMelt, aes(x = Method, y = AUC, group = Lectin)) + geom_jitter(width=0.0, height = 0, alpha=0.7, size = 2) +
  # scale_fill_manual(values = alpha(c('dodgerblue1', 'gold2', 'plum3'), 0.6)) +
  geom_path(aes(color = Lectin),
            alpha = 0.6, size = 2.5,
            lineend = 'round', linejoin = 'round') +
  # scale_colour_manual(values = lecCols) +
  scale_color_hue(l=80, c=120) + 
  theme_minimal() +
  coord_cartesian(ylim = c(0.7,1.05)) +
  # scale_x_reverse(limits = c(1, 0.7)) +  
  # scale_x_discrete(limits = c('CCARL', 'GlyNet', 'glyBERT')) +
  labs(y = '5x CV mean AUCs', y = 'Method') + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(color = c('#16a3f5', '#f5d416', '#f5165d'),
                                   face = 'bold',
                                   size = 20),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey"))

p4
pwc.p.adj
ft
p4 = p4 + stat_pvalue_manual(stat.test, y.position = 1.005, step.increase = 0.04, label = "p.ast", tip.length = 0.01, label.size = 7) +
  annotate("text", x=1.5, y=1.06, label="Friedman, p = 3.04e-0.5", size = 6)

meds =  c(median(aucMelt$AUC[aucMelt$Method == 'glyBERT']),
          median(aucMelt$AUC[aucMelt$Method == 'GlyNet']),
          median(aucMelt$AUC[aucMelt$Method == 'CCARL']))
medShape = 95
medSize = 16
# medCols = c('#16a3f5', '#f5d416', '#f5165d')
medCols = rep('red',3)

p4 + geom_point(aes(x = 1, y = meds[1]), shape = medShape, size = medSize, color = medCols[1]) + 
  geom_point(aes(x = 2, y = meds[2]), shape = medShape, size = medSize, color = medCols[2]) + 
  geom_point(aes(x = 3, y = meds[3]), shape = medShape, size = medSize, color = medCols[3])

ggsave(filename = 'lectinAUC_parallelCoord.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 4*1.5,
       height = 7.5*1.5)
ggsave(filename = 'lectinAUC_parallelCoord_wide.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 5*1.5,
       height = 7.5*1.5)
ggsave(filename = 'lectinAUC_parallelCoord_wider.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 5.5*1.5,
       height = 7.5*1.5)

#####################
# Put taxonomic performance in figure

taxo = as.data.frame(matrix(0,nrow = 6, ncol = 4))
row.names(taxo) = c('Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family')
colnames(taxo) = c('SweetTalk', 'SweetNet', 'Mohapatra et al.', 'glyBERT')
taxo = as.data.frame(t(taxo))

taxo$Domain = c(0.931, 0.942, 0.940, 0.945)
taxo$Kingdom = c(0.895, 0.910, 0.921, 0.936)
taxo$Phylum = c(0.801, 0.847, 0.808, 0.872)
taxo$Class = c(0.715, 0.757, 0.724, 0.794)
taxo$Order = c(0.533, 0.603, 0.574, 0.632)
taxo$Family = c(0.466, 0.554, 0.544, 0.557)

taxMelt = as.data.frame(matrix(0, nrow = ncol(taxo) * nrow(taxo), ncol = 3))
colnames(taxMelt) = c('Taxonomic Level', 'Method', 'Accuracy')

taxMelt$Method = rep(row.names(taxo), 6)

taxMelt$`Taxonomic Level`[1:4] = colnames(taxo)[1]
taxMelt$Accuracy[1:4] = taxo[,1]
taxMelt$`Taxonomic Level`[5:8] = colnames(taxo)[2]
taxMelt$Accuracy[5:8] = taxo[,2]
taxMelt$`Taxonomic Level`[9:12] = colnames(taxo)[3]
taxMelt$Accuracy[9:12] = taxo[,3]
taxMelt$`Taxonomic Level`[13:16] = colnames(taxo)[4]
taxMelt$Accuracy[13:16] = taxo[,4]
taxMelt$`Taxonomic Level`[17:20] = colnames(taxo)[5]
taxMelt$Accuracy[17:20] = taxo[,5]
taxMelt$`Taxonomic Level`[21:24] = colnames(taxo)[6]
taxMelt$Accuracy[21:24] = taxo[,6]

taxMelt$`Taxonomic Level` = factor(taxMelt$`Taxonomic Level`, levels = colnames(taxo))
taxMelt$Method = factor(taxMelt$Method, levels = c('glyBERT', 'SweetNet', 'Mohapatra et al.', 'SweetTalk'))

taxMelt = taxMelt[taxMelt$Method != 'Mohapatra et al.', ] # Drop pre=print values for comparison since numerous methods and architectures are used, hard to compare and numbers aren't totally clear...

p5 = ggplot(taxMelt, aes(`Taxonomic Level`, y = Accuracy, fill = Method)) + geom_bar(position = position_dodge(.75), width = .75, stat = 'identity',color = 'grey33', size=.85) + theme_minimal()
p5
dev.off()
p5 + coord_cartesian(ylim = c(0.4,1)) +
  scale_fill_manual(values = c('#16a3f5', '#f77e0d', '#ae16f5')) +
  labs(x = 'Rank of taxonomic origin', y = 'Accuracy') + 
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'bold'),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey")) 
ggsave(filename = 'taxonomyAcc.pdf', device = 'pdf', path = '~/cbk/sugarTrees/plots/',
       width = 11*1.5,
       height = 7.5*1.5)












