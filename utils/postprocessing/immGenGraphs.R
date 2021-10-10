

uup = read.delim('/Users/dmattox/cbk/sugarTrees/immGenDF.tsv', stringsAsFactors = F)
uup = uup[uup$IUPAC != '',]

dropPaths = rep(F, nrow(uup))
for (i in 1:length(unique(uup$Generative.Path))){
  p = unique(uup$Generative.Path)[i]
  hits = unique(uup$Existing.glycan[uup$Generative.Path == p])
  if (!"Immunogenic" %in% hits){
    dropPaths[uup$Generative.Path == p] = T
  }
}

uup[dropPaths,]

uup = uup[!dropPaths, ]

######################################################
# library("Rgraphviz")
# V <- letters[1:10]
# M <- 1:4
# g1 <- randomGraph(V, M, 0.2)
# g1 <- layoutGraph(g1)
# 
# graph.par(list(nodes=list(fill="lightgray", textCol="red")))
# graph.par(list(nodes=list(col="darkgreen", lty="dotted", lwd=2, fontsize=6)))
# graph.par(list(edges=list(col="lightblue", lty="dashed", lwd=3)))
# renderGraph(g1)
# 
# 
# nodes = unique(uup$IUPAC)
# nodes = rbind(nodes, c(LETTERS, letters[1:(length(nodes) - length(LETTERS))]))

######################################################
library("igraph")
library("scales")
library(ggplot2)
library(reshape2)
graph_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                     c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph <- graph.edgelist(graph_edges, directed = TRUE)
plot.igraph(graph, axes =F)
axis(1, seq(-1, 1, 0.5))


nodes = unique(uup$IUPAC)
nCnt = length(nodes)
y = 1:nCnt
x = rep(0, nCnt)
glyType = rep('', nCnt)
glyExists = rep('', nCnt)
for (i in 1:nCnt){
  x[i] = uup$Immuno.Probability[uup$IUPAC == nodes[i]][1]
  types = unique(uup$Glycan.Type[uup$IUPAC == nodes[i]])
  if ("Non-immunogenic" %in% types){
    glyType[i] = "Non-immunogenic"
  }else{
    glyType[i] = types[1]
  }
  exists = unique(uup$Existing.glycan[uup$IUPAC == nodes[i]])[1]
  glyExists[i] = exists
}
NodeList <- data.frame(nodes, x, y, glyType, glyExists)

# NodeList[NodeList$x > 0.89,]
# NodeList$y[NodeList$x > 0.93] = seq(from = 1, to = nCnt, length.out = sum(NodeList$x > 0.9))

NodeList$y[NodeList$x > 0.95 & NodeList$y <=7] = seq(from = 1, to = 6, length.out = sum(NodeList$x > 0.95 & NodeList$y <=7))

NodeList$y[NodeList$x > 0.95 & (8 <= NodeList$y) & (NodeList$y <= 19)] = seq(from = 7, to = 16.6, length.out = sum(NodeList$x > 0.95 & (8 <= NodeList$y) & (NodeList$y <= 19)))

NodeList$y[NodeList$nodes == 'LDManHep(b1-4)GlcNAc(b1-3)GalNAc'] = 7
NodeList$y[NodeList$nodes == 'GlcNAc(b1-4)GlcNAc(b1-3)GalNAc'] = 8

NodeList$y[NodeList$nodes == 'GalOS(b1-4)GlcNAc(b1-3)GalNAc'] = 13.5

NodeList$y[NodeList$nodes == 'NeuNGc(a2-3)Gal(b1-3)GalNAc'] = 18
NodeList$y[NodeList$nodes == 'NeuNAc(a2-3)Gal(b1-3)GalNAc'] = 19.5
NodeList$y[NodeList$nodes == 'GalNAc(a2-3)Gal(b1-3)GalNAc'] = 21
NodeList$y[NodeList$nodes == 'NeuNAc(a2-3)Gal(b1-3)GlcNAc'] = 17.8





# 
# NodeList$y[grepl('LDManHep', NodeList$nodes)] = 14
# 
# NodeList[NodeList$nodes %in% c('NeuNAc(a2-3)Gal(b1-3)GalNAc', 'NeuNGc(a2-3)Gal(b1-3)GalNAc', 'GalNAc(a2-3)Gal(b1-3)GalNAc'), ]
# NodeList[NodeList$nodes %in% c('NeuNAc(a2-3)Gal(b1-3)GalNAc', 'NeuNGc(a2-3)Gal(b1-3)GalNAc', 'GalNAc(a2-3)Gal(b1-3)GalNAc'), ]$y = c(28,31,34)
# 
# NodeList[NodeList$x < 0.5,]
# NodeList[NodeList$x < 0.5,]$y = c(6, 14, 24, 29, 27)

# NodeList$y[NodeList$glyType == 'Non-immunogenic'] = seq(from = -1, to = -1+(7/nrow(NodeList))*2, length.out = sum(NodeList$glyType == 'Non-immunogenic'))
# NodeList$y[NodeList$glyType == 'Generated'] = seq(from = -1+(7/nrow(NodeList))*2, to = 1, length.out = sum(NodeList$glyType == 'Generated'))




edgeLst = list()
paths = unique(uup$Generative.Path)
for (i in 1:length(paths)){
  p = paths[i]
  pathDat = uup[uup$Generative.Path == p,]
  pathDat = pathDat[nrow(pathDat):1,] # Reverse order to follow steps of path
  for (j in 2:nrow(pathDat)){
    e = paste(pathDat$IUPAC[j-1], pathDat$IUPAC[j], sep = '_')
    edgeLst[[(length(edgeLst) + 1)]] = e
  }
}
edgeLst = unlist(edgeLst)
EdgeList = as.data.frame(matrix(0, nrow = length(unique(edgeLst)), ncol =3))
colnames(EdgeList) = c('from', 'to', 'weight')
for (i in 1:nrow(EdgeList)){
  e = unique(edgeLst)[i]
  EdgeList$from[i] = strsplit(e, '_')[[1]][1]
  EdgeList$to[i] = strsplit(e, '_')[[1]][2]
  EdgeList$weight[i] = 1 + sum(edgeLst == e)
}

summary(NodeList$x)
imm2x <- function(immVal,iMin = min(NodeList$x), iMax = max(NodeList$x)){
  # Scale from immunogenicity prob value in NodeList$x to x values scaled from [-1, 1] by igraph
  x = 2*((immVal - iMin)/(iMax-iMin)) - 1
  return(x)
}


a<- graph_from_data_frame(vertices = NodeList, d= EdgeList, directed = T)

# V(a)$glyType
# V(a)$glyExists
V(a)$col = rep("", length(V(a)$glyType))
V(a)$col[V(a)$glyType == 'Non-immunogenic'] = 'blue'
V(a)$col[V(a)$glyExists == 'Immunogenic'] = 'red'
V(a)$col[V(a)$glyExists == 'Unknown'] = 'grey'
V(a)$col[V(a)$glyExists == 'Novel'] = 'black'

V(a)$cnt = rep(0, nrow(NodeList))
for (i in 1:nrow(NodeList)){
  V(a)$cnt[i] = sum(EdgeList$to == V(a)$name[i])
}


pdf(file = '/Users/dmattox/cbk/sugarTrees/plots/immGenGraph.pdf',
    width = 17,
    height = 17)
plot.igraph(a, 
            edge.width = E(a)$weight,
            edge.arrow.size = .5,
            edge.arrow.width = 1.5,
            edge.curved = T,
            vertex.label = NA,
            vertex.size = 3 + 2*(V(a)$cnt),
            vertex.label.dist = 0, # 0 - label on vertex, 1 - label above vertex, -1 - label below vertex
            vertex.color = V(a)$col)
# lines(x = imm2x(rep(0.5,2)), y =c(-1.1,1), lty = 2)
axis(1, at= imm2x(c(0.25,0.5, 0.75, 1)), labels = c(0.25,0.5, 0.75, 1), line = 0, cex.axis = 2)
title(xlab = 'Predicted probability of immunogenicity', cex.lab = 2)
dev.off()




###############
# Get lists of high prob of immunogenicity glycans

hiImm = unique(uup$IUPAC[(uup$Immuno.Probability > 0.95)])
hiCnts = rep(0,length(hiImm))
for (i in 1:length(hiImm)){
  hiCnts[i] = sum(uup$IUPAC == hiImm[i])
}
strsplit(hiImm, '[\\(,\\)]')

immDF = as.data.frame(matrix(NA, nrow = length(hiImm), ncol = 8))
colnames(immDF) = c('Glycan', 'Mono 1', 'Bond 1', 'Mono 2', 'Bond 2', 'Mono 3', 'Type', 'Freq')
for (i in 1:length(hiImm)){
  immDF$Glycan[i] = hiImm[i]
  immDF[i,2:6] = strsplit(hiImm[i], '[\\(,\\)]')[[1]]
  immDF$Freq[i] = hiCnts[i]
  immDF$Type[i] = NodeList$glyExists[NodeList$nodes == hiImm[i]]
}

immMelt = melt(data = immDF, id.vars = c('Glycan', 'Type', 'Freq'))



# 
# is_alluvia_form(as.data.frame(immDF), axes = 1:5, silent = TRUE)
# 
# immDF$M1 = factor(immDF$M1, levels = unique(immDF$M1))
# immDF$B1 = factor(immDF$B1, levels = unique(immDF$B1))
# immDF$M2 = factor(immDF$M2, levels = unique(immDF$M2))
# immDF$B2 = factor(immDF$B2, levels = unique(immDF$B2))
# immDF$M3 = factor(immDF$M3, levels = unique(immDF$M3))


ggplot(as.data.frame(immDF),
       aes(y = Freq, axis1 = `Mono 1`, axis2 = `Bond 1`, axis3 = `Mono 2`, axis4 = `Bond 2`, axis5 = `Mono 3`)) +
  geom_alluvium(aes(fill = Type), width = 1/12, ) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Mono 1", "Link 1", "Mono 2", 'Link 2', 'Mono 3'), expand = c(.05, .05)) +
  # scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_fill_manual(values = alpha(c('red','grey', 'black'), 0.5)) +
  ggtitle("Generated glycans > 95% P(Immunogenicity)") +
  theme_minimal()




hiImm = uup$IUPAC[(uup$Immuno.Probability > 0.95)]
hiCnts = rep(0,length(hiImm))
for (i in 1:length(hiImm)){
  hiCnts[i] = sum(uup$IUPAC == hiImm[i])
}
strsplit(hiImm, '[\\(,\\)]')

immDF = as.data.frame(matrix(NA, nrow = length(hiImm), ncol = 8))
colnames(immDF) = c('Glycan', 'Mono 1', 'Bond 1', 'Mono 2', 'Bond 2', 'Mono 3', 'Type', 'Freq')
for (i in 1:length(hiImm)){
  immDF$Glycan[i] = i
  immDF[i,2:6] = strsplit(hiImm[i], '[\\(,\\)]')[[1]]
  immDF$Freq[i] = hiCnts[i]
  immDF$Type[i] = NodeList$glyExists[NodeList$nodes == hiImm[i]]
}



immMelt = melt(data = immDF, id.vars = c('Glycan', 'Type', 'Freq'))

immMelt$value <- factor(immMelt$value, levels = unique(immMelt$value))
is_lodes_form(immMelt, key = "variable", value = 'value', id = 'Glycan')
immMelt$Type[immMelt$Type != 'Immunogenic'] = "Novel/unknown"

ggplot(immMelt,
       aes(x = variable, stratum = value, alluvium = Glycan,label = value, color = Type, fill = value)) +
  scale_color_manual(values = alpha(c('red','black'), 1)) +
  scale_fill_manual(values = alpha(c("#ffd900", # gal/yellow
                                     "#0080ff", # glc/blue
                                     "#ffd900",
                                     "#ffd900",
                                     "#9e1fff",
                                     'grey30',
                                     'grey30',
                                     "#0080ff",
                                     'grey30',
                                     "#0080ff",
                                     "#ff0000"), 0.8)) +  # fuc/red
  geom_flow(stat = "alluvium", lode.guidance = "frontback",size=0.75) +
  geom_stratum(width = c(rep(1/3,5),
                         1/6,1/6,
                         rep(1/3,3),
                         1/6,1/6,
                         rep(1/3,6)), color = 'black') +
  geom_label(stat = "stratum" , fill = 'ivory2', color = 'black')+#, aes(label = after_stat(stratum))) +
  theme(legend.position = "bottom") +
  ggtitle("Generated glycans > 95% P(Immunogenicity)") +
  theme_minimal()
ggsave('~/cbk/sugarTrees/plots/immGen1_alluvial_all.pdf', device = 'pdf', width = 9, height = 12, units = 'in', dpi = 300)

####################
# Break it up by starting glycan structures
immMelt1 = melt(data = immDF[immDF$`Bond 1` == 'b1-4' & immDF$`Bond 2` == 'b1-4',], id.vars = c('Glycan', 'Type', 'Freq'), measure.vars = c('Mono 1', 'Mono 2', 'Mono 3'))
immMelt2 = melt(data = immDF[immDF$`Bond 1` == 'b1-4' & immDF$`Bond 2` == 'b1-3',], id.vars = c('Glycan', 'Type', 'Freq'), measure.vars = c('Mono 1', 'Mono 2', 'Mono 3'))
immMelt3 = melt(data = immDF[immDF$`Bond 1` == 'a2-3' & immDF$`Bond 2` == 'b1-3',], id.vars = c('Glycan', 'Type', 'Freq'), measure.vars = c('Mono 1', 'Mono 2', 'Mono 3'))

immMelt1$value <- factor(immMelt1$value, levels = unique(immMelt1$value))
immMelt1$Type[immMelt1$Type != 'Immunogenic'] = "Novel/unknown"
immMelt2$value <- factor(immMelt2$value, levels = unique(immMelt2$value))
immMelt2$Type[immMelt2$Type != 'Immunogenic'] = "Novel/unknown"
immMelt3$value <- factor(immMelt3$value, levels = unique(immMelt3$value))
immMelt3$Type[immMelt3$Type != 'Immunogenic'] = "Novel/unknown"

ggplot(immMelt1,
       aes(x = variable, stratum = value, alluvium = Glycan,label = value, color = Type, fill = value)) +
  scale_color_manual(values = alpha(c('red','black'), 1)) +
  scale_fill_manual(values = alpha(c("#ffd900", # gal/yellow
                                     "#0080ff", # glc/blue
                                     "#ffd900"), 0.8)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",size=0.75) +
  geom_stratum() +
  geom_label(stat = "stratum" , fill = 'ivory2', color = 'black')+#, aes(label = after_stat(stratum))) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,16)) +
  ggtitle("Generated glycans > 95% P(Immunogenicity)") +
  theme_minimal()
ggsave('~/cbk/sugarTrees/plots/immGen1_alluvial_1.pdf', device = 'pdf', width = 9, height = 12, units = 'in', dpi = 300)


ggplot(immMelt2,
       aes(x = variable, stratum = value, alluvium = Glycan,label = value, color = Type, fill = value)) +
  scale_color_manual(values = alpha(c('red','black'), 1)) +
  scale_fill_manual(values = alpha(c("#ffd900", # gal/yellow
                                     "#ffd900",
                                     "#0080ff", # glc/blue
                                     "#0080ff", # glc/blue
                                     "#0080ff", # glc/blue
                                     "#ff0000"), 0.8)) +  # fuc/red
  geom_flow(stat = "alluvium", lode.guidance = "frontback",size=0.75) +
  geom_stratum() +
  geom_label(stat = "stratum" , fill = 'ivory2', color = 'black')+#, aes(label = after_stat(stratum))) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,16)) +
  ggtitle("Generated glycans > 95% P(Immunogenicity)") +
  theme_minimal()
ggsave('~/cbk/sugarTrees/plots/immGen1_alluvial_2.pdf', device = 'pdf', width = 9, height = 12, units = 'in', dpi = 300)


ggplot(immMelt3,
       aes(x = variable, stratum = value, alluvium = Glycan,label = value, color = Type, fill = value)) +
  scale_color_manual(values = alpha(c('red','black'), 1)) +
  scale_fill_manual(values = alpha(c("#9e1fff",
                                     "#ffd900", # gal/yellow
                                     "#0080ff"), 0.8)) +  # # glc/blue
  geom_flow(stat = "alluvium", lode.guidance = "frontback",size=0.75) +
  geom_stratum() +
  geom_label(stat = "stratum" , fill = 'ivory2', color = 'black')+#, aes(label = after_stat(stratum))) +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0,16)) +
  ggtitle("Generated glycans > 95% P(Immunogenicity)") +
  theme_minimal()
ggsave('~/cbk/sugarTrees/plots/immGen1_alluvial_3.pdf', device = 'pdf', width = 9, height = 12, units = 'in', dpi = 300)
