

uup = read.delim('/Users/dmattox/cbk/sugarTrees/immunoFiles/run2/matchmonoDF.tsv', stringsAsFactors = F)
uup = uup[uup$IUPAC != '',]


uup$isStart = uup$isStart == 'True'

# 
# dropPaths = rep(F, nrow(uup))
# for (i in 1:length(unique(uup$Generative.Path))){
#   p = unique(uup$Generative.Path)[i]
#   hits = unique(uup$Existing.glycan[uup$Generative.Path == p])
#   if (!"Immunogenic" %in% hits){
#     dropPaths[uup$Generative.Path == p] = T
#   }
# }
# 
# uup[dropPaths,]
# 
# uup = uup[!dropPaths, ]

uup[(! uup$isStart) & uup$Existing.glycan == 'Non-immunogenic',]

# uup[uup$Path == 5,]
# uup[uup$Path == 34,]
# uup[uup$Path == 54,]

unique(uup$IUPAC[(uup$isStart)])

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
library(ggalluvial)
library(stringr)
library(ggrepel)


uup$mCnt = str_count(uup$MonoIDs, ',') + 1

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
  x[i] = uup$Prob_imm[uup$IUPAC == nodes[i]][1]
  types = unique(uup$Existing.glycan[uup$IUPAC == nodes[i]])
  if ("Non-immunogenic" %in% types){
    glyType[i] = "Non-immunogenic"
  }else{
    glyType[i] = types[1]
  }
  exists = unique(uup$Existing.glycan[uup$IUPAC == nodes[i]])[1]
  glyExists[i] = exists
}
NodeList <- data.frame(nodes, x, y, glyType, glyExists)




# 
# # NodeList[NodeList$x > 0.89,]
# # NodeList$y[NodeList$x > 0.93] = seq(from = 1, to = nCnt, length.out = sum(NodeList$x > 0.9))
# 
# NodeList$y[NodeList$x > 0.95 & NodeList$y <=7] = seq(from = 1, to = 6, length.out = sum(NodeList$x > 0.95 & NodeList$y <=7))
# 
# NodeList$y[NodeList$x > 0.95 & (8 <= NodeList$y) & (NodeList$y <= 19)] = seq(from = 7, to = 16.6, length.out = sum(NodeList$x > 0.95 & (8 <= NodeList$y) & (NodeList$y <= 19)))
# 
# NodeList$y[NodeList$nodes == 'LDManHep(b1-4)GlcNAc(b1-3)GalNAc'] = 7
# NodeList$y[NodeList$nodes == 'GlcNAc(b1-4)GlcNAc(b1-3)GalNAc'] = 8
# 
# NodeList$y[NodeList$nodes == 'GalOS(b1-4)GlcNAc(b1-3)GalNAc'] = 13.5
# 
# NodeList$y[NodeList$nodes == 'NeuNGc(a2-3)Gal(b1-3)GalNAc'] = 18
# NodeList$y[NodeList$nodes == 'NeuNAc(a2-3)Gal(b1-3)GalNAc'] = 19.5
# NodeList$y[NodeList$nodes == 'GalNAc(a2-3)Gal(b1-3)GalNAc'] = 21
# NodeList$y[NodeList$nodes == 'NeuNAc(a2-3)Gal(b1-3)GlcNAc'] = 17.8




# 
# # 
# # NodeList$y[grepl('LDManHep', NodeList$nodes)] = 14
# # 
# # NodeList[NodeList$nodes %in% c('NeuNAc(a2-3)Gal(b1-3)GalNAc', 'NeuNGc(a2-3)Gal(b1-3)GalNAc', 'GalNAc(a2-3)Gal(b1-3)GalNAc'), ]
# # NodeList[NodeList$nodes %in% c('NeuNAc(a2-3)Gal(b1-3)GalNAc', 'NeuNGc(a2-3)Gal(b1-3)GalNAc', 'GalNAc(a2-3)Gal(b1-3)GalNAc'), ]$y = c(28,31,34)
# # 
# # NodeList[NodeList$x < 0.5,]
# # NodeList[NodeList$x < 0.5,]$y = c(6, 14, 24, 29, 27)
# 
# # NodeList$y[NodeList$glyType == 'Non-immunogenic'] = seq(from = -1, to = -1+(7/nrow(NodeList))*2, length.out = sum(NodeList$glyType == 'Non-immunogenic'))
# # NodeList$y[NodeList$glyType == 'Generated'] = seq(from = -1+(7/nrow(NodeList))*2, to = 1, length.out = sum(NodeList$glyType == 'Generated'))
# 



edgeLst = list()
paths = unique(uup$Path)
for (i in 1:length(paths)){
  p = paths[i]
  pathDat = uup[uup$Path == p,]
  # pathDat = pathDat[nrow(pathDat):1,] # Reverse order to follow steps of path
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
  EdgeList$weight[i] = 2 + 2*sum(edgeLst == e)
}

summary(NodeList$x)
imm2x <- function(immVal,iMin = min(NodeList$x), iMax = max(NodeList$x)){
  # Scale from immunogenicity prob value in NodeList$x to x values scaled from [-1, 1] by igraph
  x = 2*((immVal - iMin)/(iMax-iMin)) - 1
  return(x)
}

summary(NodeList$y)

####################
# Polar coordinates
# Get r,theta of polar points, then convert to x,y
polarNodeList = NodeList

summary(polarNodeList$x)
r = polarNodeList$x # Distance from center is probability of immunogenicity
summary(r)

t = seq(0,2*pi, length.out = nrow(NodeList) + 1)[1:nrow(NodeList)] # Go ALMOST all the way around

polarNodeList$x = r*cos(t)
polarNodeList$y = r*sin(t)

summary(polarNodeList$y)


####################

# a<- graph_from_data_frame(vertices = NodeList, d= EdgeList, directed = T)
a <- graph_from_data_frame(vertices = polarNodeList, d= EdgeList, directed = T)

# V(a)$glyType
# V(a)$glyExists
V(a)$col = rep("", length(V(a)$glyType))
V(a)$col[V(a)$glyType == 'Non-immunogenic'] = 'blue'
V(a)$col[V(a)$glyExists == 'Immunogenic'] = 'red'
V(a)$col[V(a)$glyExists == 'Unknown'] = 'grey'
V(a)$col[V(a)$glyExists == 'Novel'] = 'black'
V(a)$col[V(a)$glyExists == 'Unknown - Imm ano mismatch'] = 'pink'
V(a)$col[V(a)$glyExists == 'Novel - Imm ano mismatch'] = 'orange'


V(a)$cnt = rep(0, nrow(NodeList))
for (i in 1:nrow(NodeList)){
  V(a)$cnt[i] = sum(EdgeList$to == V(a)$name[i])
}



dev.off()

pdf(file = '/Users/dmattox/cbk/sugarTrees/plots/immGenGraph_run2_anomismatch.pdf',
    width = 40,
    height = 40)
plot.igraph(a, 
            edge.width = E(a)$weight,
            edge.arrow.size = .5,
            edge.arrow.width = 1.5,
            edge.curved = T,
            # vertex.label = NA,
            # vertex.label.font = 2,
            vertex.label.cex = 2,
            vertex.label.color = 'mediumpurple1',
            vertex.size = 3 + (V(a)$cnt),
            vertex.label.dist = 0, # 0 - label on vertex, 1 - label above vertex, -1 - label below vertex
            vertex.color = alpha(V(a)$col, 0.6))
# lines(x = imm2x(rep(0.5,2)), y =c(-1.1,1), lty = 2)
# axis(1, at= imm2x(c(0.25,0.5, 0.75, 1)), labels = c(0.25,0.5, 0.75, 1), line = 0, cex.axis = 2)
# title(xlab = 'Predicted probability of immunogenicity', cex.lab = 2)
points(0,0,pch=19, cex = 5)
points(x = c(c(0.25, 0.5, 0.75, 1),-1*c(0.25, 0.5, 0.75, 1)), y=rep(0,8), pch='|', cex = 3)
points(x = rep(0,8), y = c(c(0.25, 0.5, 0.75, 1),-1*c(0.25, 0.5, 0.75, 1)), pch = '-', cex = 4)
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 1000) # angles for drawing points around the circle
rad = c(0.25, 0.5, 0.75, 1)
for (i in 1:length(rad)){
  lines(x = rad[i] * cos(theta) + center_x, y = rad[i] * sin(theta) + center_y,
        lty=3,
        lwd=0.5)
}
text(x = c(0,0.25, 0.5, 0.75, 1), y=rep(0,5),labels = c(0,0.25, 0.5, 0.75, 1), adj = c(0.5,-2),cex=3, font=2)
# text(x = -1*c(0.25, 0.5, 0.75, 1), y=rep(0,5),labels = c(0.25, 0.5, 0.75, 1), adj = c(0.5,-2), cex=2)
dev.off()


sum((NodeList$x > 0.95))
sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Immunogenic')))
sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Unknown')))
sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Unknown - Imm ano mismatch')))
sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Novel')))
sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Novel - Imm ano mismatch')))



# install.packages("waffle", repos = "https://cinc.rud.is")
library(waffle)

# Vector
x_waff <- c(Immunogenic = sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Immunogenic'))),
            Unknown = sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Unknown'))),
            Unknown_AnoMismatch = sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Unknown - Imm ano mismatch'))),
            Novel = sum((NodeList$x > 0.95) & ((NodeList$glyType == 'Novel'))))

# Waffle chart
waffle(x_waff, rows = 10,
       legend_pos = "bottom",
       colors = c("red", "grey", "pink", "black"),
       title = 'Origin of glycans with predicted probablity of immunogenicity > 95%')


communities = cluster_label_prop(a)
length(communities)
communities[2]
communities[8]


#############################################
# Get lists of glycans 


hiImm = unique(uup$IUPAC[(uup$Prob_imm > 0.95) & ((uup$mCnt == 3) | (uup$mCnt == 4 & !grepl('^Fuc',uup$IUPAC)))])
hiCnts = rep(0,length(hiImm))
for (i in 1:length(hiImm)){
  hiCnts[i] = sum(uup$IUPAC == hiImm[i])
}
strsplit(hiImm, '[\\(,\\)]')

immDF = as.data.frame(matrix(NA, nrow = length(hiImm), ncol = 10))
colnames(immDF) = c('Glycan', 'Mono 1', 'Bond 1', 'Mono 2', 'Bond 2', 'Mono 3','Bond 3', 'Mono 4', 'Type', 'Freq')
for (i in 1:length(hiImm)){
  immDF$Glycan[i] = hiImm[i]
  glyBits = strsplit(hiImm[i], '[\\(,\\)]')[[1]]
  if (length(glyBits) == 5){
    glyBits = c(glyBits, c(NA, NA))
  }
  immDF[i,2:8] = glyBits
  immDF$Freq[i] = hiCnts[i]
  immDF$Type[i] = NodeList$glyExists[NodeList$nodes == hiImm[i]]
}

immDF$Type[immDF$Type == 'Unknown - Imm ano mismatch'] = 'Unknown'
immDF$Type[immDF$Type == 'Novel - Imm ano mismatch'] = 'Novel'

# If all trisaccharides:
any(apply(apply(immDF, MARGIN = 2, is.na),
      MARGIN = 2, all))
immDF = immDF[,!apply(apply(immDF, MARGIN = 2, is.na),
                      MARGIN = 2, all)]



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
       aes(y = Freq, axis1 = `Mono 1`, axis2 = `Bond 1`, axis3 = `Mono 2`, axis4 = `Bond 2`, axis5 = `Mono 3`)) + #, axis6 = `Bond 3`, axis7 = `Mono 4`)) +
  geom_alluvium(aes(fill = Type), width = 1/12, ) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Mono 1", "Link 1", "Mono 2", 'Link 2', 'Mono 3'), expand = c(.05, .05)) +#, 'Link 3', 'Mono 4'), expand = c(.05, .05)) +
  # scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_fill_manual(values = alpha(c('red','black', 'grey'), 0.5)) +
  ggtitle("Generated glycans > 95% P(Immunogenicity)") +
  theme_minimal()



# hiImm = uup$IUPAC[(uup$Immuno.Probability > 0.95)]
# hiCnts = rep(0,length(hiImm))
# for (i in 1:length(hiImm)){
#   hiCnts[i] = sum(uup$IUPAC == hiImm[i])
# }
# strsplit(hiImm, '[\\(,\\)]')
# 
# immDF = as.data.frame(matrix(NA, nrow = length(hiImm), ncol = 8))
# colnames(immDF) = c('Glycan', 'Mono 1', 'Bond 1', 'Mono 2', 'Bond 2', 'Mono 3', 'Type', 'Freq')
# for (i in 1:length(hiImm)){
#   immDF$Glycan[i] = i
#   immDF[i,2:6] = strsplit(hiImm[i], '[\\(,\\)]')[[1]]
#   immDF$Freq[i] = hiCnts[i]
#   immDF$Type[i] = NodeList$glyExists[NodeList$nodes == hiImm[i]]
# }

#############

immMelt = melt(data = immDF, id.vars = c('Glycan', 'Type', 'Freq'))

immMelt$value <- factor(immMelt$value, levels = unique(immMelt$value))
is_lodes_form(immMelt, key = "variable", value = 'value', id = 'Glycan')
immMelt$Type[immMelt$Type != 'Immunogenic'] = "Novel/unknown"

ggplot(immMelt,
       aes(x = variable, stratum = value, alluvium = Glycan,label = value, color = Type, fill = value)) +
  scale_color_manual(values = alpha(c('red','black'), 1)) +
  scale_fill_manual(values = alpha(c("#9e1fff", # neuac/Purple"
                                     "#ffd900", # gal/yellow
                                     "#ffd900", # gal/yellow
                                     "#96f2f7", # NeuGc, light blue
                                     "#0080ff", # glc/blue
                                     "#0080ff", # glc/blue
                                     "#00A651", #Rha/green
                                     "#ffd900", # gal/yellow
                                     "#ff0000", # fuc/red
                                     'grey66',
                                     'grey66',
                                     'grey66',
                                     "#00A651", #Rha/green
                                     "#00A651", #Rha/green
                                     "#0080ff", # glc/blue
                                     "white", # b-tri-ol
                                     "#ffd900", # gal/yellow
                                     "#0080ff", # glc/blue
                                     "#0080ff", # glc/blue
                                     "#00A651", #Rha/green
                                     "#ffd900", # gal/yellow
                                     "#ffd900", # gal/yellow
                                     "#ffd900", # gal/yellow
                                     "#0080ff"# glc/blue
                                     # "#ffd900", # gal/yellow
                                     ), 0.8)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",size=0.75) +
  geom_stratum(color = 'black') +
  geom_label_repel(stat = "stratum" , color = 'black', point.size = NA)+#, aes(label = after_stat(stratum))) +
  theme(legend.position = "bottom") +
  ggtitle("Generated glycans > 95% P(Immunogenicity)") +
  theme_minimal()
ggsave('~/cbk/sugarTrees/plots/immGen2_alluvial_all.pdf', device = 'pdf', width = 15, height = 17, units = 'in', dpi = 600)



