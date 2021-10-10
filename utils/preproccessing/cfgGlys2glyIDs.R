

rfuDat <- read.delim('~/cbk/sugarTrees/glyPredict/Clean_Lectin_intensities_67uni.txt', sep = '\t', stringsAsFactors = F, strip.white = T)
rfuDat$glyName = gsub('[[:blank:]]*', '', rfuDat$glyName)


emb <- read.delim('~/cbk/sugarTrees/glyPredict/CFG_glycans_535_feb9_emb.csv', sep = ',' )
row.names(emb) = emb$glyID
emb$glyID <- NULL



tests = read.delim('~/cbk/sugarTrees/glyPredict/coff_et_al_test-train_splits/single_test_train_split/training_set/training_set_ABA_14361_100ug_v5.0_DATA.csv', stringsAsFactors = F,sep = ',')

sum(tests$glycan %in% rfuDat$glyName)
nrow(tests)

tests$glycan[! tests$glycan %in% rfuDat$glyName]

# "GlcNAcb1-4-MDPLys" --> "GlcNAcb1"
rfuDat$glyName[grepl('^GlcNAcb', rfuDat$glyName)]

# "Galb1-4GlcNAcb1-6(Galb1-4GlcNAcb1-2)Mana1-6(GlcNAcb1-4Galb1-4GlcNAcb1-4(Galb1-4GlcNAcb1-2)Mana1-3)Manb1-4GlcNAcb1-4(Fuca1-6)GlcNAc-Sp21"
rfuDat$glyName[rfuDat$glyID == 509]
# improperly formatted IUPAC string for probe 509, nit included with glycan embeddings


# "Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-6(Neu5Aca2-6Galb1-4GlcNAcb1-2Mana1-3)Manb1-4GlcNAcb1-4GlcNAcb-Sp21"
rfuDat$glyName[rfuDat$glyID == 57] # Extra hyphen in CFG iupac string

rfuDat$glyName[grepl('^Gala1-3.*-Sp14', rfuDat$glyName)]

