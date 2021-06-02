getwd()

library(Biostrings)
library(rphast)
library(ape)

#Loading multiple sequence alignment onto R - uses both Biostrings and rphast
fasta <- readDNAStringSet("Primate_Ensembl_MSA.fa")
seq_names = names(fasta)
sequence = paste(fasta)

primate_msa <- msa(sequence, names = seq_names, alphabet = "ACGT", is.ordered = TRUE, offset = NULL)
primate_msa
summary.msa(primate_msa)


#Change names
names(primate_msa)
primate_msa$names <- c("Human", "Bonobo", "Chimpanzee", "Western_Lowland_Gorilla", "Sumatran_Orangutan", "Northern_White_Cheeked_Gibbon", "Sooty_Mangabey", "Drill", "Olive_Baboon", "Crab_Eating_Macaque", "Rhesus_Macaque", "Pig_Tailed_Macaque", "Vervet", "Black_Snub_Nosed_Monkey", "Golden_Snub_Nosed_Monkey", "Mas_Night_Monkey", "Common_Marmoset", "White_Faced_Capuchin", "Bolivian_Squirrel_Monkey", "Tarsier", "Gray_Mouse_Lemur", "Greater_Bamboo_Lemur", "Coquerels_Sifaka", "Bushbaby")
names(primate_msa)

#Neutral model
neutral_model <- phyloFit(primate_msa, tree = "((((((((Human, (Chimpanzee, Bonobo)), Western_Lowland_Gorilla), Sumatran_Orangutan), Northern_White_Cheeked_Gibbon), (((((Sooty_Mangabey, Drill), Olive_Baboon), ((Crab_Eating_Macaque, Rhesus_Macaque), Pig_Tailed_Macaque)), Vervet), (Black_Snub_Nosed_Monkey, Golden_Snub_Nosed_Monkey))), (Mas_Night_Monkey, (Common_Marmoset, (White_Faced_Capuchin, Bolivian_Squirrel_Monkey)))), Tarsier), ((Gray_Mouse_Lemur, (Greater_Bamboo_Lemur, Coquerels_Sifaka)), Bushbaby))", subst.mod = "SSREV", nrates = 4, EM = TRUE)
neutral_model

human_exclude <- sub.msa(primate_msa, c("Human"), keep = FALSE, refseq = NULL, start.col = 1, end.col = 236)
mod_human_ex <- phyloFit(human_exclude, tree = "((((((((Bonobo, Chimpanzee), Western_Lowland_Gorilla), Sumatran_Orangutan), Northern_White_Cheeked_Gibbon), (((((Sooty_Mangabey, Drill), Olive_Baboon), ((Crab_Eating_Macaque, Rhesus_Macaque), Pig_Tailed_Macaque)), Vervet), (Black_Snub_Nosed_Monkey, Golden_Snub_Nosed_Monkey))), (Mas_Night_Monkey, (Common_Marmoset, (White_Faced_Capuchin, Bolivian_Squirrel_Monkey)))), Tarsier), ((Gray_Mouse_Lemur, (Greater_Bamboo_Lemur, Coquerels_Sifaka)), Bushbaby))", subst.mod = "SSREV", nrates = 4, EM = TRUE)
mod_human_ex

#Non-neutrality in all primates
phyloP_NNEUT <- phyloP(neutral_model, primate_msa, method = "LRT", mode = "NNEUT")
phyloP_NNEUT

write.csv(phyloP_NNEUT, "phyloP_NNEUT.csv")

#Conservation in non-human primates
phyloP_CON <- phyloP(mod_human_ex, human_exclude, method = "LRT", mode = "CON")
phyloP_CON

write.csv(phyloP_CON, "phyloP_CON.csv")

#Acceleration in humans
phyloP_ACC <- phyloP(neutral_model, primate_msa, method = "LRT", mode = "ACC", subtree = "Human")
phyloP_ACC

write.csv(phyloP_ACC, "phyloP_ACC.csv")

#Plotting NNEUT vs CON vs ACC
phyloPtrack_NNEUT <- as.track.wig(coord = phyloP_NNEUT$coord, score = phyloP_NNEUT$score, name = "phyloP score : primate non-neutrality", col = "black", ylim = c(0, 2.5), smooth = TRUE, horiz.line = 0)
phyloPtrack_CON <- as.track.wig(coord = phyloP_CON$coord, score = phyloP_CON$score, name = "phyloP score: non-human conservation", col = "blue", ylim = c(0, 2.5), smooth = TRUE, horiz.line = 0)
phyloPtrack_ACC <- as.track.wig(coord = phyloP_ACC$coord, score = phyloP_ACC$score, name = "phyloP score: human acceleration", col = "red", ylim = c(0, 2.5), smooth = TRUE, horiz.line = 0)

pdf(file = "compare_phyloP.pdf")
plot.track(list(phyloPtrack_NNEUT, phyloPtrack_CON, phyloPtrack_ACC))
dev.off()

#Creating (neutral) tree - uses ape
pdf(file = "neutral_tree.pdf")
par(mfrow = c(1,1))
plot.phylo(read.tree(text = "((((((((Human:0.010032,(Chimpanzee:2.31831e-06,Bonobo:2.19969e-06):0.00502916):2.22361e-06,Western_Lowland_Gorilla:0.0310373):0.0101416,Sumatran_Orangutan:2.15993e-06):2.32086e-06,Northern_White_Cheeked_Gibbon:0.0156079):0.00472851,(((((Sooty_Mangabey:0.0152545,Drill:2.16723e-06):2.20145e-06,Olive_Baboon:0.0100431):0.00487776,((Crab_Eating_Macaque:2.3162e-06,Rhesus_Macaque:0.00498755):2.26277e-06,Pig_Tailed_Macaque:0.0100501):2.17268e-06):2.33237e-06,Vervet:0.0150012):0.0102755,(Black_Snub_Nosed_Monkey:2.28458e-06,Golden_Snub_Nosed_Monkey:0.0051017):0.00494024):0.0210427):0.0116621,(Mas_Night_Monkey:0.0220281,(Common_Marmoset:0.0487036,(White_Faced_Capuchin:0.0102055,Bolivian_Squirrel_Monkey:0.0154411):2.14099e-06):2.21258e-06):0.0364941):0.0647371,Tarsier:0.134945):0.000774658,((Gray_Mouse_Lemur:0.123078,(Greater_Bamboo_Lemur:0.0875855,Coquerels_Sifaka:0.037512):6.20826e-06):0.0418855,Bushbaby:0.0992905):0.000774658);"))
dev.off()