
#=====================================================
# preamble
#=====================================================

if (!require("pacman")) install.packages("pacman")

library("phylostratr")
library("data.table")




#=====================================================
# preamble
#=====================================================

f = "D:/programming/iodine/output/ABcTax.csv"

dat = fread(
  f
)

library("ape")


tree = ncbi_tree(dat$taxId)

treeDist = cophenetic.phylo(tree)

treeDist = data.table(treeDist)

treeDist$taxId = as.numeric(names(treeDist))

names(treeDist) = paste0("taxid_", names(treeDist))



fwrite(treeDist, "D:/programming/iodine/output/taxIdDistance.csv")




