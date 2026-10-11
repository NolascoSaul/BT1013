# Ejemplo de árbol filogenético con el paquete ape en R
library(ape)

string <- "(Orangutans, (Gorillas, ( (Bonobos, Chimpanzees), Humans) ));"
vert_tree <- read.tree(text = string)
plot(
  vert_tree,
  no.margin = TRUE,
  edge.width = 5,
  main = "Árbol filogenético de los primates"
)

# Visualización del árbol con phytools
library(phytools)

roundPhylogram(
  vert_tree,
  fsize = 1,
  lwd = 2
)

# Árbol sin raíz
plot(
  unroot(vert_tree),
  type = "unrooted",
  no.margin = TRUE,
  lab4ut = "axial",
  edge.width = 2,
  main = "Árbol filogenético sin raíz"
)

# Objeto Filo
tree <- read.tree(text = "(((A,B),(C,D)),E);")
plotTree(tree, offset = 1)
tiplabels()
nodelabels()
str(tree)

# Escribir y leer archivos
tree_text <- "((Anolis_ahli,Anolis_allogus),(Anolis_sagrei,Anolis_bremeri));"
write(tree_text, file = "alo.tre")
anolis_tree <- read.tree(file = "alo.tre")
plot(anolis_tree)
