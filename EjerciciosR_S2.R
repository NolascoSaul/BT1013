dna_bases <- c("A", "C", "G", "T")

# Secuencia aleatoria de DNA de tamaño n
rand_dna <- function(n) {
  paste(sample(dna_bases, n, replace = TRUE), collapse = "")
}

# Tamaño de una secuencia de DNA
size_dna <- function(dna) {
  nchar(dna)
}

# Porcentaje de cada base en una secuencia de DNA
base_percentage <- function(dna) {
  total_bases <- size_dna(dna)
  percentages <- sapply(dna_bases, function(base) {
    count <- sum(strsplit(dna, "")[[1]] == base)
    (count / total_bases) * 100
  })
  names(percentages) <- dna_bases
  percentages
}

# Transcribe de DNA a RNA
dna2rna <- function(dna) {
  rna <- chartr("T", "U", dna)
  rna
}

# Traducción de RNA a proteínas
rna2protein <- function(rna) {
  codon_table <- c(
    # Para algunos casos, la "N" cae en la tercera posición de un codón
    # y no importa qué letra sea, el resultado es el mismo
    "UUU" = "F", "UUC" = "F", "UUA" = "L", "UUG" = "L",
    "CUU" = "L", "CUC" = "L", "CUA" = "L", "CUG" = "L", "CUN" = "L",
    "AUU" = "I", "AUC" = "I", "AUA" = "I", "AUG" = "M",
    "GUU" = "V", "GUC" = "V", "GUA" = "V", "GUG" = "V", "GUN" = "V",
    "UCU" = "S", "UCC" = "S", "UCA" = "S", "UCG" = "S", "UCN" = "S",
    "CCU" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P", "CCN" = "P",
    "ACU" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T", "ACN" = "T",
    "GCU" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A", "GCN" = "A",
    "UAU" = "Y", "UAC" = "Y", "UAA" = "STOP", "UAG" = "STOP",
    "CAU" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "AAU" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "GAU" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "UGU" = "C", "UGC" = "C", "UGA" = "STOP", "UGG" = "W",
    "CGU" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R", "CGN" = "R",
    "AGU" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GGU" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G", "GGN" = "G"
  )

  protein <- c()
  for (i in seq(1, nchar(rna) - 2, by = 3)) {
    codon <- substr(rna, i, i + 2)

    if (grepl("-", codon)) next

    aminoacid <-
      if (codon %in% names(codon_table)) codon_table[[codon]] else "X"

    if (aminoacid == "STOP") break

    protein <- c(protein, aminoacid)
  }

  paste(protein, collapse = "")
}

# Hebra inversa
reverse_strand <- function(forward_strand) {
  rev_strand <- paste(rev(strsplit(forward_strand, "")[[1]]), collapse = "")
  rev_strand
}

# Hebra complementaria
complementary_strand <- function(forward_strand) {
  complement <- chartr("ACGT", "TGCA", forward_strand)
  complement
}

# Hebra complementaria inversa
rev_complementary_strand <- function(complementary_strand) {
  rev_comp <- reverse_strand(complementary_strand)
  rev_comp
}

# Ejemplo de uso
dna <- rand_dna(15)
rna <- dna2rna(dna)
comp_dna <- complementary_strand(dna)

cat("DNA original:                ", dna, "\n")
cat("Tamaño:                      ", size_dna(dna), "\n")
cat("Porcentaje de bases:\n")
print(base_percentage(dna))
cat("RNA transcrito:              ", rna, "\n")
cat("Proteína traducida:          ", rna2protein(rna), "\n")
cat("Hebra inversa:               ", reverse_strand(dna), "\n")
cat("Hebra complementaria:        ", comp_dna, "\n")
cat("Hebra complementaria inversa:", rev_complementary_strand(comp_dna), "\n")
