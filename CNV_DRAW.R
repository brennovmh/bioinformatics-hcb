library(GenVisR)

chromosome <- "chr3"

# limites aproximados do chr3 em hg38
chr_end <- 198295559

# deleção real
del_start <- 197676964
del_end   <- 197686472

set.seed(1)

# aumenta bastante a densidade ao longo de todo o cromossomo
n_points <- 8000
coordinate <- sort(sample(1:chr_end, size = n_points, replace = FALSE))

# CN ~2 no cromossomo inteiro e CN ~1 na deleção
cn <- ifelse(coordinate >= del_start & coordinate <= del_end,
             rnorm(length(coordinate), mean = 1.0, sd = 0.06),
             rnorm(length(coordinate), mean = 2.0, sd = 0.06))

data <- data.frame(
  chromosome = chromosome,
  coordinate = coordinate,
  cn = cn,
  stringsAsFactors = FALSE
)

# segmento marcando explicitamente a deleção
dataSeg <- data.frame(
  chromosome = "chr3",
  start      = del_start,
  end        = del_end,
  segmean    = 1,
  stringsAsFactors = FALSE
)

cnView(
  data,
  z = dataSeg,
  chr = "chr3",
  genome = "hg38",
  ideogram_txtSize = 3
)
