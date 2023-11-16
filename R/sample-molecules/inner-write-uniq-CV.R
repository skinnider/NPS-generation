setwd("~/git/deepmet")
options(stringsAsFactors = FALSE)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-write-uniq-CV.R')
grid = read.delim("sh/grids/CV/write-uniq-CV.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(data.table)

# set up input files
cv_folds = seq_len(10)
input_files = map_chr(cv_folds, ~ gsub('\\$', as.character(.x), 
                                       args$input_file))

# read them all
dats = map(input_files, fread)

# merge all of the sampled datasets across CV folds
uniq_smiles = map(dats, 'smiles') %>% 
  Reduce(union, .)
mat = matrix(0, nrow = length(uniq_smiles), ncol = length(cv_folds),
             dimnames = list(uniq_smiles, as.character(cv_folds)))
for (fold_idx in seq_along(dats)) {
  dat = dats[[fold_idx]]
  mat[dat$smiles, fold_idx] = dat$size
}
# fill NAs with zeroes (freq=0)
mat[is.na(mat)] = 0

# read the training set SMILES, too
for (fold_idx in seq_along(dats)) {
  cv_file = gsub('\\$', as.character(fold_idx), args$cv_file)
  cv_dat = read.csv(cv_file) %>% 
    # limit to SMILES already in the matrix (no point adding rows with just NA)
    filter(smiles %in% uniq_smiles)
  # censor these values
  mat[cv_dat$smiles, as.character(fold_idx)] = NA
}

# optionally, normalize by total sampling frequency
if (args$summary_fn == 'fp10k') {
  mat = 10e3 * mat / colSums(mat, na.rm = TRUE)
}

# calculate mean / sum
if (args$summary_fn == 'freq_sum') {
  sums = rowSums(mat, na.rm = TRUE)
  df = data.frame(smiles = rownames(mat), size = sums)
} else {
  means = rowMeans(mat, na.rm = TRUE)
  df = data.frame(smiles = rownames(mat), size = means)
}
# arrange by size
df %<>%
  arrange(desc(size)) %>% 
  filter(size > 0)

# add in metadata (mass and formula)
meta = dats %>% 
  map(~ dplyr::select(.x, 'smiles', 'mass', 'formula')) %>% 
  bind_rows() %>% 
  distinct()
## keep only one row per SMILES
dup = which(duplicated(meta$smiles))
meta %<>% extract(-dup, )
# merge into the frequency data frame
df %<>% left_join(meta, by = 'smiles') %>% 
  dplyr::select(smiles, mass, formula, size)

# write
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
fwrite(df, args$output_file)

# gzip
system(paste('gzip --force', args$output_file))
