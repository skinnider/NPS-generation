# Consolidate unique masses across all CV folds.
setwd("~/git/deepmet")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-write-uniq-CV.R')
parser$add_argument('--allocation', type = 'character', default = "root")
args = parser$parse_args()

library(tidyverse)
library(magrittr)
library(maslib)
detect_system()

# load grid functions
source("R/functions/grids.R")

# set up grid
grid = data.frame(
  database = 'HMDB4',
  conditional = FALSE,
  pretrain = c('ChEMBL', 'PubChem', 'GDB17'),
  enum_factor = 30,
  downsample_nonlipids = c(TRUE, FALSE, TRUE),
  downsample_n = c(1000, 250, 1000),
  include_lipids = c(FALSE, TRUE, TRUE),
  k = 10,
  cv_fold = '$', 
  
  # architecture - fixed
  rnn_type = 'LSTM',
  embedding_size = 128,
  hidden_size = 1024,
  n_layers = 3,
  dropout = 0,
  batch_size = 64,
  learning_rate = 0.001,
  
  # 1b molecules at a time  
  sample_idx = 1
) %>% 
  tidyr::crossing(
    # different ways of summarizing across folds
    summary_fn = c('freq_sum', 'freq_avg', 'fp10k')
  )

# define input files
input_dirnames = pmap_chr(grid, function(...) {
  current = tibble(...)
  current %>%
    dplyr::select(database:learning_rate) %>%
    map2(., names(.), ~ {
      paste0(.y, '=', .x)
    }) %>%
    paste0(collapse = '-')
})
input_dirs = file.path(base_dir, "CV", "sample", input_dirnames)
input_files = with(grid, file.path(input_dirs, 
                                   paste0('unique-masses-novel-', 
                                          sample_idx, '.tsv')))

# define output files
output_files = map2_chr(input_files, grid$summary_fn, ~ gsub('\\$', .y, .x)) %>% 
  gsub("\\.gz$", "", .)

# define CV files (train data)
cv_dir = file.path(base_dir, "CV", "inputs")
cv_filenames = pmap_chr(grid, function(...) {
  current = tibble(...)
  current %>%
    ## replace enum factor with 0
    mutate(enum_factor = 0) %>% 
    dplyr::select(database, enum_factor, downsample_nonlipids, downsample_n, 
                  include_lipids, k, cv_fold) %>%
    map2(., names(.), ~ {
      paste0(.y, '=', .x)
    }) %>%
    paste0(collapse = '-') %>%
    paste0('train-', ., '.csv')
})
cv_files = file.path(cv_dir, cv_filenames)

# now, check for which parameters are already complete
grid0 = grid %>% 
  mutate(input_file = input_files,
         cv_file = cv_files,
         output_file = output_files) %>% 
  sample_frac(1) %>% 
  mutate(keep = map_lgl(input_file, ~ {
    file = .x
    files = map_chr(seq_len(10), ~ gsub('\\$', as.character(.x), file))
    all(file.exists(files))
  })) %>% 
  filter(keep, !file.exists(output_file %>% paste0('.gz')))

# write the grid that still needs to be run
grid_file = "sh/grids/CV/write-uniq-CV.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = TRUE)
write.table(grid0, grid_file, quote = FALSE, row.names = FALSE, sep = "\t")

# write the sh file dynamically
sh_file = '~/git/deepmet/sh/CV/write-uniq-CV.sh'
write_sh(job_name = 'write-uniq-CV',
         sh_file = sh_file,
         grid_file = grid_file,
         inner_file = 'R/sample-molecules/inner-write-uniq-CV.R',
         time = 12,
         mem = 220,
         env = '/scratch/st-ljfoster-1/decoy-generation/env-r')

# finally, run the job on whatever system we're on
submit_job(grid0, sh_file, allocation = allocation)
