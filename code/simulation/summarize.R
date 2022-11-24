OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}

devtools::load_all("../pkg")

summarizeSimulation(path = '../../output/simulation/data',
                    dest = '../../output/simulation/summary',
                    workers = 16)
