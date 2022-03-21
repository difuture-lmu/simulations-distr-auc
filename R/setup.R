FILES        = "R/helper.R"
BASE_SEED    = 31415L

# Set different tried out values:
L2SENS       = seq(0.01, 0.09, 0.02)
EPSILON      = seq(0.1, 0.5, 0.1)
DELTA        = seq(0.1, 0.5, 0.1)

# Set TEST variable in `../simulation.R` to control if a test or the full benchmark is run:
if (TEST) {
  REPETITIONS = 10L
} else {
  REPETITIONS  = 10000L
}
