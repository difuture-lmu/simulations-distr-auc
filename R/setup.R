FILES        = "R/helper.R"
BASE_SEED    = 31415L

# Set different tried out values:
L2SENS       = c(0.01, 0.1, 0.2, 0.3, 0.4)
EPSILON      = c(0.1, 0.5, 1, 5, 10)

DELTA        = c(0.00001, 0.0001, 0.001, 0.01, 0.1)



# Set TEST variable in `../simulation.R` to control if a test or the full benchmark is run:
if (TEST) {
  REPETITIONS = 10L
} else {
  REPETITIONS  = 10000L
}
