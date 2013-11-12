trajectories <- setClass(
  "trajectories",
  representation(
    time = "numeric",
    sigma = "matrix",
    G = "matrix"
  )
)