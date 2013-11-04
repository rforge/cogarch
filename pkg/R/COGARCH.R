COGARCH <- setClass(
  "COGARCH",
  representation(
    time = "numeric",
    sigma = "matrix",
    G = "matrix"
  )
)