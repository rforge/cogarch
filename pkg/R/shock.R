shock <- setClass(
  "shock",
  representation(
    timeinc = "numeric",
    inc = "matrix",
    Time = "numeric"
  )
)