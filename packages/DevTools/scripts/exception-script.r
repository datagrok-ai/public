#name: ExceptionScriptR
#language: r
#input: int a
#output: int out

if (a == 0) {
  stop("exception")
} else {
  a <- a + 1
}
out <- a
