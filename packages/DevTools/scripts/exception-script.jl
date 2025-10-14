#name: ExceptionScriptJulia
#language: julia
#input: int a
#output: int out

if a == 0
    error("exception")
else
    a += 1
end
out = a
