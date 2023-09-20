#name: BSA
#description: Body Surface Area
#language: r
#input: double height
#input: double weight
#output: double bsa
#test: BSA(175.145, 77.30) == 1.93926593357154

bsa <- sqrt(weight * height / 3600)
