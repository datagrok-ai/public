#name: BSA
#description: Body Surface Area
#language: r
#input: double height
#input: double weight
#output: double bsa

bsa <- sqrt(weight * height / 3600)
