#name: rxodeCommand
#language: r
#input: double inputSD = 1
#output: graphics plot

require("ggplot2")

ds <- data.frame(
  X =  rnorm(200, 0, inputSD)
)

ggplot() +
  geom_density(data = ds, aes(x = X), col = "Red", fill = "Red", alpha = 0.2)+
  theme_bw()
