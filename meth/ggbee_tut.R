library(ggbeeswarm)
library(ggplot2)

data(iris)

#raw
ggplot(iris,aes(x=Species, y=Sepal.Length)) +
  geom_beeswarm()

#grouped
ggplot(iris,aes(x=Species, y=Sepal.Length, colour=Species)) +
  geom_beeswarm()
