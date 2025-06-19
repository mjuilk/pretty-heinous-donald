# Load the necessary library
library(ggplot2)

# Sample data
set.seed(123)
data <- data.frame(
  x = rnorm(1000),
  y = rnorm(1000),
  z = rnorm(1000)
)

# Plot using stat_summary_hex
ggplot(data, aes(x = x, y = y, z = z)) +
  stat_summary_hex(aes(fill = after_stat(value)), bins = 30, fun = mean) +
  scale_fill_gradientn(colors = c("blue", "yellow", "red")) + # Customize the color gradient
  theme_minimal() +
  labs(
    title = "Hexagonal Binning with Z Value Color Gradient",
    x = "X-axis",
    y = "Y-axis",
    fill = "Mean Z"
  )
