##this script will take the output of the Takara CogentAP dry run output, produce the knee plot and calcualte the barcode distribution inflection.

#input file: demultiplexed_fastqs_counts_all.estimated.csv from the dry run output folder
#ouput: knee_plot_with_inflection.png; calculated inflection point  

# Load necessary libraries
source("/etc/rstudio/Rprofile.site")
library(ggplot2)
library(dplyr)
library(kneedle)
library(plotly)

# Read the CSV file (assuming the first column contains barcode counts)
barcode_counts <- read.csv("./demu_dry/demultiplexed_fastqs/demultiplexed_fastqs_counts_all.estimated.csv", header = F, row.names = 1)
#remove undetermined 
barcode_counts <- barcode_counts %>% filter(V2 != "Non_sample")


# Sort the barcode counts in descending order
sorted_counts <- sort(barcode_counts[, 2], decreasing = TRUE)

# Create a data frame for plotting
plot_data <- data.frame(
  Rank = seq_along(sorted_counts),
  Count = sorted_counts
)

# Calculate the inflection point
# Use the point of maximum curvature on the log-transformed data
log_counts <- log10(plot_data$Count)
log_rank <- log10(plot_data$Rank)

# Calculate the first and second derivatives
first_deriv <- diff(log_counts) / diff(log_rank)
second_deriv <- diff(first_deriv) / diff(log_rank[-1])

# Find the point of maximum curvature (inflection point)
inflection_index <- which.max(abs(second_deriv)) + 1  # +1 to adjust for diff() offset
inflection_rank <- plot_data$Rank[inflection_index]
inflection_count <- plot_data$Count[inflection_index]

# Calculate the knee (inflection) point using Kneedle method

kneedle_result <- kneedle(c(1:length(sorted_counts)),sorted_counts)




# Generate the knee plot with inflection point
knee_plot <- ggplot(plot_data, aes(x = Rank, y = Count)) +
  geom_line(color = "blue") +
  geom_vline(xintercept = inflection_rank, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = inflection_rank, y = max(plot_data$Count), 
           label = paste("Inflection Point:", inflection_rank), 
           hjust = 1.1, vjust = 1.5, color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Knee Plot with Inflection Point",
    x = "Rank (log10)",
    y = "Barcode Count (log10)"
  ) +
  theme_minimal()

# Display the plot
interactive_knee <- ggplotly(knee_plot)

interactive_knee

# save the plot to a file
ggsave("knee_plot_with_inflection.png", knee_plot, width = 6, height = 4)

# Print the inflection point
cat("Inflection Point caculated by deriv:\n")
cat("Rank:", inflection_rank, "\n")
cat("Count:", inflection_count, "\n")
cat("\n")

# Print the inflection point
cat("Inflection Point caculated by kneedle:\n")
cat("rank:", kneedle_result[1], "\n")
cat("count:", kneedle_result[2], "\n")