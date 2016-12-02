setwd("~/Developer/Projects/ITKPhase/src")

library(ggplot2)
library(reshape2)

# inputFile = '../data/harp/congruence.csv'
# outputFile = '../data/harp/congruence.png'
inputFile = '../data/swi/congruence.csv'
outputFile = '../data/swi/congruence.png'

data <- read.csv(inputFile)

# # Remove indices that were cropped from HARP image 
# data <- data[!(data$Index < 34),]
# data <- data[!(data$Index > 150),]

data <- melt(data, id.vars = 'Index')

names(data)[2] <- 'Group'

theme_set(theme_gray(base_size = 20))

ggplot(data, aes(x = Index, y = value)) +
  geom_line(aes(color = Group)) +
  labs(x = 'Index', y = 'Pixel Intensity')

ggsave(outputFile, width = 8, height = 4)
