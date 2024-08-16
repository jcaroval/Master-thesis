######Packages#######
library(dplyr)
library(ggplot2) 
library(ggpubr)
library(tidyr) 
library(igraph) 
library(vegan) 
library(corrplot) 
library(ggbreak)
library(effects)
library(stargazer)
library(abdiv)
library(pheatmap) 

######Input of the dataset######
diversity_data <- read.csv(file.choose(), sep = ";")
diversity_data <- diversity_data %>%
  mutate(Bush = factor(Bush))
diversity_data <- diversity_data %>%
  mutate(Rock = factor(Rock))
diversity_data <- diversity_data %>%
  mutate(Cactus = factor(Cactus))
diversity_data <- diversity_data %>%
  mutate(Habitat = factor(Habitat))
diversity_data <- diversity_data %>%
  mutate(Sampling.site = factor(Sampling.site))
str(diversity_data)

######Distance Matrix and Grid Estimation######
#Define coordinates of the sampling spots (done with a 100x100 raster in Powerpoint based on the scheme that was provided from the sampling site)
coord <- data.frame(
  Point = c("A", "B", "C",  "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"),
  X = c(6.1, 6.2, 7.8, 9.8, 5.6, 3.8, 3.0, 4.8, 8.8, 2.4, 1.4, 6.4, 7.2, 7.6, 10.0, 4.7),
  Y = c(3.8, 4.6, 8.7, 10.1, 5.1, 2.6, 4.6, 7.8, 7.0, 9.1, 8.1, 6.5, 2.8, 2.1, 8.2, 6.4))

#Input of the distance matrix (based on the measurements taken at the EPT)
distance_matrix <- matrix(c(0.00,	1.80,	NA,	6.30,	1.70,	2.40,	3.10,	NA,	2.30,	NA,	NA,	2.90,	1.50,	3.20,	4.20,	4.50,
                            1.80,	0.00,	4.60,	4.80,	2.40,	3.90,	4.00,	NA,	1.20,	NA,	NA,	NA,	2.80,	3.70,	NA,	4.10,
                            NA,	4.60,	0.00,	2.30,	NA,	NA,	NA,	1.80,	NA,	3.20,	NA,	3.00,	NA,	NA,	5.30,	2.90,
                            6.30,	4.80,	2.30,	0.00,	NA,	NA,	NA,	4.00,	5.10,	NA,	NA,	3.40,	NA,	NA,	4.20,	6.10,
                            1.70,	2.40,	NA,	NA,	0.00,	1.90,	1.70,	NA,	3.50,	5.50,	4.20,	2.90,	2.60,	4.90,	5.10,	3.10,
                            2.40,	3.90,	NA,	NA,	1.90,	0.00,	2.20,	NA,	NA,	NA,	5.20,	4.80,	2.50,	4.60,	NA,	NA,
                            3.10,	4.00,	NA,	NA,	1.70,	2.20,	0.00,	NA,	NA,	4.80,	3.00,	4.00,	4.00,	NA,	NA,	2.80,
                            NA,	NA,	1.80,	4.00,	NA,	NA,	NA,	0.00,	NA,	2.00,	2.70,	NA,	NA,	NA,	NA,	NA,
                            2.30,	1.20,	NA,	5.10,	3.50,	NA,	NA,	NA,	0.00,	NA,	NA,	2.50,	2.80,	3.00,	2.00,	5.50,
                            NA,	NA,	3.20,	NA,	5.50,	NA,	4.80,	2.00,	NA,	0.00,	2.00,	NA,	NA,	NA,	NA,	2.50,
                            NA,	NA,	NA,	NA,	4.20,	5.20,	3.00,	2.70,	NA,	2.00,	0.00,	5.00,	NA,	NA,	NA,	1.80,
                            2.90,	NA,	3.00,	3.40,	2.90,	4.80,	4.00,	NA,	2.50,	NA,	5.00,	0.00,	NA,	NA,	NA,	3.20,
                            1.50,	2.80,	NA,	NA,	2.60,	2.50,	4.00,	NA,	2.80,	NA,	NA,	NA,	0.00,	2.10,	NA,	NA,
                            3.20,	3.70,	NA,	NA,	4.90,	4.60,	NA,	NA,	3.00,	NA,	NA,	NA,	2.10,	0.00,	4.70,	NA,
                            4.20,	NA,	5.30,	4.20,	5.10,	NA,	NA,	NA,	2.00,	NA,	NA,	NA,	NA,	4.70,	0.00,	NA,
                            4.50,	4.10,	2.90,	6.10,	3.10,	NA,	2.80,	NA,	5.50,	2.50,	1.80,	3.20,	NA,	NA,	NA,	0.00), nrow = 16, byrow = TRUE,
                          dimnames = list(coord$Point, coord$Point))
#Visualize and print incomplete distance matrix
plot(1, type="n", xlim=c(0.5, ncol(distance_matrix) + 0.5), 
     ylim=c(0.5, nrow(distance_matrix) + 0.5), xlab="", ylab="", axes = FALSE)
for (i in 1:nrow(distance_matrix)) {
  for (j in 1:ncol(distance_matrix)) {
    value <- distance_matrix[i, j]
    if (!is.na(value)) {
      if (value != round(value, 1)) {  # Überprüfe, ob der Wert gerundet wurde
        text(j, nrow(distance_matrix) - i + 1, round(value, 1), col="#3C93CD", cex = 1)
      } else {
        text(j, nrow(distance_matrix) - i + 1, round(value, 1), col="black", cex = 1)
      }
    } else {
      text(j, nrow(distance_matrix) - i + 1, "NA", col="#D55E00", cex = 1)
    }
  }
}
axis(1, at = 1:ncol(distance_matrix), labels = colnames(distance_matrix), tick = FALSE, pos = 1)
axis(2, at = 1:nrow(distance_matrix), labels = rev(rownames(distance_matrix)), las = 1, tick = FALSE, pos = 0.7)
title(main = "Incomplete distance matrix")
#Additionally exported as a svg-file and further optimized using Inkscape

# Convert the given coordinates to an igraph object
graph <- graph_from_data_frame(as.data.frame(t(combn(coord$Point, 2))))
V(graph)$X <- coord$X
V(graph)$X <- coord$X

#Predict the missing distances based on the igraph
for (i in 1:length(coord$Point)) {
  for (j in 1:length(coord$Point)) {
    if (is.na(distance_matrix[i, j]) && i != j) {
      # Calculate the euclidian distance between the coordinates of the sampling sites
      distance <- sqrt((coord$X[i] - coord$X[j])^2 + (coord$Y[i] - coord$Y[j])^2)
      distance_matrix[i, j] <- distance
      distance_matrix[j, i] <- distance
    }
  }
}

#Visualize and print complete distance matrix
plot(1, type="n", xlim=c(0.5, ncol(distance_matrix) + 0.5), 
     ylim=c(0.5, nrow(distance_matrix) + 0.5), xlab="", ylab="", axes = FALSE)
for (i in 1:nrow(distance_matrix)) {
  for (j in 1:ncol(distance_matrix)) {
    value <- distance_matrix[i, j]
    if (!is.na(value)) {
      if (value != round(value, 1)) {  # Überprüfe, ob der Wert gerundet wurde
        text(j, nrow(distance_matrix) - i + 1, round(value, 1), col="#3C93CD", cex = 1)
      } else {
        text(j, nrow(distance_matrix) - i + 1, round(value, 1), col="black", cex = 1)
      }
    } else {
      text(j, nrow(distance_matrix) - i + 1, "NA", col="#F49E3F", cex = 1)
    }
  }
}
axis(1, at = 1:ncol(distance_matrix), labels = colnames(distance_matrix), tick = FALSE, pos = 1)
axis(2, at = 1:nrow(distance_matrix), labels = rev(rownames(distance_matrix)), las = 1, tick = FALSE, pos = 0.7)
title(main = "Complete distance matrix")
#Additionally exported as a svg-file and further optimized using Inkscape

# Export the modified distance matrix as complete_distance_matrix
write.csv(distance_matrix, file = "compl_dist_matrix.csv", row.names = TRUE)

# MDS
mds_result <- cmdscale(distance_matrix, k = 2) # k = 2 für 2D-Darstellung
print(mds_result)

# Replace the X- and Y-coordinates with the results of the MDS
coord$X <- mds_result[, 1]
coord$Y <- mds_result[, 2]

# Look at the newly defined data
print(coord)
# Since the plot had the wrong ordination, the X values were multiplied by (-1).
# In addition, to look at positive values only, a specified value was added to the X and Y values.
coord$X <- -coord$X 
coord$X <- coord$X+6
coord$Y <- coord$Y+5

##stats#
cactus <- c(3, 6.4, 7.9, 6.4, 2, 5.2, 3.1, 4, 6.3, 4.8, 2.2, 1.5, 5.2, 3.1, 2.4, 5.9, 7.9, 2.5, 3.2, 7.3, 6.6)
rock <- c(2, 7.5, 7.7, 8.7, 4, 5.2, 6.4, 8.3)
bush <- c(7.5, 6.1, 6.5, 4.8, 5.1, 4.2, 9.6, 2.5, 5.5, 5.9, 6.7, 7.7, 6.6, 3.1, 4.1, 5.5, 5.6, 3.9, 2.4, 3.5, 5.1, 1.9, 1.2, 5.2, 2.4, 2, 6.7, 8.4)

mean_cactus <- mean(cactus)
mean_rock <- mean(rock)
mean_bush <- mean(bush)
median_cactus <- median(cactus)
median_rock <- median(rock)
median_bush <- median(bush)
sd_cactus <- sd(cactus)
sd_rock <- sd(rock)
sd_bush <- sd(bush)
iqr_cactus <- IQR(cactus)
iqr_rock <- IQR(rock)
iqr_bush <- IQR(bush)

test <- list(
  cactus = list(mean = mean_cactus, median = median_cactus, sd = sd_cactus, iqr = iqr_cactus),
  rock = list(mean = mean_rock, median = median_rock, sd = sd_rock, iqr = iqr_rock),
  bush = list(mean = mean_bush, median = median_bush, sd = sd_bush, iqr = iqr_bush)
)
str(test)

data <- data.frame(
  value = c(cactus, rock, bush),
  group = factor(rep(c("cactus", "rock", "bush"), c(length(cactus), length(rock), length(bush))))
)

shapiro_test_cactus <- shapiro.test(cactus)
shapiro_test_rock <- shapiro.test(rock)
shapiro_test_bush <- shapiro.test(bush)

shapiro <- list(
  cactus = shapiro_test_cactus$p.value,
  rock = shapiro_test_rock$p.value,
  bush = shapiro_test_bush$p.value
)

str(shapiro)

anova_result <- aov(value ~ group, data = data)
summary(anova_result)

posthoc <- TukeyHSD(anova_result)
posthoc

##measured & calculated stats
measured <- c(1.80, 6.30, 1.70, 2.40, 3.10, 2.30, 2.90, 1.50, 3.20, 4.20, 4.50, 4.60, 4.80, 
              2.40, 3.90, 4.00, 1.20, 2.80, 3.70, 4.10, 2.30, 1.80, 3.40, 3.20, 3.00, 4.00, 
              5.10, 5.30, 2.90, 4.20, 6.10, 1.90, 1.70, 5.20, 4.80, 2.50, 4.60, 3.50, 5.50, 
              4.20, 2.90, 2.60, 4.90, 5.10, 3.10, 2.20, 4.80, 3.00, 4.00, 4.00, 2.00, 2.70, 
              2.80, 2.50, 2.80, 3.00, 2.00, 5.50, 2.00, 3.20, 2.10, 1.80, 4.70, 2.50, 5.00)

calculated <- c(5.2, 4.2, 6.5, 6.4, 3.5, 5.9, 5.9, 1.9, 5.2, 4.2, 7.3, 6.3, 2.0, 6.4, 5.9, 6.6, 
                6.5, 9.6, 8.7, 7.5, 8.6, 7.7, 8.3, 2.8, 5.3, 6.7, 6.6, 8.4, 3.9, 3.7, 6.3, 5.2, 
                7.9, 4.1, 2.1, 5.5, 6.4, 5.2, 1.4, 6.7, 7.5, 4.8, 7.9, 8.7, 7.7, 7.9, 8.6, 8.6, 
                3.8, 4.6, 4.0, 6.1, 4.4, 5.2, 5.6)

mean_measured <- mean(measured)
mean_calculated <- mean(calculated)
median_measured <- median(measured)
median_calculated <- median(calculated)
sd_measured <- sd(measured)
sd_calculated <- sd(calculated)
iqr_measured <- IQR(measured)
iqr_calculated <- IQR(calculated)

list(
  measured = list(mean = mean_measured, median = median_measured, sd = sd_measured, iqr = iqr_measured),
  calculated = list(mean = mean_calculated, median = median_calculated, sd = sd_calculated, iqr = iqr_calculated)
)

shapiro_test_measured <- shapiro.test(measured)
shapiro_test_calculated <- shapiro.test(calculated)

list(
  measured = shapiro_test_measured$p.value,
  calculated = shapiro_test_calculated$p.value
)

t_test_result <- t.test(measured, calculated)
t_test_result

data <- data.frame(
  value = c(measured, calculated),
  group = factor(rep(c("measured", "calculated"), c(length(measured), length(calculated))))
)

kruskal_test_result <- kruskal.test(value ~ group, data = data)
kruskal_test_result

##all-distances statistics
all_dist <- c(1.80, 6.30, 1.70, 2.40, 3.10, 2.30, 2.90, 1.50, 3.20, 4.20, 4.50, 4.60, 4.80, 
  2.40, 3.90, 4.00, 1.20, 2.80, 3.70, 4.10, 2.30, 1.80, 3.40, 3.20, 3.00, 4.00, 
  5.10, 5.30, 2.90, 4.20, 6.10, 1.90, 1.70, 5.20, 4.80, 2.50, 4.60, 3.50, 5.50, 
  4.20, 2.90, 2.60, 4.90, 5.10, 3.10, 2.20, 4.80, 3.00, 4.00, 4.00, 2.00, 2.70, 
  2.80, 2.50, 2.80, 3.00, 2.00, 5.50, 2.00, 3.20, 2.10, 1.80, 4.70, 2.50, 5.00,
  5.2, 4.2, 6.50, 6.4, 3.5, 5.9, 5.9, 1.9, 4.2, 7.3, 6.3, 5.2, 2, 2.8, 6.4, 5.9,
  6.6, 7.5, 8.6, 6.5, 9.6, 8.7, 7.7, 8.3, 5.3, 6.7, 6.6, 2.1, 5.5, 6.4, 5.2, 1.4, 
  6.7, 7.5, 8.4, 3.9, 3.7, 6.3, 4.1, 5.2, 7.9, 7.9, 8.6, 8.6, 3.8, 4.6, 4, 5.2, 
  5.6, 6.1, 4.4, 4.8, 7.9, 8.7, 7.7)

mean(all_dist)
sd(all_dist)
median(all_dist)
IQR(all_dist)

# Grid visualization
g <- ggplot(coord, aes(x = X, y = Y)) +
  geom_polygon(data = data.frame(x = c(1.75, 10.25, 11.5, 0.25, 1), y = c(0.15, 0, 10, 9.9, 3)), # polygon was drawn based on the sketch of EPT
               aes(x = x, y = y), 
               fill = "#cdbda7", alpha = 0.2, 
               color = "black",
               linetype = "dashed") +
  geom_point(color = "#847a6c", shape = 16, size = 3, alpha = 0) +
  geom_text(aes(label = Point), size = 7, hjust = 0.5, vjust = 0.5) +
  labs(x = "[m]", y = "[m]", title = "EPT Sampling Grid", subtitle = "23°31'36''S    67°51'55''W    3592 masl") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 11, by = 1 ))+
  scale_y_continuous(breaks = seq(0, 11, by = 1))+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5,
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 14))

######Abundance Total ######
#create subset of data
totabundance <- diversity_data[,c(1,23:29)]
colnames(totabundance) <- c("Sampling.site", "Cephalobidae", "Panagrolaimidae", "Alaimidae", "Aphelenchoididae", "Qudsianematidae", "Anguinidae", "unidentified")

relabundance <- diversity_data[,c(1,31:37)]
colnames(relabundance) <- c("Sampling.site", "Cephalobidae", "Panagrolaimidae", "Alaimidae", "Aphelenchoididae", "Qudsianematidae", "Anguinidae", "unidentified")

col <- c("#999999","#CC79A7" , "#56B4E9","#6EED80",  "#F0E442","#E69F00", "#0072B2" )
#pivot-longer
totlong <- pivot_longer(totabundance, -Sampling.site, names_to = "Taxonomy", values_to = "count")
rellong <- pivot_longer(relabundance, -Sampling.site, names_to = "Taxonomy", values_to = "relabundance")

totlong$Sampling.site <- factor(totlong$Sampling.site, levels = c("EPT0I", "EPT0J", "EPT0F", "EPT0O", "EPT0D", "EPT0G", "EPT0E", "EPT0C", "EPT0A", "EPT0B", "EPT0H", "EPT0M", "EPT0P", "EPT0K", "EPT0N", "EPT0L"))
rellong$Sampling.site <- factor(rellong$Sampling.site, levels = c("EPT0I", "EPT0J", "EPT0F", "EPT0O", "EPT0D", "EPT0G", "EPT0E", "EPT0C", "EPT0A", "EPT0B", "EPT0H", "EPT0M", "EPT0P", "EPT0K", "EPT0N", "EPT0L"))

rellong$Taxonomy <- factor(rellong$Taxonomy, levels = c("unidentified", "Alaimidae", "Aphelenchoididae","Panagrolaimidae", "Anguinidae", "Qudsianematidae", "Cephalobidae"))
totlong$Taxonomy <- factor(totlong$Taxonomy, levels = c("unidentified", "Alaimidae", "Aphelenchoididae","Panagrolaimidae", "Anguinidae", "Qudsianematidae", "Cephalobidae"))
# Plot
tot <- ggplot(totlong, aes(x = Sampling.site, y = count))+
  geom_bar(aes(fill = Taxonomy), stat = "identity")+
  scale_fill_manual(values = col)+
  theme_minimal()+
  scale_x_discrete(name = "Sampling spot")+
  scale_y_continuous(name = "Count [n]",
                     breaks = seq(0,310, by = 50),
                     limits = c(0, 310),
                     expand = c(0,0.5))
rel <- ggplot(rellong, aes(x = Sampling.site, y = relabundance))+
  geom_bar(aes(fill = Taxonomy), stat = "identity")+
  scale_fill_manual(values = col)+
  theme_minimal()+
  scale_x_discrete(name = "Sampling spot")+
  scale_y_continuous(name = "Relative abundance [%]")
hab <- ggplot(data=rellong, aes(x = Sampling.site, y = relabundance))+
  scale_fill_manual(values = col)+
  theme_minimal()+
  scale_x_discrete(name = "Sampling spot")+
  scale_y_continuous(name = "Habitat",
                     breaks = seq(0.01, 1.01, by = 2))+
  theme(axis.title.x = element_text(margin = margin(t = 10))) #the icons were added manually after printing the graph

combined_habitat <-  ggarrange(
  tot + rremove("xlab"), rel +rremove("xlab"), hab, labels = c("a", "b", "c"),
  ncol = 1, nrow = 3,
  common.legend = TRUE, legend = "right")

#Statistics Abundance Total#
range_coverage <- range(diversity_data$total.extract)
mean_coverage <- mean(diversity_data$total.extract)
sd_coverage <- sd(diversity_data$total.extract)
median_coverage <- median(diversity_data$total.extract)
iqr_coverage <- IQR(diversity_data$total.extract)

print(mean_coverage)
print(sd_coverage)
print(median_coverage)
print(iqr_coverage)
print(range_coverage)

#the visualization was optimized using PowerPoint.
######Abundance Habitat ######
desired_order <- c("All.habitats", "Cactus&Bush", "Bush", "Bush&Rock", "Cactus", "Rock", "none") #sorted by the mean (highest to lowest)
diversity_data$Habitat <- factor(diversity_data$Habitat, levels = desired_order)
my_comparisons <- list (c("All.habitats","Bush"), c("All.habitats", "Bush&Rock"), c("All.habitats", "Cactus"), c("All.habitats", "Cactus&Bush"), c("All.habitats","none"), c("All.habitats","Rock"), c("Bush", "Bush&Rock"), c("Bush", "Cactus"), c("Bush", "Cactus&Bush"), c("Bush", "none"), c("Bush","Rock"), c("Bush&Rock", "Cactus"), c("Bush&Rock", "Cactus&Bush"), c("Bush&Rock","none"), c("Bush&Rock","Rock"), c("Cactus","Cactus&Bush"), c("Cactus","none"), c("Cactus","Rock"), c("Cactus&Bush", "none"), c("Cactus&Bush","Rock"), c("none","Rock"))
df_summary <- diversity_data %>%
  group_by(Habitat)%>%
  summarize(total.extract = mean(total.extract, na.rm = TRUE))
ggplot(data = df_summary, aes(x = Habitat, y = total.extract))+
  geom_bar(stat = "identity", position = position_dodge(), fill = "#457992", alpha = 1)+
  geom_jitter(data=diversity_data, fill = "orange", 
              shape = 21, size = 5, alpha = 0.8, width = 0.1)+
  stat_summary(data=diversity_data, fun.data = mean_cl_boot, #mean_cl_boot instead of mean_cl_normal because with the bootstaps you can obtain cl's for population mean without assuming normality
               geom = "errorbar",
               aes(width = 0.4),
               size = 0.5,
               alpha = 1,
               color = "darkblue")+
  scale_y_continuous(breaks = c(seq(0,300, by = 25)))+
  scale_x_discrete(labels = c("CactBushRock", "CactBush", "Bush", "BushRock", "Cact", "Rock", "none"))+
  ylab("Count [n]")+
  xlab("Habitat")+
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14, margin = margin(t=10)),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.15, color = "gray"),
        panel.grid.minor.y = NULL)

#statistics
#wilcox test nemabundance barplots
group1 <- diversity_data$Sampling.site[diversity_data$Habitat == "All.habitats"]
group2 <- diversity_data$Sampling.site[diversity_data$Habitat == "Cactus&Bush"]
group3 <- diversity_data$Sampling.site[diversity_data$Habitat == "Bush"]
group4 <- diversity_data$Sampling.site[diversity_data$Habitat == "Bush&Rock"]
group5 <- diversity_data$Sampling.site[diversity_data$Habitat == "Cactus"]
group6 <- diversity_data$Sampling.site[diversity_data$Habitat == "Rock"]
group7 <- diversity_data$Sampling.site[diversity_data$Habitat == "none"]

combined_data <- data.frame(
  variable = c(group1, group2, group3, group4, group5, group6, group7),
  group = factor(c(
    rep("All.habitats", length(group1)),
    rep("Cactus&Bush", length(group2)),
    rep("Bush", length(group3)),
    rep("Bush&Rock", length(group4)),
    rep("Cactus", length(group5)),
    rep("Rock", length(group6)),
    rep("none", length(group7))
  ))
)

kruskal.test(variable ~group, data = combined_data)

######Diversity (barcoding data)######
##alpha-div##
y_axis_limits <- c(0, 5)
#families
totabundance_pool <- diversity_data[,c(1,23:28)]
df_pool <- as.data.frame(totabundance_pool)

mat_pool <- as.matrix(df_pool[, -1])
richness_pool <- apply(mat_pool, 1, richness)
shannon_pool <- apply(mat_pool, 1, shannon)
invsimpson_pool <- apply(mat_pool, 1, invsimpson)

div_data_pool <- data.frame(Sampling.site = LETTERS[1:16], Richness = richness_pool, InvSimpson = invsimpson_pool)

fam <- ggplot(div_data_pool, aes(x = Sampling.site)) +
  geom_point(aes(y = InvSimpson), color = "darkorange", alpha = 0.8, shape = 16, size = 4) +
  geom_point(aes(y = Richness), color = "#457992", alpha = 0.8, shape = 16, size = 4) +
  labs(title = "Diversity Measures Family",
       x = "Sampling Spot",
       y = "Diversity Index",
       color = "Measure") +
  theme_minimal()+
  ylim(y_axis_limits) +
  scale_color_manual(values = c("#457992", "darkorange"),
                     labels = c("Richness", "Inverse Simpson"))

##Genera
totabundance_gen <- diversity_data[,c(1,14,15,17)]
df_single <- as.data.frame(totabundance_single)

mat_single <- as.matrix(df_single[, -1])
richness_single <- apply(mat_single, 1, richness)
shannon_single <- apply(mat_single, 1, shannon)
invsimpson_single <- apply(mat_single, 1, invsimpson)

div_data_single <- data.frame(Sampling.site = LETTERS[1:16], Richness = richness_single, InvSimpson = invsimpson_single)

gen <- ggplot(div_data_single, aes(x = Sampling.site)) +
  geom_point(aes(y = InvSimpson), color = "darkorange", alpha = 0.8, shape = 16, size = 4) +
  geom_point(aes(y = Richness), color = "#457992", alpha = 0.8, shape = 16, size = 4) +
  labs(title = "Diversity Measures Genus",
       x = "Sampling Spot",
       y = "Diversity Index",
       color = "Measure") +
  theme_minimal()+
  ylim(y_axis_limits) +
  scale_color_manual(values = c("#457992", "darkorange"),
                     labels = c("Richness", "Inverse Simpson"))

##species
totabundance_species <- diversity_data[,c(1,39:41)]
df_species <- as.data.frame(totabundance_species)

mat_species <- as.matrix(df_species[, -1])
richness_species <- apply(mat_species, 1, richness)
shannon_species <- apply(mat_species, 1, shannon)
invsimpson_species <- apply(mat_species, 1, invsimpson)

div_data_species <- data.frame(Sampling.site = LETTERS[1:16], Richness = richness_species, InvSimpson = invsimpson_species)

spec <- ggplot(div_data_species, aes(x = Sampling.site)) +
  geom_point(aes(y = InvSimpson), color = "darkorange", alpha = 0.8, shape = 16, size = 4) +
  geom_point(aes(y = Richness), color = "#457992", alpha = 0.8, shape = 16, size = 4) +
  labs(title = "Diversity Measures Species",
       x = "Sampling Spot",
       y = "Diversity Index",
       color = "Measure") +
  theme_minimal()+
  ylim(y_axis_limits) +
  scale_color_manual(values = c("#457992", "darkorange"),
                     labels = c("Richness", "Inverse Simpson"))


##haplotypes
totabundance_hap <- diversity_data[,c(1,42:44)]
df_hap <- as.data.frame(totabundance_hap)

mat_hap <- as.matrix(df_hap[, -1])
richness_hap <- apply(mat_hap, 1, richness)
shannon_hap <- apply(mat_hap, 1, shannon)
invsimpson_hap <- apply(mat_hap, 1, invsimpson)

div_data_hap <- data.frame(Sampling.site = LETTERS[1:16], Richness = richness_hap, InvSimpson = invsimpson_hap)

hap <- ggplot(div_data_hap, aes(x = Sampling.site)) +
  geom_point(aes(y = InvSimpson), color = "darkorange", alpha = 0.8, shape = 16, size = 4) +
  geom_point(aes(y = Richness), color = "#457992", alpha = 0.8, shape = 16, size = 4) +
  labs(title = "Diversity Measures Haplotypes",
       x = "Sampling Spot",
       y = "Diversity Index",
       color = "Measure") +
  theme_minimal()+
  ylim(y_axis_limits) +
  scale_color_manual(values = c("#457992", "darkorange"),
                     labels = c("Richness", "Inverse Simpson"))

ggarrange(fam, gen, spec, hap, nrow = 2, ncol = 2)

richness_data <- data.frame(
  Sampling.site = rep(LETTERS[1:16], 4),
  Richness = c(richness_pool, richness_single, richness_species, richness_hap),
  Measure = rep(c("Family", "Genus", "Species", "Haplotypes"), each = 16)
)

invsimpson_data <- data.frame(
  Sampling.site = rep(LETTERS[1:16], 4),
  InvSimpson = c(invsimpson_pool, invsimpson_single, invsimpson_species, invsimpson_hap),
  Measure = rep(c("Family", "Genus", "Species", "Haplotypes"), each = 16)
)

richness_data$Measure <- factor(richness_data$Measure, levels = c("Family", "Genus", "Species", "Haplotypes"))
invsimpson_data$Measure <- factor(invsimpson_data$Measure, levels = c("Family", "Genus", "Species", "Haplotypes"))

##stats
richness_summary <- richness_data %>%
  group_by(Measure) %>%
  summarise(
    Mean = mean(Richness, na.rm = TRUE),
    SD = sd(Richness, na.rm = TRUE),
    Median = median(Richness, na.rm = TRUE),
    IQR = IQR(Richness, na.rm = TRUE)
  )

invsimpson_summary <- invsimpson_data %>%
  group_by(Measure) %>%
  summarise(
    Mean = mean(InvSimpson, na.rm = TRUE),
    SD = sd(InvSimpson, na.rm = TRUE),
    Median = median(InvSimpson, na.rm = TRUE),
    IQR = IQR(InvSimpson, na.rm = TRUE)
  )

richness_summary
invsimpson_summary

shapiro_tests <- richness_data %>%
  group_by(Measure) %>%
  summarise(p_value = shapiro.test(Richness)$p.value)

print(shapiro_tests)

shapiro_tests_invsimpson <- invsimpson_data %>%
  group_by(Measure) %>%
  summarise(p_value = shapiro.test(InvSimpson)$p.value)

print(shapiro_tests_invsimpson)

kruskal_test_richness <- kruskal.test(Richness ~ Measure, data = richness_data)
print(kruskal_test_richness)

kruskal_test_invsimpson <- kruskal.test(InvSimpson ~ Measure, data = invsimpson_data)
print(kruskal_test_invsimpson)

dunn_test_richness <- dunnTest(Richness ~ Measure, data = richness_data)
print(dunn_test_richness)

dunn_test_invsimpson <- dunnTest(InvSimpson ~ Measure, data = invsimpson_data)
print(dunn_test_invsimpson)

##beta-div##
##Families
increasing <- colorRampPalette(c("#f7f7f7","#fee090","#fdae61"))(300)

totabundance_pool <- diversity_data[,c(1,23:28)]
beta_pool <- as.data.frame(totabundance_pool)
beta_counts_pool <- beta_pool[,-1]
jaccard_pool <- vegdist(beta_counts_pool, method = "jaccard", binary = T)
print(jaccard_pool)

pheatmap(as.matrix(jaccard_pool),
         display_numbers = TRUE,
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 20,
         border_color = "white",
         angle_col = "0",
         color = increasing,
         fontsize_number = 15)


##Genera
totabundance_single <- diversity_data[,c(1,14:19)]
beta_single <- totabundance_single
beta_counts_single <- beta_single[,-1]
jaccard_single <- vegdist(beta_counts_single, method = "jaccard", binary = T)
print(jaccard_single)

pheatmap(as.matrix(jaccard_single),
         display_numbers = TRUE,
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 20,
         border_color = "white",
         angle_col = "0",
         color = increasing,
         fontsize_number = 15)

##Species
beta_species <- as.data.frame(totabundance_species)
beta_counts_species <- beta_species[,-1]
jaccard_species <- vegdist(beta_counts_species, method = "jaccard", binary = T)
print(jaccard_species)

pheatmap(as.matrix(jaccard_species),
         display_numbers = TRUE,
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 20,
         border_color = "white",
         angle_col = "0",
         color = increasing,
         fontsize_number = 15)

##haplotypes
beta_hap <- as.data.frame(totabundance_hap)
beta_counts_hap <- beta_hap[,-1]
jaccard_hap <- vegdist(beta_counts_hap, method = "jaccard", binary = T)
print(jaccard_hap)

pheatmap(as.matrix(jaccard_hap),
         display_numbers = TRUE, 
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 20,
         border_color = "white",
         angle_col = "0",
         color = increasing,
         fontsize_number = 15)

######Nucleotide Diversity######
pi <- read.table(file.choose(), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

##Statistics
min_value <- min(pi$avg_pi)
max_value <- max(pi$avg_pi)
mean_value <- mean(pi$avg_pi)
sd_value <- sd(pi$avg_pi)
median_value <- median(pi$avg_pi)
iqr_value <- IQR(pi$avg_pi)

summary_stats <- pi %>%
  summarize(
    Min = min(avg_pi),
    Max = max(avg_pi),
    Mean = mean(avg_pi),
    SD = sd(avg_pi),
    Median = median(avg_pi),
    IQR = IQR(avg_pi)
  )
print(summary_stats)


grouped_summary <- pi %>%
  # filter(avg_pi > 0) %>%
  group_by(pop) %>%
  summarize(
    Min = min(avg_pi),
    Max = max(avg_pi),
    Mean = mean(avg_pi),
    SD = sd(avg_pi)
  )
print(grouped_summary)

normality_tests <- pi %>%
  group_by(pop) %>%
  summarize(
    p_value = shapiro.test(avg_pi)$p.value
  )
print(normality_tests)

kruskal_result <- kruskal.test(avg_pi ~ pop, data = pi)
print(kruskal_result)

anova_result <- aov(avg_pi ~ pop, data = pi)

shapiro_test <- shapiro.test(residuals(anova_result))
print(shapiro_test) #####keine Normalverteilung --> keine ANOVA, sondern Kruskal Wallis/Dunn Test

kruskal_result <- kruskal.test(avg_pi ~ pop, data = pi_factor)
print(kruskal_result) #####keine Signifikanz, also kein post-hoc (kein Dunn-Test)

pi_factor <- pi %>%
  mutate(pop = as.factor(pop))
dunn_result <- dunnTest(avg_pi ~ pop, data = pi_factor, method = "bonferroni")
print(dunn_result$res)

manual_significance <- data.frame(
  group1 = c("EPT-0A", "EPT-0C", "EPT-0D"),
  group2 = c("EPT-0N", "EPT-0N", "EPT-0N"),
  y_position = c(max(pi_factor$avg_pi) + 0.020, max(pi_factor$avg_pi) + 0.005, max(pi_factor$avg_pi) -0.015),
  annotations = c("*", "**", "**")
)

##Plot
nuc_div_a <- ggplot(pi, aes(x = pop, y = avg_pi)) +
  geom_violin(fill = "grey", size = 0.2, alpha = 0.4, trim = FALSE, scale = "width") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2, size = 1, color = "#457992") +
  stat_summary(fun = mean, geom = "point", color = "#457992", size = 3, shape = 16) +
  labs(x = "Sampling spot", y = "Average nucleotide diversity") +
  scale_y_continuous(breaks = c(seq (0, 0.26, by = 0.02)))+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    panel.grid.major.y = element_line(size = 0.5, color = "grey80"),
    panel.grid.minor.y = element_line(size = 0.25, color = "grey90")
  ) +
  geom_signif(
    comparisons = list(c("EPT-0A", "EPT-0N"), c("EPT-0C", "EPT-0N"), c("EPT-0D", "EPT-0N")),
    annotations = manual_significance$annotations,
    y_position = manual_significance$y_position,
    textsize = 7
  )

nuc_div_b <- ggplot(pi_factor, aes(x = pop, y = avg_pi)) +
  stat_summary(fun = mean, geom = "point", color = "#457992", size = 4, shape = 16) +
  labs(x = "Sampling spot", y = "Average nucleotide diversity") +
  coord_cartesian(ylim = c(0.02, 0.032)) +
  scale_y_continuous(breaks = seq(0, 0.038, by = 0.002)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid.major.y = element_line(size = 0.5, color = "grey80"),
    panel.grid.minor.y = element_line(size = 0.25, color = "grey90")
  )

ggarrange(print(nuc_div_a), print(nuc_div_b), nrow = 2, ncol = 1, heights = c(0.7, 0.3), labels = c("a", "b"))

##### Nucleotide divergence between populations (dxy)
dxy <- read.table(file.choose(), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

dxy_summary <- dxy %>%
  group_by(pop1, pop2) %>%
  summarise(avg_dxy = mean(avg_dxy, na.rm = TRUE), .groups = 'drop')

dxy_summary_reversed <- dxy_summary %>%
  rename(pop1 = pop2, pop2 = pop1)

dxy_summary_combined <- bind_rows(dxy_summary, dxy_summary_reversed)

dxy_matrix <- dxy_summary_combined %>%
  spread(key = pop2, value = avg_dxy)

dxy_matrix <- as.data.frame(dxy_matrix)
rownames(dxy_matrix) <- dxy_matrix$pop1
dxy_matrix$pop1 <- NULL

dxy_long <- melt(as.matrix(dxy_matrix))

pheatmap(as.matrix(dxy_matrix), 
         display_numbers = TRUE,
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 20,
         border_color = "white",
         angle_col = "0",
         color = increasing, 
         fontsize_number = 15,
         number_format = "%.4g")

diversity_Laura <- c(#0.0402963,
                     #0.0408627,
                     #0.106631,
                     #0.206227,
                     #0.138237,
                     #0.221082,
                     #0.107891,
                     0.0545159,
                     0.0167081,
                     0.0488625,
                     0.00241247)
diversity_Laura_sex <- c( 0.0545159,
                     0.0167081,
                     0.0488625,
                     0.00241247)
diversity_Laura_asex <- c(0.0402963,
                    0.0408627,
                    0.106631,
                    0.206227,
                    0.138237,
                    0.221082,
                    0.107891)
diversity_Julia <- pi[, 5]

dataframe <- factor(c(rep("A. emmatus", length(diversity_Julia)),
                      rep("Panagrolaimidae (sexual)", length(diversity_Laura_sex)),
                      rep("Panagrolaimidae (asexual)", length(diversity_Laura_asex))))
diversity <- c(diversity_Julia, diversity_Laura_sex, diversity_Laura_asex)
data <- data.frame(species = dataframe, diversity = diversity)

summary_stats <- data %>%
  group_by(species) %>%
  summarize(mean = mean(diversity),
            sd = sd(diversity),
            median = median(diversity),
            IQR = IQR(diversity),
            .groups = 'drop')
summary_stats

diversity_combined <- c(diversity_Julia, diversity_Laura_sex, diversity_Laura_asex)
group <- factor(c(rep("Julia", length(diversity_Julia)), rep("Laura_sex", length(diversity_Laura_sex)), rep("Laura_asex", length(diversity_Laura_asex))))

kruskal_test_result <- kruskal.test(diversity_combined ~ group)

p_value <- kruskal_test_result$p.value
p_value

dunn_result <- dunnTest(diversity_combined, group, method = "bonferroni")
print(dunn_result$res)

ggplot(data, aes(x= species, y= diversity, fill = species))+
  geom_boxplot(color="black", alpha = 0.8)+
  geom_signif(comparisons=list(c("A. emmatus","Panagrolaimidae (sexual)")), annotations = paste0("NS"), y_position = 0.24, textsize =5)+
  geom_signif(comparisons=list(c("A. emmatus","Panagrolaimidae (asexual)")), annotations = paste0("***"), y_position = 0.225, textsize =7)+
  labs(#title = paste("p-value:", format(p_value, digits = 3)),
       x = "Species",
       y = "Average nucleotide diversity")+
  scale_fill_manual(values = c("A. emmatus" = "darkorange", 
                               "Panagrolaimidae (sexual)" = "#457992", 
                               "Panagrolaimidae (asexual)" = "#457992")) +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )

######Correlation (Dependence on Habitat and Nematodes)######
#####Families
D = data.frame(lapply(diversity_data[,c(3:5,23:28)], as.integer))
colnames(D) <- c("Cactus", "Bush", "Rock", "Cephalobidae", "Panagrolaimidae", "Alaimidae", "Aphelenchoididae", "Qudsianematidae", "Anguinidae")
M <- cor(D, method = "spearman")
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp <- cor.test(mat[,i], mat[,j],...)
      p.mat[i,j] <- p.mat [j,i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat}
p.mat <- cor.mtest(D)
p.mat
corrplot(M, method = "circle", mar = c(0,0,2,0), diag = FALSE, outline = TRUE, bg = "transparent", 
         addgrid.col = "transparent", type = "lower", tl.col="black", tl.srt=90,
         p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.col = "white", pch.cex = 2, col = COL2("RdBu",200))
#####Genera
D = data.frame(lapply(diversity_data[,c(3:5,14:19)], as.integer))
colnames(D) <- c("Cactus", "Bush", "Rock", "Acrobeles", "Acrobeloides", "Aphelenchoides", "Stegelletina", "Ditylenchus", "Panagrolaimus")
M <- cor(D, method = "spearman")
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp <- cor.test(mat[,i], mat[,j],...)
      p.mat[i,j] <- p.mat [j,i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat}
p.mat <- cor.mtest(D)
p.mat
corrplot(M, method = "circle", mar = c(0,0,2,0), diag = FALSE, outline = TRUE, bg = "transparent", 
         addgrid.col = "transparent", type = "lower", tl.col="black", tl.srt=90,
         p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.col = "white", pch.cex = 2, col = COL2("RdBu",200))

#####Species
D = data.frame(lapply(diversity_data[,c(3:5,39:41)], as.integer))
colnames(D) <- c("Cactus", "Bush", "Rock", "A. emmatus", "A. mariannae", "complex")
M <- cor(D, method = "spearman")
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp <- cor.test(mat[,i], mat[,j],...)
      p.mat[i,j] <- p.mat [j,i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat}
p.mat <- cor.mtest(D)
p.mat
corrplot(M, method = "circle", mar = c(0,0,2,0), diag = FALSE, outline = TRUE, bg = "transparent", 
         addgrid.col = "transparent", type = "lower", tl.col="black", tl.srt=90,
         p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.col = "white", pch.cex = 2, col = COL2("RdBu",200))


#####Haplotypes
D = data.frame(lapply(diversity_data[,c(3:5,42:44)], as.integer))
colnames(D) <- c("Cactus", "Bush", "Rock", "Hap A", "Hap B", "Hap C")
M <- cor(D, method = "spearman")
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp <- cor.test(mat[,i], mat[,j],...)
      p.mat[i,j] <- p.mat [j,i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat}
p.mat <- cor.mtest(D)
p.mat
corrplot(M, method = "circle", mar = c(0,0,2,0), diag = FALSE, outline = TRUE, bg = "transparent", 
         addgrid.col = "transparent", type = "lower", tl.col="black", tl.srt=90,
         p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.col = "white", pch.cex = 2, col = COL2("RdBu",200))

######NMDS######
#Subset the dataframe on which to base the ordination (dataframe 1)
data_single <- diversity_data[,c(14:19)]
data_pool <- diversity_data[,c(23:28)]
data_species <- diversity_data[,c(39:41)]
data_hap <- diversity_data[,c(42:44)]
data_hap_clean <- data_hap[rowSums(data_hap) > 0, ]

#Identify the columns that contain descriptive/environmental data (dataframe 2)
data_2 <- diversity_data[1:5]
data_2_noGKM <- data_2[-c(7, 11, 13), ]

NMDS_single <- metaMDS(data_single, distance = "bray", k = 2, trace=FALSE, trymax = 100)
NMDS_pool <- metaMDS(data_pool, distance = "bray", k = 2, trace=FALSE, trymax = 100)
NMDS_spec <- metaMDS(data_species, distance = "bray", k = 2, trace=FALSE, trymax = 100)
NMDS_hap <- metaMDS(data_hap_clean, distance = "bray", k = 2, trace=FALSE, trymax = 100) 
stressplot(NMDS_single)
stressplot(NMDS_pool)
stressplot(NMDS_spec)
stressplot(NMDS_hap)


#Visualization of the NMDS
x_limits_single <- c(-1.1, 1.5)
y_limits_single <- c(-1.6, 0.5)
x_limits_pool <- c(-1.5, 1.5)
y_limits_pool <- c(-1.6, 1)
x_limits_spec <- c(-1,1.5)
y_limits_spec <- c(-2,1)
x_limits_hap <- c(-1,1.5)
y_limits_hap <- c(-0.5,1)
stress_value_single <- NMDS_single$stress
stress_value_pool <- NMDS_pool$stress
stress_value_spec <- NMDS_spec$stress
stress_value_hap <- NMDS_hap$stress
shape=c(4,4,4,4,4,4)
txt=c("Rock", "Bush", "Cactus")
par(cex.axis = 1.4)
par(cex.lab = 1.4)
par(cex.main = 1.4)
par(mar = c(4.5,5,2,1)) #margin bottom, left, top, right

##Families##
plot(NMDS_pool$species, pch = shape[data_2$Sampling.site],
     cex=1.2, main = "EPT community composition pool", xlab = "NMDS1", ylab = "NMDS2",
     xlim = x_limits_pool, ylim = y_limits_pool, 
     cex.lab = 1.5, 
     cex.axis = 1.3,
     cex.main = 1.5,
     bty = "n")
abline(v = 0, lty = "dashed", col = "grey", lwd = 2)
abline(h = 0, lty = "dashed", col = "grey", lwd = 2)
grid(nx = NULL, ny = NULL, col = "grey", lty = "solid", equilogs = TRUE, lwd = 0.3)
ordihull(NMDS_pool, groups = data_2$Bush, draw = "polygon", lty = c(0,1), col = c("transparent","#E69F00"), alpha = 0.5)
ordihull(NMDS_pool, groups = data_2$Cactus, draw = "polygon", lty = c(0,1), col = c("transparent", "#33a02c"), alpha = 0.5)
ordihull(NMDS_pool, groups = data_2$Rock, draw = "polygon", lty = c(0,1), col = c("transparent", "#999999"), alpha = 0.5)
points(NMDS_pool,display="species", pch=4, col="black", cex=2)
text(NMDS_pool, display="species",labels=c("Cephalobidae", "Panagrolaimidae", "Alaimidae","Aphelenchoididae","Qudsianematidae","Anguinidae"), pos=1,offset=0.5,cex=2)
text(x = x_limits_pool[1], y = y_limits_pool[2], labels = paste("stress:", round(stress_value_pool, 3)),
     adj = c(0, 1), cex = 2)  # Vergrößert die Schriftgröße des Stress-Textes
legend('bottomright', txt, pch = c(16, 16, 16), col=c("#999999","#E69F00","#33a02c"), cex=1.4, pt.cex = 2, bty= "y")

##Genera##
plot(NMDS_single$species, pch = shape[data_2$Sampling.site],
     cex=1.2, main = "EPT community composition Rhabditida", xlab = "NMDS1", ylab = "NMDS2",
     xlim = x_limits_single, ylim = y_limits_single,
     cex.lab = 1.5,   
     cex.axis = 1.3,  
     cex.main = 1.5, 
     bty = "n")
abline(v = 0, lty = "dashed", col = "grey", lwd = 2)
abline(h = 0, lty = "dashed", col = "grey", lwd = 2)
grid(nx = NULL, ny = NULL, col = "grey", lty = "solid", equilogs = TRUE, lwd = 0.3)
ordihull(NMDS_single, groups = data_2$Bush, draw = "polygon", lty = c(0,1), col = c("transparent","#E69F00"), alpha = 0.5)
ordihull(NMDS_single, groups = data_2$Cactus, draw = "polygon", lty = c(0,1), col = c("transparent", "#33a02c"), alpha = 0.5)
ordihull(NMDS_single, groups = data_2$Rock, draw = "polygon", lty = c(0,1), col = c("transparent", "#999999"), alpha = 0.5)
points(NMDS_single,display="species", pch=4, col="black", cex=2)
text(NMDS_single, display="species",labels=c("Acrobeles", "Acrobeloides", "Aphelenchoides","Stegelletina","Ditylenchus","Panagrolaimus"), pos=1,offset=0.5,cex=2)
text(x = x_limits_single[1], y = y_limits_single[2], labels = paste("stress:", round(stress_value_single, 3)),
     adj = c(0, 1), cex = 2)

##species
plot(NMDS_spec$species, pch = shape[data_2$Sampling.site],
     cex=1.2, main = "EPT community composition Rhabditida", xlab = "NMDS1", ylab = "NMDS2",
     xlim = x_limits_spec, ylim = y_limits_spec,
     cex.lab = 1.5,  
     cex.axis = 1.3, 
     cex.main = 1.5, 
     bty = "n")
abline(v = 0, lty = "dashed", col = "grey", lwd = 2)
abline(h = 0, lty = "dashed", col = "grey", lwd = 2)
grid(nx = NULL, ny = NULL, col = "grey", lty = "solid", equilogs = TRUE, lwd = 0.3)
ordihull(NMDS_spec, groups = data_2$Bush, draw = "polygon", lty = c(0,1), col = c("transparent","#E69F00"), alpha = 0.5)
ordihull(NMDS_spec, groups = data_2$Cactus, draw = "polygon", lty = c(0,1), col = c("transparent", "#33a02c"), alpha = 0.5)
ordihull(NMDS_spec, groups = data_2$Rock, draw = "polygon", lty = c(0,1), col = c("transparent", "#999999"), alpha = 0.5)
points(NMDS_spec,display="species", pch=4, col="black", cex=2)
text(NMDS_spec, display="species",labels=c("A. emmatus", "A. mariannae", "complex"), pos=1,offset=0.5,cex=2)
text(x = x_limits_spec[1], y = y_limits_spec[2], labels = paste("stress:", round(stress_value_spec, 3)),
     adj = c(0, 1), cex = 2)

### Haplotypes
plot(NMDS_hap$species, pch = shape[data_2_noGKM$Sampling.site],
     cex=1.2, main = "EPT community composition hap", xlab = "NMDS1", ylab = "NMDS2",
     xlim = x_limits_hap, ylim = y_limits_hap, 
     cex.lab = 1.5,  
     cex.axis = 1.3, 
     cex.main = 1.5, 
     bty = "n")
abline(v = 0, lty = "dashed", col = "grey", lwd = 2)
abline(h = 0, lty = "dashed", col = "grey", lwd = 2)
grid(nx = NULL, ny = NULL, col = "grey", lty = "solid", equilogs = TRUE, lwd = 0.3)
ordihull(NMDS_hap, groups = data_2_noGKM$Bush, draw = "polygon", lty = c(0,1), col = c("transparent","#E69F00"), alpha = 0.5)
ordihull(NMDS_hap, groups = data_2_noGKM$Cactus, draw = "polygon", lty = c(0,1), col = c("transparent", "#33a02c"), alpha = 0.5)
ordihull(NMDS_hap, groups = data_2_noGKM$Rock, draw = "polygon", lty = c(0,1), col = c("transparent", "#999999"), alpha = 0.5)
points(NMDS_hap,display="species", pch=4, col="black", cex=2)
text(NMDS_hap, display="species",labels=c("Hap A", "Hap B", "Hap C"), pos=1,offset=0.5,cex=2)
text(x = x_limits_hap[1], y = y_limits_hap[2], labels = paste("stress:", round(stress_value_hap, 3)),
     adj = c(0, 1), cex = 2)  # Vergrößert die Schriftgröße des Stress-Textes

#The plots were visually further optimized using Inkscape. Namely, the "transparent" hull layers of the absent data were removed completely from the plot. 
#This was done because even when set to "transparent", they remained visible and thereby interfered with the data visualization (presence data only). 
#The removal of the absent data hulls did not change anything to the remaining plot.

##Statistical Analysis##
##Genera##
fit_bray_combined <- adonis2(data_single ~ Bush + Rock + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_single ~ Bush + Cactus + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_single ~ Rock + Bush + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_single ~ Rock + Cactus + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_single ~ Cactus + Rock + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_single ~ Cactus + Bush + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
##Families##
fit_bray_combined <- adonis2(data_pool ~ Bush + Rock + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_pool ~ Bush + Cactus + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_pool ~ Rock + Bush + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_pool ~ Rock + Cactus + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_pool ~ Cactus + Rock + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_pool ~ Cactus + Bush + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 

# Species
fit_bray_combined <- adonis2(data_species ~ Bush + Rock + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_species ~ Bush + Cactus + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_species ~ Rock + Bush + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_species ~ Rock + Cactus + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_species ~ Cactus + Rock + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_species ~ Cactus + Bush + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined 

#Hap
fit_bray_combined <- adonis2(data_hap_clean ~ Bush + Rock + Cactus, data=data_2_noGKM, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_hap_clean ~ Bush + Cactus + Rock, data=data_2_noGKM, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_hap_clean ~ Rock + Bush + Cactus, data=data_2_noGKM, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_hap_clean ~ Rock + Cactus + Bush, data=data_2_noGKM, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_hap_clean ~ Cactus + Rock + Bush, data=data_2_noGKM, permutations = 999, method = "bray")
fit_bray_combined 
fit_bray_combined <- adonis2(data_hap_clean ~ Cactus + Bush + Rock, data=data_2_noGKM, permutations = 999, method = "bray")
fit_bray_combined 
#The plots were visually further optimized using Inkscape. Namely, the "transparent" hull layers of the absent data were removed completely from the plot. 
#This was done because even when set to "transparent", they remained visible and thereby interfered with the data visualization (presence data only). 
#The removal of the absent data hulls did not change anything to the remaining plot.


######GLMs######
##Families##
glm_Alaimidae <- glm(Alaimidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) 
glm_Anguinidae <- glm(Anguinidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) 
glm_Aphelenchoididae <- glm(Aphelenchoididae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) 
glm_Cephalobidae <- glm(Cephalobidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
glm_Panagrolaimidae <- glm(Panagrolaimidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) 
glm_Qudsianematidae <- glm(Qudsianematidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_Alaimidae)
summary(glm_Anguinidae)
summary(glm_Aphelenchoididae)
summary(glm_Cephalobidae)
summary(glm_Panagrolaimidae)
summary(glm_Qudsianematidae)
##Genera##
glm_Acrobeles <- glm(Acrobeles_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
glm_Acrobeloides <- glm(Acrobeloides_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
glm_Aphelenchoides <- glm(Aphelenchoides_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) 
glm_Stegelletina <- glm(Stegelletina_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
glm_Ditylenchus <- glm(Ditylenchus_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) 
glm_Panagrolaimus <- glm(Panagrolaimus_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_Acrobeles)
summary(glm_Ditylenchus)
summary(glm_Aphelenchoides)
summary(glm_Acrobeloides)
summary(glm_Stegelletina)
summary(glm_Panagrolaimus)
#Concatenate all data into one table
#here, issues occured when trying to display all at once, which is why two datasets were created each. They were then manually concatenated into one table each.
#Families
stargazer(glm_Alaimidae, glm_Anguinidae, glm_Aphelenchoididae, 
          type = "latex", 
          out = "models_summary_part1.html",
          title = "GLM-families-1",
          align = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("Alaimidae", "Anguinidae", "Aphelenchoididae"),
          covariate.labels = c("Bush", "Cactus", "Rock"))

stargazer(glm_Cephalobidae, glm_Panagrolaimidae, glm_Qudsianematidae, 
          type = "latex", 
          out = "models_summary_part2.html",
          title = "GLM-families-2",
          align = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("Cephalobidae", "Panagrolaimidae", "Qudsianematidae"),
          covariate.labels = c("Bush", "Cactus", "Rock"))
#Genera
stargazer(glm_Acrobeles, glm_Acrobeloides, glm_Aphelenchoides,
          type = "html", 
          out = "models_summary_individuals_part1.html",
          title = "GLM-genera-1",
          align = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("Acrobeles", "Acrobeloides", "Aphelenchoides"),
          covariate.labels = c("Bush", "Cactus", "Rock"))

stargazer(glm_Stegelletina, glm_Ditylenchus, glm_Panagrolaimus,
          type = "html", 
          out = "GLM-genera-2",
          title = "Zusammenfassung mehrerer GLMs",
          align = TRUE,
          dep.var.labels = c("Cehalobidae", "Alaimidae", "Panagrolaimidae"),
          column.labels = c("Stegelletina", "Ditylenchus", "Panagrolaimus"),
          covariate.labels = c("Bush", "Cactus", "Rock"))

#Plot the GLMs that are build on at least 20 individuals
#the plots were then saved as SVG-files and were concatenated into one file using Inkscape.
effCephalobidae <- allEffects(glm_Cephalobidae)
effQudsianematidae <- allEffects(glm_Qudsianematidae)
effAcrobeles <- allEffects(glm_Acrobeles)
plot(effCephalobidae)
plot(effQudsianematidae)
plot(effAcrobeles)

glm_A.emmatus <- glm(A.emmatus ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_A.emmatus)
glm_A.mariannae <- glm(A.mariannae ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_A.mariannae)
glm_complex <- glm(complex ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_complex)

stargazer(glm_A.emmatus, glm_A.mariannae, glm_complex, 
          type = "latex", 
          out = "models_summary_part1.html",
          title = "GLM-families-1",
          align = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("A.emmatus", "A.mariannae", "complex"),
          covariate.labels = c("Bush", "Cactus", "Rock"))

glm_hapA <- glm(hap_A ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_hapA)
glm_hapB <- glm(hap_B ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_hapB)
glm_hapC <- glm(hap_C ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
summary(glm_hapC)

stargazer(glm_hapA, glm_hapB, glm_hapC, 
          type = "latex", 
          out = "models_summary_part1.html",
          title = "GLM-families-1",
          align = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("Hap A", "Hap B", "Hap C"),
          covariate.labels = c("Bush", "Cactus", "Rock"))

##
effemmatus <- allEffects(glm_A.emmatus)
plot(effemmatus)
effmariannae <- allEffects(glm_A.mariannae)
plot(effmariannae)
effhapC <- allEffects(glm_hapC)
plot(effhapC)




######Isolation-by-Distance & Isolation-by-Environment######
##IBD##
##Input of the distance matrices
# Genetic distance matrix
dist_matrix_genetic <- as.matrix(dxy_matrix)
dist_matrix_genetic_dist <- as.dist(dist_matrix_genetic)

##IBD##
#Input of the Distance matrix (reduced, since no Acrobeles emmatus was found in sampling spots G, K, and M)
dist_matrix_grid <- matrix(c(
  0, 1.8, 5.2, 6.3, 1.7, 2.4, 4.2, 2.3, 6.5, 2.9, 3.2, 4.2, 4.5,
  1.8, 0, 4.6, 4.8, 2.4, 3.9, 3.5, 1.2, 5.9, 1.9, 3.7, 5.2, 4.1,
  5.2, 4.6, 0, 2.3, 4.2, 7.3, 1.8, 2, 3.2, 3, 8.3, 5.3, 2.9,
  6.3, 4.8, 2.3, 0, 6.5, 9.6, 4, 5.1, 7.5, 3.4, 8.3, 4.2, 6.1,
  1.7, 2.4, 4.2, 6.5, 0, 1.9, 2.8, 3.5, 5.5, 2.9, 4.9, 5.1, 3.1,
  2.4, 3.9, 7.3, 9.6, 1.9, 0, 5.3, 6.7, 6.6, 4.8, 4.6, 8.4, 3.9,
  4.2, 3.5, 1.8, 4, 2.8, 5.3, 0, 4.1, 2, 2.1, 6.4, 5.2, 1.4,
  2.3, 1.2, 2, 5.1, 3.5, 6.7, 4.1, 0, 6.7, 2.5, 3, 2, 5.5,
  6.5, 5.9, 3.2, 7.5, 5.5, 6.6, 2, 6.7, 0, 4.8, 8.7, 7.7, 2.5,
  2.9, 1.9, 3, 3.4, 2.9, 4.8, 2.1, 2.5, 4.8, 0, 4.6, 4, 3.2,
  3.2, 3.7, 8.3, 8.3, 4.9, 4.6, 6.4, 3, 8.7, 4.6, 0, 4.7, 5.2,
  4.2, 5.2, 5.3, 4.2, 5.1, 8.4, 5.2, 2, 7.7, 4, 4.7, 0, 5.6,
  4.5, 4.1, 2.9, 6.1, 3.1, 3.9, 1.4, 5.5, 2.5, 3.2, 5.2, 5.6, 0
), nrow = 13, byrow = TRUE, dimnames = list(
  c("A", "B", "C", "D", "E", "F", "H", "I", "J", "L", "N", "O", "P"),
  c("A", "B", "C", "D", "E", "F", "H", "I", "J", "L", "N", "O", "P")
))
dist_matrix_grid_dist <- as.dist(dist_matrix_grid)

# Mantel test
mantel_test_IBD <- mantel(dist_matrix_genetic, dist_matrix_grid, method = "spearman", permutations = 999)

##IBB##
#Input of the barrier dataframe
environment_data <- data.frame(
  Sampling.site = c("EPT0A", "EPT0B", "EPT0C", "EPT0D", "EPT0E", "EPT0F", "EPT0H", "EPT0I", "EPT0J", "EPT0L", "EPT0N", "EPT0O", "EPT0P"),
  Bush = c(0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1),
  Rock = c(0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0),
  Cactus = c(1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0)
)
environment_dist <- vegdist(environment_data[, -1], method = "jaccard")
mantel_test_IBE <-mantel(dist_matrix_genetic, environment_dist, method = "spearman", permutations = 999)

#Plot the results
genetic_dist_vector <- dist_matrix_genetic[upper.tri(dist_matrix_genetic)]
geo_dist_vector <- dist_matrix_grid[upper.tri(dist_matrix_grid)]
environment_dist_matrix <- as.matrix(environment_dist)
environment_dist_vector <- environment_dist_matrix[upper.tri(environment_dist_matrix)]
ibd_data <- data.frame(GeographicDistance = geo_dist_vector, GeneticDistance = genetic_dist_vector)
ibe_data <- data.frame(EnvironmentDistance = environment_dist_vector, GeneticDistance = genetic_dist_vector)
p_value_ibd <- mantel_test_IBD$signif
p_value_ibe <- mantel_test_IBE$signif

plot_ibd <- ggplot(ibd_data, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_smooth(method = "lm", color = "darkorange", fill = "grey70", se = TRUE) +
  geom_point(color = "#456992", size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", color = "darkorange", se = FALSE) +
  labs(title = "Isolation-by-Distance",
       x = "Geographical Distance",
       y = "Genetic Distance") +
  annotate("text", x = 8.8, y = 0.0224, label = paste("p-value:", format(p_value_ibd, digits = 3)), size = 5)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

plot_ibe <- ggplot(ibe_data, aes(x = EnvironmentDistance, y = GeneticDistance)) +
  geom_smooth(method = "lm", color = "darkorange", fill = "grey70", se = TRUE) +
  geom_point(color = "#456992", size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", color = "darkorange", se = FALSE) +
  labs(title = "Isolation-by-Environment",
       x = "Environmental Distance",
       y = "Genetic Distance") +
  annotate("text", x = 0.85, y = 0.0224, label = paste("p-value:", format(p_value_ibe, digits = 3)), size = 5) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

ggarrange(plot_ibd, plot_ibe,
          labels = c("a", "b"),
          ncol = 1, nrow = 2)

######UCE Coverage######
coverage <- read.table(file.choose(), header = FALSE, sep = "\t", stringsAsFactors = FALSE)

##statistics
range_coverage <- range(coverage$V3)
mean_coverage <- mean(coverage$V3)
sd_coverage <- sd(coverage$V3)
median_coverage <- median(coverage$V3)
iqr_coverage <- IQR(coverage$V3)

print(mean_coverage)
print(sd_coverage)
print(median_coverage)
print(iqr_coverage)
print(range_coverage)

# Create the violin plot (plotting all coverages at all positions for all UCEs in one violin plot)
coverage_plot <- ggplot(coverage, aes(x = "", y = V3)) +
  geom_jitter(size = 0.01, color = "grey40", alpha = 0.05, width = 1)+
  geom_violin(fill = "#457992", color = "black", alpha = 0.8) +
  geom_hline(yintercept = 8, linetype = "solid", color = "darkorange", size = 0.5)+
  scale_y_break(c(1430,7570), space = 1)+
  geom_hline(yintercept = 50, linetype = "solid", color = "darkorange", size = 0.5)+
  labs(title = "Coverage per position (all UCEs)", x = "Position", y = "Coverage") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))

####per UCE
results <- coverage %>%
  group_by(V1) %>%
  summarise(
    mean_coverage = mean(V3, na.rm = TRUE),
    median_coverage = median(V3, na.rm = TRUE),
    sd_coverage = sd(V3, na.rm = TRUE),
    iqr_coverage=IQR(V3, na.rm = TRUE)
  )
results_sorted_median <- results %>%
  arrange(median_coverage)
results_sorted_mean <- results %>%
  arrange(mean_coverage)

results_sorted_median$V1 <- factor(results_sorted_median$V1, levels = results_sorted_median$V1)
results_sorted_mean$V1 <- factor(results_sorted_mean$V1, levels = results_sorted_mean$V1)

#####plot_each
median <- ggplot(results_sorted_median, aes(x = as.factor(V1), y = median_coverage)) +
  coord_flip()+
  geom_bar(stat = "identity", fill = "#457992") +
  geom_errorbar(aes(ymin = median_coverage - iqr_coverage, ymax = median_coverage + iqr_coverage), width = 0.2) +
  labs(title = "Median Coverage per UCE with IQR", x = "UCE", y = "Median Coverage") +
  theme_minimal()+
  theme(
    axis.text.y = element_blank())

mean <- ggplot(results_sorted_median, aes(x = as.factor(V1), y = mean_coverage)) +
  coord_flip()+
  geom_bar(stat = "identity", fill = "#457992") +
  geom_errorbar(aes(ymin = mean_coverage - sd_coverage, ymax = mean_coverage + sd_coverage), width = 0.2) +
  labs(title = "Mean Coverage per UCE with SD", x = "UCE", y = "Mean Coverage") +
  theme_minimal()+
  theme(
    axis.text.y = element_blank())
png("coverage_plot.png", width = 12, height = 8, units = "in", res = 400)
ggarrange(print(coverage_plot), print(median), print(mean), labels = c("a", "b", "c"), nrow = 1, ncol = 3)
dev.off()