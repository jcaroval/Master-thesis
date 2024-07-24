######Packages#######
library(dplyr)
library(ggplot2) 
library(ggpubr) #for significance in boxplots
library(tidyr) #for manipulation of the dataset in Abundance EPT
library(igraph) #for grid (graph_from_data_frame)
library(vegan) #for metaMDS
library(corrplot) #for correlation plots
library(ggbreak) #for abundance plot (break axis)
library(effects)
library(stargazer)
library(abdiv)
library(gridExtra) #if using the second figure from alphadiv
library(pheatmap) #betadiv
library(reshape2)
library(ade4)

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

#the visualization was optimized using PowerPoint.
######Abundance Habitat ######
desired_order <- c("All.habitats", "Cactus&Bush", "Bush", "Bush&Rock", "Cactus", "Rock", "none") #sorted by the mean (highest to lowest)
diversity_data$Habitat <- factor(diversity_data$Habitat, levels = desired_order)
my_comparisons <- list (c("All.habitats","Bush"), c("All.habitats", "Bush&Rock"), c("All.habitats", "Cactus"), c("All.habitats", "Cactus&Bush"), c("All.habitats","none"), c("All.habitats","Rock"), c("Bush", "Bush&Rock"), c("Bush", "Cactus"), c("Bush", "Cactus&Bush"), c("Bush", "none"), c("Bush","Rock"), c("Bush&Rock", "Cactus"), c("Bush&Rock", "Cactus&Bush"), c("Bush&Rock","none"), c("Bush&Rock","Rock"), c("Cactus","Cactus&Bush"), c("Cactus","none"), c("Cactus","Rock"), c("Cactus&Bush", "none"), c("Cactus&Bush","Rock"), c("none","Rock"))
df_summary <- diversity_data %>%
  group_by(Habitat)%>%
  summarize(total.extract = mean(total.extract, na.rm = TRUE))
ggplot(data = df_summary, aes(x = Habitat, y = total.extract))+
  #ggtitle(label = "Nematode abundance per habitat")+
  geom_bar(stat = "identity", position = position_dodge(), fill = "#457992", alpha = 1)+
  #stat_summary(fun.data = mean_cl_boot, 
  #            geom = "pointrange",
  #           size = 2,
  #          shape = 4,
  #         alpha = 0.5)+
  geom_jitter(data=diversity_data, fill = "orange", 
              shape = 21, size = 5, alpha = 0.8, width = 0.1)+
  stat_summary(data=diversity_data, fun.data = mean_cl_boot, #mean_cl_boot instead of mean_cl_normal because with the bootstaps you can obtain cl's for population mean without assuming normality
               geom = "errorbar",
               aes(width = 0.4),
               size = 0.5,
               alpha = 1,
               color = "darkblue")+
  #geom_signif(comparisons=my_comparisons, map_signif_level = TRUE)+ #shows that all comparisons are insignificant. Therefore, only the left- and rightmost boxplots are be displayed in the final figure.
  geom_signif(comparisons=list(c("All.habitats","none")), map_signif_level = TRUE, y_position = 290, textsize = 5)+
  #stat_compare_means(data = diversity_data, comparisons = my_comparisons)
  #geom_signif(comparisons=list(c("Cactus","Bush&Rock")), map_signif_level = FALSE, y_position = 85)+ #could be displayed as well (<0.1)
  scale_y_continuous(breaks = c(seq(0,300, by = 25)))+
  #scale_y_continuous(breaks = c(seq(0, 100, by = 20), 115, seq(200, 350, by = 50)), 
  #                   limits = c(0, 330)) +
  #scale_y_break(c(110,200), scales = 0.3)+
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

#Display all the p-values:
ggplot(data = diversity_data, aes(x = Habitat, y = total.extract))+
  ggtitle(label = "Nematode abundance per habitat")+
  geom_boxplot(aes(fill = Habitat),
               linewidth = 0.2, alpha = 0.6, linetype = "blank")+
  stat_compare_means(comparisons = my_comparisons)

#Then add the p-values manually according to the comparisons
pvalues <- "All_habitats Cactus_Bush Bush Bush_Rock Cactus Rock none
All_habitats 1 1 0.8 0.67 0.33 0.67 1
Cactus_Bush 1 1 0.8 0.67 0.33 0.67 1
Bush 0.8 0.8 1 0.53 0.73 0.27 0.4
Bush_Rock 0.67 0.67 0.53 1 0.095 0.33 0.67
Cactus 0.33 0.33 0.73 0.095 1 0.43 0.33
Rock 0.67 0.67 0.27 0.33 0.43 1 0.67
none 1 1 0.4 0.67 0.33 0.67 1" 
pvalues <- read.table(text = pvalues, header = TRUE, row.names = 1)
print(pvalues)
######Diversity (Families & Genera)######
##alpha-div##
##Families
df_pool <- as.data.frame(totabundance_pool)

# Wandle die numerischen Spalten in eine Matrix um
mat_pool <- as.matrix(df_pool[, -1])  # Die erste Spalte 'Sampling.site' ist nicht relevant für die Diversitätsberechnung
richness_pool <- apply(mat_pool, 1, richness)
shannon_pool <- apply(mat_pool, 1, shannon)
invsimpson_pool <- apply(mat_pool, 1, invsimpson)

# Erstelle ein Dataframe für die Ergebnisse
div_data_pool <- data.frame(Sampling.site = LETTERS[1:16], Richness = richness_pool, InvSimpson = invsimpson_pool)

#####für Abbildung genutzt
# Shannon-Diversität plotten
fam <- ggplot(div_data_pool, aes(x = Sampling.site)) +
  #geom_line(aes(y = Shannon), color = "blue") +
  #geom_point(aes(y = Shannon), color = "blue") +
  #geom_line(aes(y = InvSimpson), color = "red") +
  geom_point(aes(y = InvSimpson), color = "#e31a1c", alpha = 0.7, shape = 16, size = 3) +
  #geom_line(aes(y = richness), color = "blue") +
  geom_point(aes(y = Richness), color = "#1f78b4", alpha = 0.7, shape = 16, size = 3) +
  labs(title = "Diversity Measures Family",
       x = "Sampling Spot",
       y = "Diversity Index",
       color = "Measure") +
  theme_minimal()+
  scale_color_manual(values = c("blue", "red"),
                     labels = c("Richness", "Inverse Simpson"))
##Genera
df_single <- as.data.frame(totabundance_single)

# Wandle die numerischen Spalten in eine Matrix um
mat_single <- as.matrix(df_single[, -1])  # Die erste Spalte 'Sampling.site' ist nicht relevant für die Diversitätsberechnung
richness_single <- apply(mat_single, 1, richness)
shannon_single <- apply(mat_single, 1, shannon)
invsimpson_single <- apply(mat_single, 1, invsimpson)
# Erstelle ein Dataframe für die Ergebnisse
div_data_single <- data.frame(Sampling.site = LETTERS[1:16], Richness = richness_single, InvSimpson = invsimpson_single)

#####für Abbildung genutzt
# Shannon-Diversität plotten
gen <- ggplot(div_data_single, aes(x = Sampling.site)) +
  #geom_line(aes(y = Shannon), color = "blue") +
  #geom_point(aes(y = Shannon), color = "blue") +
  #geom_line(aes(y = InvSimpson), color = "red") +
  geom_point(aes(y = InvSimpson), color = "#e31a1c", alpha = 0.7, shape = 16, size = 3) +
  #geom_line(aes(y = richness), color = "blue") +
  geom_point(aes(y = Richness), color = "#1f78b4", alpha = 0.7, shape = 16, size = 3) +
  labs(title = "Diversity Measures Genus",
       x = "Sampling Spot",
       y = "Diversity Index",
       color = "Measure") +
  theme_minimal()+
  scale_color_manual(values = c("#1f78b4", "#e31a1c"),
                     labels = c("Richness", "Inverse Simpson"))
ggarrange(fam, gen, nrow = 2, ncol = 1)

##beta-div##
##Families
increasing <- colorRampPalette(c("skyblue","beige","orange"))(300)

totabundance_pool <- diversity_data[,c(1,23:28)]
beta_pool <- as.data.frame(totabundance_pool)
beta_counts_pool <- beta_pool[,-1]
bray_curtis_pool <- vegdist(beta_counts_pool, method = "bray", binary = T)
jaccard_pool <- vegdist(beta_counts_pool, method = "jaccard", binary = T)
print(bray_curtis_pool)
print(jaccard_pool)

pheatmap(as.matrix(jaccard_pool),
         display_numbers = TRUE, #change to bray-curtis for second figure
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
bray_curtis_single <- vegdist(beta_counts_single, method = "bray", binary = T)
jaccard_single <- vegdist(beta_counts_single, method = "jaccard", binary = T)
print(bray_curtis_single)
print(jaccard_single)

pheatmap(as.matrix(jaccard_single),
         display_numbers = TRUE, #change to bray-curtis for second figure
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
#comparisons <- list(
#  c("EPT-0A", "EPT-0N"),
#  c("EPT-0C", "EPT-0N"),
#  c("EPT-0D", "EPT-0N"),
#  c("EPT-0F", "EPT-0N"),
#  c("EPT-0I", "EPT-0N"),
#  c("EPT-0O", "EPT-0N"))
ggplot(pi, aes(x = pop, y = avg_pi)) +
  geom_violin(fill = "grey", size = 0.2, alpha = 0.4, trim = FALSE, scale = "width") +
  #geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.2, size = 1, color = "#457992") +
  stat_summary(fun = mean, geom = "point", color = "#457992", size = 3, shape = 16) +
  #geom_point(aes(color = I("#457992")), width = 0.2, height = 0, size = 2, alpha = 0.4) +
  labs(x = "Sampling spot", y = "Average nucleotide diversity") +
  scale_y_continuous(breaks = c(seq (0, 0.26, by = 0.02)))+
  #  scale_y_break(c(0.15, 0.22))+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    panel.grid.major.y = element_line(size = 0.5, color = "grey80"),
    panel.grid.minor.y = element_line(size = 0.25, color = "grey90")
  )
### Significances: A-N *, C-N **, D-N **, F-N *, I-N *, N-O *

##### Nucleotide divergence between populations (dxy)
dxy <- read.table(file.choose(), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

dxy_summary <- dxy %>%
  group_by(pop1, pop2) %>%
  summarise(avg_dxy = mean(avg_dxy, na.rm = TRUE), .groups = 'drop')

# Berechne auch den umgekehrten Vergleich (pop2, pop1), um sicherzustellen, dass die Matrix symmetrisch ist
dxy_summary_reversed <- dxy_summary %>%
  rename(pop1 = pop2, pop2 = pop1)

# Kombiniere die beiden Datensätze
dxy_summary_combined <- bind_rows(dxy_summary, dxy_summary_reversed)

# Erstelle eine Matrix aus den Daten
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
##Comparison of the nucleotide diversities between different approaches##
diversity_Laura <- c(0.0402963,
                     0.0408627,
                     0.106631,
                     0.206227,
                     0.138237,
                     0.221082,
                     0.107891,
                     0.0545159,
                     0.0167081,
                     0.0488625,
                     0.00241247)
diversity_Julia <- pi[, 5]

dataframe <- factor(c(rep("A. emmatus (n=2093)", length(diversity_Julia)),
                      rep("Panagrolaimidae (n=11)", length(diversity_Laura))))
diversity <- c(diversity_Julia, diversity_Laura)
data <- data.frame(species = dataframe, diversity = diversity)

summary_stats <- data %>%
  group_by(species) %>%
  summarize(mean = mean(diversity),
            sd = sd(diversity),
            .groups = 'drop')
summary_stats

# Durchführung eines t-Tests, um die Signifikanz zu überprüfen
t_test_result <- t.test(diversity_Julia, diversity_Laura)

# P-Wert aus dem Test extrahieren
p_value <- t_test_result$p.value
p_value

ggplot(data, aes(x= species, y= diversity))+
  #geom_jitter(alpha = 0.4, width = 0.2, color = "grey")+
  geom_boxplot(fill="#457992", color="black", alpha = 0.7)+
  #geom_point(data = summary_stats, aes (x = species, y = mean), shape = 4, color = "orange", size = 4)+
  geom_signif(comparisons=list(c("A. emmatus (n=2093)","Panagrolaimidae (n=11)")), annotations = paste0("**"), y_position = 0.226, textsize = 7)+
  labs(title = paste("p-value:", format(p_value, digits = 3)),
       x = "Species",
       y = "Average nucleotide diversity")+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )##
######Correlation (Dependence on Habitat and Nematodes)######
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

corrplot(M, method = "circle", title = "Correlation of nematodes at the EPT", mar = c(0,0,2,0), diag = FALSE, outline = TRUE, bg = "transparent", 
         addgrid.col = "transparent", type = "lower", tl.col="black", tl.srt=90,
         p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.col = "white", pch.cex = 2, col = COL2("RdBu",200))
######NMDS######
#Subset the dataframe on which to base the ordination (dataframe 1)
data_single <- diversity_data[,c(14:19)]
data_pool <- diversity_data[,c(23:28)]
#Identify the columns that contain descriptive/environmental data (dataframe 2)
data_2 <- diversity_data[1:5]

NMDS_single <- metaMDS(data_single, distance = "bray", k = 2, trace=FALSE, trymax = 100)
NMDS_pool <- metaMDS(data_pool, distance = "bray", k = 2, trace=FALSE, trymax = 100)
stressplot(NMDS_single)
stressplot(NMDS_pool)

#Visualization of the NMDS
x_limits_single <- c(-1.1, 1.5)
y_limits_single <- c(-1.6, 0.5)
x_limits_pool <- c(-1.5, 1.5)
y_limits_pool <- c(-1.6, 1)
stress_value_single <- NMDS_single$stress
stress_value_pool <- NMDS_pool$stress
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
     cex.lab = 1.5,   # Vergrößert die Schriftgröße der Achsenbeschriftungen
     cex.axis = 1.3,  # Vergrößert die Schriftgröße der Achsenticks
     cex.main = 1.5,  # Vergrößert die Schriftgröße des Haupttitels
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
     cex.lab = 1.5,   # Vergrößert die Schriftgröße der Achsenbeschriftungen
     cex.axis = 1.3,  # Vergrößert die Schriftgröße der Achsenticks
     cex.main = 1.5,  # Vergrößert die Schriftgröße des Haupttitels
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
#The plots were visually further optimized using Inkscape. Namely, the "transparent" hull layers of the absent data were removed completely from the plot. 
#This was done because even when set to "transparent", they remained visible and thereby interfered with the data visualization (presence data only). 
#The removal of the absent data hulls did not change anything to the remaining plot.

##Statistical Analysis##
##Genera##
fit_bray_combined <- adonis2(data_single ~ Bush + Rock + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)
fit_bray_combined <- adonis2(data_single ~ Bush + Cactus + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*), Rock significant (*)
fit_bray_combined <- adonis2(data_single ~ Rock + Bush + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)
fit_bray_combined <- adonis2(data_single ~ Rock + Cactus + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)
fit_bray_combined <- adonis2(data_single ~ Cactus + Rock + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (**)
fit_bray_combined <- adonis2(data_single ~ Cactus + Bush + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)
##Families##
fit_bray_combined <- adonis2(data_pool ~ Bush + Rock + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)
fit_bray_combined <- adonis2(data_pool ~ Bush + Cactus + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)
fit_bray_combined <- adonis2(data_pool ~ Rock + Bush + Cactus, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)
fit_bray_combined <- adonis2(data_pool ~ Rock + Cactus + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (**)
fit_bray_combined <- adonis2(data_pool ~ Cactus + Rock + Bush, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (**)
fit_bray_combined <- adonis2(data_pool ~ Cactus + Bush + Rock, data=data_2, permutations = 999, method = "bray")
fit_bray_combined #Bush significant (*)

######GLMs######
##Families##
glm_Alaimidae <- glm(Alaimidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Anguinidae <- glm(Anguinidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Aphelenchoididae <- glm(Aphelenchoididae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Cephalobidae <- glm(Cephalobidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
glm_Panagrolaimidae <- glm(Panagrolaimidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Qudsianematidae <- glm(Qudsianematidae.count ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)

##Genera##
glm_Acrobeles <- glm(Acrobeles_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data)
glm_Acrobeloides <- glm(Acrobeloides_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Aphelenchoides <- glm(Aphelenchoides_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Stegelletina <- glm(Stegelletina_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Ditylenchus <- glm(Ditylenchus_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden
glm_Panagrolaimus <- glm(Panagrolaimus_single ~ Bush + Cactus + Rock, family = poisson(link = "log"), data = diversity_data) #kann nicht geplotted werden

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
#the plots were then saves as SVG-files and were concatenated into one file using Inkscape.
effCephalobidae <- allEffects(glm_Cephalobidae)
effQudsianematidae <- allEffects(glm_Qudsianematidae)
effAcrobeles <- allEffects(glm_Acrobeles)
plot(effCephalobidae)
plot(effQudsianematidae)
plot(effAcrobeles)

######Isolation-by-Distance & Isolation-by-Barrier######
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
mantel_test_IBD <- mantel(dist_matrix_genetic, dist_matrix_grid, method = "pearson", permutations = 999)

##IBB##
#Input of the barrier dataframe
barrier_data <- data.frame(
  Sampling.site = c("EPT0A", "EPT0B", "EPT0C", "EPT0D", "EPT0E", "EPT0F", "EPT0H", "EPT0I", "EPT0J", "EPT0L", "EPT0N", "EPT0O", "EPT0P"),
  Bush = c(0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1),
  Rock = c(0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0),
  Cactus = c(1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0)
)
barrier_dist <- vegdist(barrier_data[, -1], method = "jaccard")
mantel_test_IBB <-mantel(dist_matrix_genetic, barrier_dist, method = "pearson", permutations = 999)

#Plot the results
genetic_dist_vector <- dist_matrix_genetic[upper.tri(dist_matrix_genetic)]
geo_dist_vector <- dist_matrix_grid[upper.tri(dist_matrix_grid)]
barrier_dist_vector <- barrier_dist_matrix[upper.tri(barrier_dist_matrix)]
ibd_data <- data.frame(GeographicDistance = geo_dist_vector, GeneticDistance = genetic_dist_vector)
ibb_data <- data.frame(BarrierDistance = barrier_dist_vector, GeneticDistance = genetic_dist_vector)
p_value_ibd <- mantel_test_IBD$signif
p_value_ibb <- mantel_test_IBB$signif
barrier_dist_matrix <- as.matrix(barrier_dist)

plot_ibd <- ggplot(ibd_data, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_point(color = "blue", size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
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

plot_ibb <- ggplot(ibb_data, aes(x = BarrierDistance, y = GeneticDistance)) +
  geom_point(color = "blue", size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Isolation-by-Barrier",
       x = "Barrier Distance",
       y = "Genetic Distance") +
  annotate("text", x = 0.85, y = 0.0224, label = paste("p-value:", format(p_value_ibb, digits = 3)), size = 5)+
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

ggarrange(plot_ibd, plot_ibb,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)
