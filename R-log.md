<h1>Aim</h1>



The commands were run on the R version 4.3.2.

<h1>Packages</h1>
<ul>• <tt>ade4</tt> | Dray and Dufour, 2007</ul>
<ul>• <tt>dplyr (V1.1.4)</tt> | Wickham et al., 2023</ul>
<ul>• <tt>ggplot2</tt> | Wickham, 2016</ul>
<ul>• <tt>ggpubr (V0.6.0)</tt> | Kassambara, 2023</ul>
<ul>• <tt>tidyr (V1.3.1)</tt> | Wickham et al., 2024</ul>
<ul>• <tt>igraph (V2.0.3)</tt> | Csardi and Neupusz, 2006</ul>
<ul>• <tt>vegan (V2.6-6.1</tt> | Oksanen et al., 2024</ul>
<ul>• <tt>corrplot (V0.92)</tt> | Wei and Simko, 2021</ul>
<ul>• <tt>ggbreak</tt> | Xu et al., 2021</ul>
<ul>• <tt>effects</tt> | Fox and Weisberg, 2019 and Fox, 2003</ul>
<ul>• <tt>stargazer (V5.2.3)</tt> | Hlavak, 2022</ul>
<ul>• <tt>abdiv (V0.2.0)</tt> | Bittinger, 2020</ul>
<ul>• <tt>gridExtra (V2.3)</tt> | Auguie, 2017</ul>
<ul>• <tt>pheatmap (V1.0.12)</tt> | Kolde, 2019</ul>
<ul>• <tt>reshape2</tt> | Wickham, 2007</ul>
<ul>• <tt></tt> | </ul>

<h1>1 Basic preparations</h1>
<h2>1.1 Input of the dataset</h2>

The habitat characteristics need to be transformed into factors

```
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
```

<h2>1.2 Definition of color palettes</h2>

```
col <- c("#999999","#CC79A7" , "#56B4E9","#6EED80",  "#F0E442","#E69F00", "#0072B2" )
```

<h1>2 Sampling Grid and Distance Matrix</h1>
<h2>2.1 Define coordinates of the sampling spot</h2>

This was done by a 100x100 raster in PowerPoint based on the scheme that was provided from the sampling site. A coordinate system was then created for the coordinates.

```
coord <- data.frame(
  Point = c("A", "B", "C",  "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"),
  X = c(6.1, 6.2, 7.8, 9.8, 5.6, 3.8, 3.0, 4.8, 8.8, 2.4, 1.4, 6.4, 7.2, 7.6, 10.0, 4.7),
  Y = c(3.8, 4.6, 8.7, 10.1, 5.1, 2.6, 4.6, 7.8, 7.0, 9.1, 8.1, 6.5, 2.8, 2.1, 8.2, 6.4))
```

<h2>2.2 Input of the distance matrix</h2>

```
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
```
<h2>2.3 Visualization of the incomplete distance matrix</h2>

```
pdf("incompl_distmatrix.pdf", width = 10, height = 8, pointsize = 12)
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
dev.off()
```
In addition, the distance matrix was exported as an svg file. It was then further optimized using InkScape.

<h2>2.4 Convert the given coordinates to an igraph object</h2>

```
graph <- graph_from_data_frame(as.data.frame(t(combn(coord$Point, 2))))
V(graph)$X <- coord$X
V(graph)$X <- coord$X
```
<h2>2.5 Predict missing distances based on Euclidian Distance</h2>

```
for (i in 1:length(coord$Point)) {
  for (j in 1:length(coord$Point)) {
    if (is.na(distance_matrix[i, j]) && i != j) {
      distance <- sqrt((coord$X[i] - coord$X[j])^2 + (coord$Y[i] - coord$Y[j])^2)
      distance_matrix[i, j] <- distance
      distance_matrix[j, i] <- distance
    }
  }
}
```

<h2>2.6 Visualize the completed distance matrix</h2>

```
pdf("compl_distmatrix.pdf", width = 10, height = 8, pointsize = 12)
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
dev.off()
```
In addition, the distance matrix was exported as an svg file. It was then further optimized using InkScape.

<h2>2.7 Export the completed distance matrix as a CSV file</h2>

```
write.csv(distance_matrix, file = "compl_dist_matrix.csv", row.names = TRUE)
```

<h2>2.8 Plot the sampling spots A-P according to the completed distance matrix</h2>

```
mds_result <- cmdscale(distance_matrix, k = 2) # k = 2 für 2D-Darstellung
coord$X <- mds_result[, 1]
coord$Y <- mds_result[, 2]
coord$X <- -coord$X 
coord$X <- coord$X+6
coord$Y <- coord$Y+5
```

<h2>2.9 Visualize the grid</h2>

```
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

ggsave("EPT_Sampling_Grid.jpg", plot = g, width = 8, height = 8, units = "in", dpi = 800)                         
```
<h1>3 Abundance</h1>

<h2>3.1 Total abundance</h2>

<h3>3.1.1 Create subsets of the dataset</h3>

```
totabundance <- diversity_data[,c(1,23:29)]
colnames(totabundance) <- c("Sampling.site", "Cephalobidae", "Panagrolaimidae", "Alaimidae", "Aphelenchoididae", "Qudsianematidae", "Anguinidae", "unidentified")

relabundance <- diversity_data[,c(1,31:37)]
colnames(relabundance) <- c("Sampling.site", "Cephalobidae", "Panagrolaimidae", "Alaimidae", "Aphelenchoididae", "Qudsianematidae", "Anguinidae", "unidentified")
```

<h3>3.1.2 Modify the subsets</h3>

```
totlong <- pivot_longer(totabundance, -Sampling.site, names_to = "Taxonomy", values_to = "count")
rellong <- pivot_longer(relabundance, -Sampling.site, names_to = "Taxonomy", values_to = "relabundance")

totlong$Sampling.site <- factor(totlong$Sampling.site, levels = c("EPT0I", "EPT0J", "EPT0F", "EPT0O", "EPT0D", "EPT0G", "EPT0E", "EPT0C", "EPT0A", "EPT0B", "EPT0H", "EPT0M", "EPT0P", "EPT0K", "EPT0N", "EPT0L"))
rellong$Sampling.site <- factor(rellong$Sampling.site, levels = c("EPT0I", "EPT0J", "EPT0F", "EPT0O", "EPT0D", "EPT0G", "EPT0E", "EPT0C", "EPT0A", "EPT0B", "EPT0H", "EPT0M", "EPT0P", "EPT0K", "EPT0N", "EPT0L"))

rellong$Taxonomy <- factor(rellong$Taxonomy, levels = c("unidentified", "Alaimidae", "Aphelenchoididae","Panagrolaimidae", "Anguinidae", "Qudsianematidae", "Cephalobidae"))
totlong$Taxonomy <- factor(totlong$Taxonomy, levels = c("unidentified", "Alaimidae", "Aphelenchoididae","Panagrolaimidae", "Anguinidae", "Qudsianematidae", "Cephalobidae"))
```

<h3>3.1.3 Visualize the total abundance, relative abundance, and habitat composition</h3>

```
col <- c("#999999","#CC79A7" , "#56B4E9","#6EED80",  "#F0E442","#E69F00", "#0072B2" )

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

ggsave("Abundance_habitat.jpg", plot = combined_habitat, width = 12, height = 10, units = "in", dpi = 800)
```

<h2>3.2 Abundance per habitat</h2>

```
desired_order <- c("All.habitats", "Cactus&Bush", "Bush", "Bush&Rock", "Cactus", "Rock", "none") #sorted by the mean (highest to lowest)
diversity_data$Habitat <- factor(diversity_data$Habitat, levels = desired_order)
my_comparisons <- list (c("All.habitats","Bush"), c("All.habitats", "Bush&Rock"), c("All.habitats", "Cactus"), c("All.habitats", "Cactus&Bush"), c("All.habitats","none"), c("All.habitats","Rock"), c("Bush", "Bush&Rock"), c("Bush", "Cactus"), c("Bush", "Cactus&Bush"), c("Bush", "none"), c("Bush","Rock"), c("Bush&Rock", "Cactus"), c("Bush&Rock", "Cactus&Bush"), c("Bush&Rock","none"), c("Bush&Rock","Rock"), c("Cactus","Cactus&Bush"), c("Cactus","none"), c("Cactus","Rock"), c("Cactus&Bush", "none"), c("Cactus&Bush","Rock"), c("none","Rock"))
my_colors <- c('All.habitats' = '#F0E442', 
               'Cactus&Bush' = '#6a3d9a', 
               'Bush' = '#E69F00', 
               'Bush&Rock' = '#0072B2', 
               'Cactus' = '#33a02c', 
               'Rock' = '#999999', 
               'none' = '#56B4E9')
p <- ggplot(data = diversity_data, aes(x = Habitat, y = total.extract))+
  ggtitle(label = "Nematode abundance per habitat")+
  geom_boxplot(aes(fill = Habitat),
               linewidth = 0.2, alpha = 0.7, linetype = "blank")+
  stat_summary(fun.data = mean_cl_boot, #mean_cl_boot instead of mean_cl_normal because with the bootstaps you can obtain cl's for population mean without assuming normality
               geom = "errorbar",
               aes(width = 0.4),
               size = 0.5,
               alpha = 0.5)+
  stat_summary(fun.data = mean_cl_boot, 
               geom = "pointrange",
               size = 1,
               shape = 4,
               alpha = 0.5)+
  geom_point(aes(fill = Habitat), 
             shape = 21, size = 3.5, alpha = 0.5)+
  #geom_signif(comparisons=my_comparisons, map_signif_level = TRUE)+ #shows that all comparisons are insignificant. Therefore, only the left- and rightmost boxplots are be displayed in the final figure.
  geom_signif(comparisons=list(c("All.habitats","none")), map_signif_level = TRUE, y_position = 300)+
  geom_signif(comparisons=list(c("Cactus","Bush&Rock")), map_signif_level = FALSE, y_position = 85)+ #could be displayed as well (<0.1)
  scale_fill_manual(values = my_colors)+
  scale_y_continuous(breaks = c(seq(0, 100, by = 20), 115, seq(200, 350, by = 50)), 
                     limits = c(0, 330)) +
  scale_y_break(c(110,200), scales = 0.3)+
  scale_x_discrete(labels = c("CactBushRock", "CactBush", "Bush", "BushRock", "Cact", "Rock", "none"))+
  ylab("Count [n]")+
  xlab("Habitat")+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 14, margin = margin(t=10)),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.15, color = "gray"),
        panel.grid.minor.y = element_line(linewidth = 0.05, color = "gray"))
ggsave(p, filename = "Nematode abundance boxplot.jpg", width = 10, height = 6, dpi=800)
p
```

<h1>4 Diversity</h1>

