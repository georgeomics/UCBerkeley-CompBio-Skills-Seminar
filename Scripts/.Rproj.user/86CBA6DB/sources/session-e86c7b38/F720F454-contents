rm(list=ls())
library(tess3r)
library(maps)
library(rworldmap)
library(Rgraphviz)

# EXPLORING DATA
# grab data 
data(data.at)
genotype = data.at$X
coordinates = data.at$coord

# plot coordinates on a map
plot(coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

# RUNNING TESS
# 14 minutes for 25,000 SNPs, 100 individuals
tess3.obj <- tess3(X = genotype, 
                   coord = coordinates, 
                   K = 1:8, 
                   method = "projected.ls", 
                   ploidy = 1, 
                   openMP.core.num = 4) 

# CROSS-VALIDATION PLOTS
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# retrieve tess3 Q matrix for K = 5 clusters 
q.matrix <- qmatrix(tess3.obj, K = 5)

# STRUCTURE-LIKE BAR PLOTS
my.colors <- c("red","yellow","green","blue","violet")
my.palette <- CreatePalette(my.colors, 5)
barplot(q.matrix, border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 

# STRUCTURE-LIKE PIE CHARTS
plot(coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

# Loop through each row of q.matrix and coordinates
for (i in 1:nrow(q.matrix)) {
  # Call the pieGlyph function for each row with respective coordinates
  pieGlyph(q.matrix[i, ], coordinates[i, 1], coordinates[i, 2], 
           edges = 200, radius = 0.5, density = NULL, angle = 45, 
           col = my.colors, border = NULL, lty = NULL, 
           main = paste("Ancestry Coefficients - Point", i))
}


# INTERPOLATED MAP OF ANCESTRY COEFFICIENTS
plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4, 
     col.palette = my.palette)


# GENOME SCANS FOR SELECTION
# retrieve tess3 results for K = 5 
p.values <- pvalue(tess3.obj, K = 5)
hist(p.values, col = "lightblue") 

# Benjamini-Hochberg algorithm (controls false-discovery rate)
L = length(p.values)
fdr.level = 1e-4
w = which(sort(p.values) < fdr.level * (1:L)/L)
candidates = order(p.values)[w]
length(candidates)

# Manhattan plot (outliers highlighted in blue)
plot(p.values, main = "Manhattan plot", 
     xlab = "Locus id", 
     ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates, -log10(p.values)[candidates], 
       pch = 19, cex = .2, col = "blue")



