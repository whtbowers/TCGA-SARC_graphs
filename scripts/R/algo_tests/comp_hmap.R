library(ComplexHeatmap)
library(circlize)

set.seed(123)

mat <- cbind(rbind(matrix(rnorm(16, -1), 4), matrix(rnorm(32, 1), 8)),
             rbind(matrix(rnorm(24,1), 4), matrix(rnorm(48, -1), 8)))

# Permute rows and columns
mat <- mat[sample(nrow(mat), nrow(mat)), sample(ncol(mat), ncol(mat))]

rownames(mat) <- paste0("R", 1:12)
colnames(mat) <- paste0("c", 1:10)

# Initial basic heatmap
Heatmap(mat)

# Change colours
mat2 <- mat
mat2[1,1] = 100000
Heatmap(mat2, col = colorRamp2(c(-3,0,3), c("green", "white", "red")),
        cluster_rows = FALSE, cluster_columns = FALSE)

# SUPER GAY HEATMAP
Heatmap(mat, col = rev(rainbow(10)))

# If discrete values, colours should be specified as named vector
discrete_mat = matrix(sample(1:4, 100, replace = TRUE), 10, 10)
colors = structure(circlize::rand_color(4), names = c("1", "2", "3", "4"))
Heatmap(discrete_mat, col = colors)

# For character matrix (which surely exactly same as above almost?)
discrete_mat = matrix(sample(letters[1:4], 100, replace = TRUE), 10, 10)
colors = structure(circlize::rand_color(4), names = letters[1:4])
Heatmap(discrete_mat, col = colors)

# If NA values
mat_with_na = mat
mat_with_na[sample(c(TRUE, FALSE), nrow(mat)*ncol(mat), replace = TRUE, prob = c(1,9))] <- NA
Heatmap(mat_with_na, na_col = "orange", clustering_distance_rows = "pearson")

# 2 maps in one fig comparing color spaces
# Defaults LAB colourspace, but others available
# RGB, LAB, XYZ, sRGB, LUV
f1 = colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"))
f2 = colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
Heatmap(mat, col = f1, column_title = "LAB color space") +
  Heatmap(mat, col = f2, column_title = "RGB color space")

# Title scalebar
Heatmap(mat, name = "foo") # Simple way
Heatmap(mat, heatmap_legend_param = list(title = "legend")) # More complex, but more control w/ heatmap_legent_param

# Axis titles
Heatmap(
  mat,
  name = "foo",
  column_title = "I am a column title",
  row_title = "I am a row title"
)

# Bottom column title
Heatmap(
  mat,
  name = "foo",
  column_title = "I am a column title",
  column_title_side = "bottom"
)

# Bold top column title
Heatmap(mat, name = "foo", column_title = "I am a big column title", 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"))

# Rotated row title
Heatmap(mat, name = "foo", row_title = "row title", row_title_rot = 0)

## CLUSTERING
# Predefined methods euclidean or pearson

# Columns clustered only, column dendrogram hidden
Heatmap(
  mat, 
  name = "foo",
  cluster_rows = FALSE,
  show_column_dend = FALSE
  )

# Row dendrogram on right instead of left
Heatmap(
  mat, 
  name = "foo", 
  row_dend_side = "right"
  )

