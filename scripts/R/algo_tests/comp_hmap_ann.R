library(ComplexHeatmap)
library(circlize)

set.seed(123)

###############
# ANNOTATIONS #
###############

## Column annotations

df <- data.frame(type = c(rep("a", 5), rep("b", 5)))
ha <- HeatmapAnnotation(df = df)
ha # Display annotation information
ha <- HeatmapAnnotation(
  df = df, 
  col = list(type = c("a" = "red", "b" = "blue"))
  )

# If continuous, should be mapping function to colour
ha <- HeatmapAnnotation(
  df = data.frame(age = sample(1:202, 10)),
  col = list(age = colorRamp2(c(0,20), c("white", "red")))
)
ha
draw(ha, 1:10)

# Attach column annotations

df = data.frame(type = c(rep("a", 5), rep("b", 5)),
                age = sample(1:20, 10))

ha1 <- HeatmapAnnotation(
  df = df,
  col = list(
    type = c("a" = "red", "b" = "blue"),
    age = colorRamp2(c(0, 20), c("white", "red")))
)

ha2 <- HeatmapAnnotation(
  df = data.frame(age = sample(1:20, 10)),
  col = list(age = colorRamp2(c(0,20), c("white", "red")))
    )

mat = matrix(rnorm(80,2), 8, 10)
mat = rbind(mat, matrix(rnorm(40, -2), 4, 10))
rownames(mat) <- paste0("R", 1:12)
colnames(mat) <- paste0("C", 1:10)

Heatmap(
  mat, 
  top_annotation = ha1, 
  bottom_annotation = ha2
  )

## Row annotation
df = data.frame(type = c(rep("a", 6), rep("b", 6)))

ha = HeatmapAnnotation(
  df = df, 
  col = list(type = c("a" = "red", "b" = "blue")),
  which = "row", 
  width = unit(1, "cm")
  )

draw(ha, 1:12)

# Row annotation with heatmap
ha1 <- rowAnnotation(
  df = df, 
  col = list(type = c("a" = "red", "b" = "blue")), 
  width = unit(1, "cm")
  )

ha2 <- HeatmapAnnotation(
  df = df,
  col = list(
    type = c("a" = "red", "b" = "blue"),
    age = colorRamp2(c(0, 20), c("white", "red")))
)

structure(as.vector(type = as.character(coltab[,1])), cols = as.vector(as.character(coltab[,2])))

Heatmap(
  mat, 
  top_annotation = ha2
  ) + ha1
