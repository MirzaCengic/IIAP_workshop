library(spocc)

bbox_peru <- c(-82.59846, -67.38535, -20.18324, 1.793963)


out <- occ(query = "Bambusa vulgaris", from = "gbif", #geometry = bbox_peru,
           has_coords = TRUE)
# Figure out a way how to expalain variable types
out

class(out)

str(out)


