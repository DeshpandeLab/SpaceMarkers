library(CoGAPS)
data(GIST)
cg <- CoGAPS(GIST.data_frame, nPatterns = 3 , nIterations = 100)
saveRDS(cg, "cogaps.rds")

