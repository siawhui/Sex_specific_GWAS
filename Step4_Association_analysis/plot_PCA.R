args <- commandArgs(trailingOnly = TRUE)

path <- args[1]
elbow <- args[2]

library(ggplot2)

eigenvals <- read.table(paste0(path, "/GWAS_PCA.eigenval"), header = FALSE, col.names = "eigenvalue")

eig_data <- data.frame(PC = 1:dim(eigenvals)[1], eigenvalue = eigenvals)

p <- ggplot(eig_data, aes(x = PC, y = eigenvalue)) +
  geom_line() + 
  geom_point() +
  xlab("Principal Component") +
  ylab("Eigenvalue") +
  ggtitle("Scree Plot") +
  theme_minimal() +
  geom_point(data = eig_data[eig_data$PC == elbow,], aes(x = PC, y = eigenvalue), 
             color = "red", size = 2)

ggsave(paste0(path,"/PCA_scree.png"), plot = p, width = 5, height = 4)
ggsave(paste0(path,"/PCA_scree.pdf"), plot = p, width = 5, height = 4)