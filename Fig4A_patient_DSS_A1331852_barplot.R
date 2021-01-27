
# Plot bar plot of A-1331852 DSS in patient ex vivo samples (Figure 4A)

# load libraries
library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(wesanderson)

# read patient DSS
dss <- read_excel("Patient_samples_DSS.xlsx", skip = 1) %>% dplyr::select(1:5)
colnames(dss) <- gsub("\\.\\.\\..", "", colnames(dss))
colnames(dss)[1] <- "Patient"

data <- dss
data$Subtype <- "Other AML"
data$Subtype[data$Patient%in%c("AML-3", "AML-4")] <- "AML with erythroid differentiation"  
data$Subtype[data$Patient%in%c("AML-1", "AML-2")] <- "AML with megakaryoblastic differentiation"
data$Subtype[data$Patient%in%c("AML-18", "AML-22")] <- "JAK2-mutated AML"
data$Subtype[data$Patient%in%c("AML-5")] <- "EPOR-mutated AML"
data$Subtype[data$Patient %in% c("Healthy-1", "Healthy-2")] <- "Healthy bone marrow"
data$Type <- "real"

data2 <- data
data2$`A-1331852` <- -1
data2$Type <- "fake"

data <- rbind(data, data2)

pal <- wes_palette("Darjeeling1")
#pal <- c("#ED2891", "#B2509E", "#D49DC7", "#C1A72F", "#E8C51D")

pdf("Fig4A_patient_DSS_A1331852_barplot.pdf", height = 3, width = 9)
ggplot(data, aes(x = reorder(Patient, -`A-1331852`), y = `A-1331852`, fill = Subtype, color = Type)) +
  geom_bar(stat = "identity") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0, face = "italic"),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("AML with erythroid differentiation" = pal[1],
                               "AML with megakaryoblastic differentiation" = pal[4],
                               "Other AML" = "grey70",
                               "JAK2-mutated AML" = pal[3],
                               "EPOR-mutated AML" = pal[2],
                              "Healthy bone marrow" = pal[5])) +
  scale_color_manual(values = c(NA, NA)) +
  ylab("DSS (A-1331852)") +
  labs(fill = "") +
  xlab("Patient") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(color = F)
dev.off()


