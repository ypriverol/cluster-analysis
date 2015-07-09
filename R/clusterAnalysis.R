# cluster load
library(plyr)
library(seqinr)
library(ggplot2)
library(gridExtra)
library(corrplot)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

#Open the cluster result files
clusterAnalysis <- read.table(file = "/Users/yperez/work/ms_work/cluster-work/cluster_analysis_150407.contaminant.csv", sep ="\t", header = TRUE, comment.char = "")
save(clusterAnalysis, file="clusterAnalysis.rda")

fileDir = "report"

if (!file.exists(fileDir)){
    dir.create(file.path(fileDir))
}

png("report/distributionPeptides.png",width =  800)

highQuality = TRUE

# Get the number of clusters with one, two, three and four peptides.
numberPep <- nrow(clusterAnalysis)
numberSecond <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_SEQ_SECOND))
numberSecond <- nrow(numberSecond)
numberThird <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_SEQ_THIRD))
numberThird <- nrow(numberThird)
numberFour <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_SEQ_FOUR))
numberFour <- nrow(numberFour)
values <- c(numberPep, numberSecond, numberThird, numberFour)
colsVal <- c("Number of Clusters", "With Two Peptides","With Three Peptides","With four or more")
df = data.frame(values, colsVal)
colnames(df) <- c("values", "cols")

# plot
# Very basic bar graph
ggplot(data=df, aes(y=df$values, x = df$cols, fill=cols)) + geom_bar(stat = "identity") + ylab("Clusters Count") + xlab("Number of Peptides") + ggtitle("Number of clusters with one, two, three and four peptides") + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

dev.off()

# Check distribution of count and peptides

png("report/distributionPepCount.png",width =  800)

pepSEqCount <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_COUNT))
plot1 <- ggplot(pepSEqCount, aes(x=PEP_COUNT, fill=IS_PEP_CONTAMINANT)) + geom_bar() + ylab("Best Peptides Count") + xlab("Peptide Count") + ggtitle("Number of count for peptides for Best Hit") + scale_x_log10()

pepSEqSecondCount <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_COUNT_SECOND))
plot2 <- ggplot(pepSEqSecondCount, aes(x=PEP_COUNT_SECOND, fill=IS_PEP_CONTAMINANT_SECOND)) + geom_bar() + ylab("Second best Peptides Count") + xlab("Peptide Count") + ggtitle("Number of count for peptides for Second Best Hit") + scale_x_log10()

pepSEqThirdCount <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_COUNT_THIRD))
plot3 <- ggplot(pepSEqThirdCount, aes(x=PEP_COUNT_THIRD, fill=IS_PEP_CONTAMINANT_THIRD)) + geom_bar() + ylab("Third best Peptides Count") + xlab("Peptide Count") + ggtitle("Number of count for peptides for Third best Hit") + scale_x_log10()

pepSEqFourCount <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_COUNT_FOUR))
plot4 <- ggplot(pepSEqFourCount, aes(x=PEP_COUNT_FOUR, fill=IS_PEP_CONTAMINANT_FOUR)) + geom_bar() + ylab("Four best Peptides Count") + xlab("Peptide Count") + ggtitle("Number of count for peptides for Four Best Hit") + scale_x_log10()

multiplot(plot1, plot2, plot3, plot4, cols = 2)

dev.off()

png("report/peptideContamminant.png",width =  800)

#contaminant aanalysis
pepNotSecYes <- subset(clusterAnalysis, is.na(clusterAnalysis$PEP_CONTAMINANT))
pepNotSecYes <- subset(pepNotSecYes, !is.na(pepNotSecYes$PEP_CONTAMINANT_SECOND))
pepNotSecYes <- nrow(pepNotSecYes)

pepNotThirdYes <- subset(clusterAnalysis, is.na(clusterAnalysis$PEP_CONTAMINANT))
pepNotThirdYes <- subset(pepNotThirdYes, !is.na(pepNotThirdYes$PEP_CONTAMINANT_THIRD))
pepNotThirdYes <- nrow(pepNotThirdYes)

pepNotFourYes <- subset(clusterAnalysis, is.na(clusterAnalysis$PEP_CONTAMINANT))
pepNotFourYes <- subset(pepNotFourYes, !is.na(pepNotFourYes$PEP_CONTAMINANT_FOUR))
pepNotFourYes <- nrow(pepNotFourYes)

pepContaminat <- subset(clusterAnalysis, (clusterAnalysis$IS_PEP_CONTAMINANT == "true"))
pepContaminat <- subset(pepContaminat, !is.na(pepContaminat$PEP_COUNT_SECOND))
pepContaminat <- subset(pepContaminat, is.na(pepContaminat$PEP_CONTAMINANT_SECOND))
pepContaminat <- nrow(pepContaminat)

pepContaminat2 <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_CONTAMINANT))
pepContaminat2 <- subset(pepContaminat2, !is.na(pepContaminat2$PEP_COUNT_THIRD))
pepContaminat2 <- subset(pepContaminat2, is.na(pepContaminat2$PEP_CONTAMINANT_THIRD))
pepContaminat2 <- nrow(pepContaminat2)

pepContaminat3 <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_CONTAMINANT))
pepContaminat3 <- subset(pepContaminat3, !is.na(pepContaminat3$PEP_COUNT_FOUR))
pepContaminat3 <- subset(pepContaminat3, is.na(pepContaminat3$PEP_CONTAMINANT_FOUR))
pepContaminat3 <- nrow(pepContaminat3)

values <- c(pepNotSecYes, pepNotThirdYes, pepNotFourYes, pepContaminat, pepContaminat2, pepContaminat3)
colsVal <- c("Pep Not/Sec Yes", "Pep Not/Third Yes","Pep Not/Four Yes","Pep Yes/ Sec No", "Pep Yes/ Third No", "Pep Yes/ Four No")
df = data.frame(values, colsVal)
colnames(df) <- c("values", "cols")
ggplot(data=df, aes(y=df$values, x = df$cols, fill=cols)) + geom_bar(stat = "identity") + ylab("Clusters Count") + xlab("Contaminant Combination") + ggtitle("Contaminants vs Peptide Rank") + ylim(0,1000000) + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

dev.off()

png("report/similarity.png",width =  800)

#Analysis of similarity Scores
pepSecondSimilarity <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_SIMILARITY_SECOND))

plot1 <- ggplot(pepSecondSimilarity, aes(x=PEP_SIMILARITY_SECOND, fill=IS_PEP_CONTAMINANT_SECOND)) + geom_bar() + geom_vline(xintercept = 90) + ylab("Second best Peptides Count") + xlab("Aligmenent Score") + ggtitle("Similarity Score histogram for Second best Hit") + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

pepThirdSimilarity <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_SIMILARITY_THIRD))

plot2<- ggplot(pepThirdSimilarity, aes(x=PEP_SIMILARITY_THIRD, fill=IS_PEP_CONTAMINANT_THIRD)) + geom_bar() + geom_vline(xintercept = 90) + ylab("Third best Peptides Count") + xlab("Aligmenent Score") + ggtitle("Similarity Score histogram for Third best Hit") + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

pepFourSimilarity <- subset(clusterAnalysis, !is.na(clusterAnalysis$PEP_SIMILARITY_FOUR))

plot3 <- ggplot(pepFourSimilarity, aes(x=PEP_SIMILARITY_FOUR, fill=IS_PEP_CONTAMINANT_FOUR)) + geom_bar() + geom_vline(xintercept = 90) + ylab("Four best Peptides Count") + xlab("Aligmenent Score") + ggtitle("Similarity Score histogram for Four best Hit") + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

multiplot(plot1, plot2, plot3, cols = 2)

dev.off()

png("report/correlation.png",width =  800)

df <- clusterAnalysis[, !sapply(clusterAnalysis, is.factor)]
df <- subset(df, select=-c(NUM_PROJECTS,NUM_PROJECTS_HIGHEST))
datMy.scale<- scale(df,center=TRUE,scale=TRUE);
corMatMy <- cor(datMy.scale, use="pairwise.complete.obs")
#compute the correlation matrix
plot(corMatMy)
corrplot(corMatMy, order = "hclust")

dev.off()

pepSecondSimilarity <- subset(pepSecondSimilarity, pepSecondSimilarity$PEP_SIMILARITY_SECOND >= 90)
pepSecondSimilarity <- nrow(pepSecondSimilarity)

