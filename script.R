researchTable <- data.frame(read.table("/home/john/Database.csv", header=TRUE, sep=","))

#групи жени - 1, 2, 3, 4
groupTable <- table(researchTable$group)
groupPercentage <- round(prop.table(groupTable)*100, 2) 
pie(groupTable, labels = groupPercentage,  main = "Таргет трупи в %", col = rainbow(n = length(groupTable)))
groups <- unique(researchTable$group)
legend(x = 'bottomleft', legend = groups, cex = 0.8, fill = rainbow(length(groupTable)))

barplotRegion <- barplot(height = groupTable, col = "seagreen", main = "Брой жени във всяка група ", las = 1)
frequencies <- tabulate(researchTable$group)
text(x = barplotRegion, y = frequencies - 3, label = frequencies, pos = 3, cex = 1, col = "black")

getMode <- function(values) {
  uniqueValues <- unique(values)
  uniqueValues[which.max(tabulate(match(values,
                                        uniqueValues)))]
}
getDescriptiveWithHisto <- function(values, xlabArg, mainArg) {
  print(summary(values)) #min, max median, mean 
  print(var(values)) #дисперсия
  print(sd(values)) #стандартно отклонение
  print(getMode(values)) #мода
  h<-hist(values, breaks=10, col="red", xlab=xlabArg,
          main=mainArg)
  xfit<-seq(min(values),max(values),length=40)
  yfit<-dnorm(xfit,mean=mean(values),sd=sd(values))
  yfit <- yfit*diff(h$mids[1:2])*length(values)
  lines(xfit, yfit, col="blue", lwd=2)

  shapiro.test(values) 
  #From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. 
  #In other words, we can assume the normality.
}

groups.group1 <- researchTable[which(researchTable$group == 1), ]
groups.group2 <- researchTable[which(researchTable$group == 2), ]
groups.group3 <- researchTable[which(researchTable$group == 3), ]
groups.group4 <- researchTable[which(researchTable$group == 4), ]

group1CHOL <- groups.group1$CHOL.mmol.l 
group2CHOL <- groups.group2$CHOL.mmol.l
group3CHOL <- groups.group3$CHOL.mmol.l
group4CHOL <- groups.group4$CHOL.mmol.l

group1GLUC <- groups.group1$GLUC.mmol.l
group2GLUC <- groups.group2$GLUC.mmol.l
group3GLUC <- groups.group3$GLUC.mmol.l
group4GLUC <- groups.group4$GLUC.mmol.l

group1Tg <- groups.group1$Tg.mmol.l
group2Tg <- groups.group2$Tg.mmol.l
group3Tg <- groups.group3$Tg.mmol.l
group4Tg <- groups.group4$Tg.mmol.l

getDescriptiveWithHisto(group1CHOL, "CHOL", "CHOL group 1") # p-value = 0.02293 < 0.005 we assume abnormal distribution of the data
getDescriptiveWithHisto(group2CHOL, "CHOL", "CHOL group 2") # p-value = 0.2117 > 0.05 normal distribution
getDescriptiveWithHisto(group3CHOL, "CHOL", "CHOL group 3") # p-value = 0.344 > 0.05 normal distribution
getDescriptiveWithHisto(group4CHOL, "CHOL", "CHOL group 4") # 0.1965

getDescriptiveWithHisto(group1GLUC, "GLUC", "GLUC group 1") # 0.2761
getDescriptiveWithHisto(group2GLUC, "GLUC", "GLUC group 2") # 0.04414
getDescriptiveWithHisto(group3GLUC, "GLUC", "GLUC group 3") # 0.2576
getDescriptiveWithHisto(group4GLUC, "GLUC", "GLUC group 4") #0.0194

getDescriptiveWithHisto(group1Tg, "Tg", "Tg group 1") # 0.1082
getDescriptiveWithHisto(group2Tg, "Tg", "Tg group 2") # 0.01227
getDescriptiveWithHisto(group3Tg, "Tg", "Tg group 3") # 0.05004
getDescriptiveWithHisto(group4Tg, "Tg", "Tg group 4") # 0.050002



#group vs gluc
kruskal.test(GLUC.mmol.l ~ group, data = researchTable) # p-value is < 0.05 so we can conclude that there are significant differences between the treatment groups.
boxplot(researchTable$GLUC.mmol.l ~ researchTable$group, names = c("1", "2", "3", "4"), xlab = "Groups", ylab = "GLUC",
        main = "GLUC / Groups", col = rainbow(length(groupTable))) 
wilcox.test(groups.group1$GLUC.mmol.l, groups.group4$GLUC.mmol.l)
wilcox.test(groups.group2$GLUC.mmol.l, groups.group4$GLUC.mmol.l)
wilcox.test(groups.group3$GLUC.mmol.l, groups.group4$GLUC.mmol.l)

wilcox.test(groups.group1$CHOL.mmol.l, groups.group4$CHOL.mmol.l)

#group vs chol
kruskal.test(CHOL.mmol.l ~ group, data = researchTable) # p-value = 0.3941 > 0.05 so we can conclude that there are not any significant differences between the treatment groups.
boxplot(researchTable$CHOL.mmol.l ~ researchTable$group, names = c("1", "2", "3", "4"), xlab = "Groups", ylab = "CHOL",
        main = "CHOL / Groups", col = rainbow(length(groupTable))) 

#group vs tg
kruskal.test(Tg.mmol.l ~ group, data = researchTable) # p-value = 0.00001172 < 0.05 so we can conclude that there are significant differences between the treatment groups.
boxplot(researchTable$Tg.mmol.l ~ researchTable$group, names = c("1", "2", "3", "4"), xlab = "Groups", ylab = "Tg",
        main = "Tg / Groups", col = rainbow(length(groupTable))) 

#group 1
plot(groups.group1$GLUC.mmol.l, groups.group1$CHOL.mmol.l, main = "GLUC / CHOl Group 1", xlab = "GLUC", ylab = "CHOL")
cor(groups.group1$GLUC.mmol.l, groups.group1$CHOL.mmol.l, method = "spearman")

model <- lm(GLUC.mmol.l ~ CHOL.mmol.l, data = groups.group1)
summary(model)

#group 4
plot(groups.group4$GLUC.mmol.l, groups.group4$CHOL.mmol.l, main = "GLUC / CHOl Group 4", xlab = "GLUC", ylab = "CHOL")
cor(groups.group4$GLUC.mmol.l, groups.group4$CHOL.mmol.l, method = "spearman")

model <- lm(GLUC.mmol.l ~ CHOL.mmol.l, data = groups.group4)
summary(model)

#group 2
plot(groups.group2$GLUC.mmol.l, groups.group2$Tg.mmol.l, main = "GLUC / Tg Group 2", xlab = "GLUC", ylab = "Tg")
cor(groups.group2$GLUC.mmol.l, groups.group2$Tg.mmol.l, method = "spearman")

model <- lm(GLUC.mmol.l ~ Tg.mmol.l, data = groups.group2)
summary(model)

#group 1
plot(groups.group1$GLUC.mmol.l, groups.group1$Tg.mmol.l, main = "GLUC / Tg Group 1", xlab = "GLUC", ylab = "Tg")
cor(groups.group1$GLUC.mmol.l, groups.group1$Tg.mmol.l, method = "spearman")

model <- lm(GLUC.mmol.l ~ Tg.mmol.l, data = groups.group1)
summary(model)

#tg vs chol
#group 1
plot(groups.group1$CHOL.mmol.l, groups.group1$Tg.mmol.l, main = "CHOL / Tg Group 1", xlab = "GLUC", ylab = "Tg")
cor(groups.group1$CHOL.mmol.l, groups.group1$Tg.mmol.l, method = "spearman")

model <- lm(GLUC.mmol.l ~ Tg.mmol.l, data = groups.group1)
summary(model)

#group 2
plot(groups.group2$CHOL.mmol.l, groups.group2$Tg.mmol.l, main = "CHOL / Tg Group 2", xlab = "GLUC", ylab = "Tg")
cor(groups.group2$CHOL.mmol.l, groups.group2$Tg.mmol.l, method = "spearman")

model <- lm(CHOL.mmol.l ~ Tg.mmol.l, data = groups.group2)
summary(model)
