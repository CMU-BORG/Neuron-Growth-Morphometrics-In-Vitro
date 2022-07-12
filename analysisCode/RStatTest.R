
#####
# This performs the statistical analysis based on
# A S Liao, W Cui, V Webster-Wood, Y J Zhang
# Semi-automated quantitative evaluation of neuron developmental morphology in vitro using the change-point test
# Submitted to: Neuroinformatics 2022

# This script takes the morphometric results to conduct the Anderson-Darling, Kruskal-Wallis, and Dunn (with Bonferroni correction) tests as described in Liao et al. (2022):
# Morphometrics Analyzed:
# -total length per cell (sum of all neurite lengths)
# -average tortuosity per cell
# -number of end points per cell (degree)
# -number of neurites per cell
# -number of change points per cell
# -average segment length (distance between change points) per cell
# -average relative turning angle between change points per cell

#####
# Clear plots
if(!is.null(dev.list())) dev.off()

# Remove all existing R objects
rm(list=ls())

#call relevant libraries
library(dplyr)
library(nortest)
library(dunn.test)

#read the csv with the morphometrics
data <- read.csv("G:/Shared drives/CMU BORG - Neuron Culture and Quantification/NeuronQuantPaper_PublishCode/Publishable Code/Jupyter/allFeatures_fig7/perCellMetrics(den10,rem8,remMiss).csv")


#Summary Statistics
#http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r#multiple-pairwise-comparison-between-groups
#https://dplyr.tidyverse.org/reference/summarise.html

my_data <- data
my_data$div <- ordered(my_data$div, levels = c(0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0))
basicStats_degree <- group_by(my_data,div) %>%
                      summarise(
                        count = n(),
                        mean = mean(degree, na.rm = TRUE),
                        sd = sd(degree, na.rm = TRUE),
                        median = median(degree, na.rm = TRUE),
                        IQR = IQR(degree, na.rm = TRUE),
                        Q1 = quantile(degree,1/4),
                        Q3 = quantile(degree,3/4)
                      )

basicStats_numNeurites <- group_by(my_data,div) %>%
                            summarise(
                              count = n(),
                              mean = mean(numNeurites, na.rm = TRUE),
                              sd = sd(numNeurites, na.rm = TRUE),
                              median = median(numNeurites, na.rm = TRUE),
                              IQR = IQR(numNeurites, na.rm = TRUE),
                              Q1 = quantile(numNeurites,1/4),
                              Q3 = quantile(numNeurites,3/4)
                            )

basicStats_totCPs <- group_by(my_data,div) %>%
                      summarise(
                        count = n(),
                        mean = mean(totCPs, na.rm = TRUE),
                        sd = sd(totCPs, na.rm = TRUE),
                        median = median(totCPs, na.rm = TRUE),
                        IQR = IQR(totCPs, na.rm = TRUE),
                        Q1 = quantile(totCPs,1/4),
                        Q3 = quantile(totCPs,3/4)
                      )

basicStats_avgAbsRelTurnAngle <- group_by(my_data,div) %>%
                                  summarise(
                                    count = n(),
                                    mean = mean(avgAbsRelTurnAngle, na.rm = TRUE),
                                    sd = sd(avgAbsRelTurnAngle, na.rm = TRUE),
                                    median = median(avgAbsRelTurnAngle, na.rm = TRUE),
                                    IQR = IQR(avgAbsRelTurnAngle, na.rm = TRUE),
                                    Q1 = quantile(avgAbsRelTurnAngle,1/4),
                                    Q3 = quantile(avgAbsRelTurnAngle,3/4)
                                  )

basicStats_avgSegLen <- group_by(my_data,div) %>%
                          summarise(
                            count = n(),
                            mean = mean(avgSegLen, na.rm = TRUE),
                            sd = sd(avgSegLen, na.rm = TRUE),
                            median = median(avgSegLen, na.rm = TRUE),
                            IQR = IQR(avgSegLen, na.rm = TRUE),
                            Q1 = quantile(avgSegLen,1/4),
                            Q3 = quantile(avgSegLen,3/4)
                          )

basicStats_totLen <- group_by(my_data,div) %>%
                      summarise(
                        count = n(),
                        mean = mean(totLen, na.rm = TRUE),
                        sd = sd(totLen, na.rm = TRUE),
                        median = median(totLen, na.rm = TRUE),
                        IQR = IQR(totLen, na.rm = TRUE),
                        Q1 = quantile(totLen,1/4),
                        Q3 = quantile(totLen,3/4)
                      )

basicStats_avgTort <- group_by(my_data,div) %>%
                        summarise(
                          count = n(),
                          mean = mean(avgTort, na.rm = TRUE),
                          sd = sd(avgTort, na.rm = TRUE),
                          median = median(avgTort, na.rm = TRUE),
                          IQR = IQR(avgTort, na.rm = TRUE),
                          Q1 = quantile(avgTort,1/4),
                          Q3 = quantile(avgTort,3/4)
                        )

#Anderson-Darling Test for Normality
#https://www.r-bloggers.com/2021/11/anderson-darling-test-in-r-quick-normality-check/

adTest_avgAbsRelTurnAngle005 <- ad.test(data$avgAbsRelTurnAngle_0.5)
adTest_avgAbsRelTurnAngle010 <- ad.test(data$avgAbsRelTurnAngle_1.0)
adTest_avgAbsRelTurnAngle015 <- ad.test(data$avgAbsRelTurnAngle_1.5)
adTest_avgAbsRelTurnAngle020 <- ad.test(data$avgAbsRelTurnAngle_2.0)
adTest_avgAbsRelTurnAngle030 <- ad.test(data$avgAbsRelTurnAngle_3.0)
adTest_avgAbsRelTurnAngle040 <- ad.test(data$avgAbsRelTurnAngle_4.0)
adTest_avgAbsRelTurnAngle060 <- ad.test(data$avgAbsRelTurnAngle_6.0)

adTest_degree005 <- ad.test(data$degree_0.5)
adTest_degree010 <- ad.test(data$degree_1.0)
adTest_degree015 <- ad.test(data$degree_1.5)
adTest_degree020 <- ad.test(data$degree_2.0)
adTest_degree030 <- ad.test(data$degree_3.0)
adTest_degree040 <- ad.test(data$degree_4.0)
adTest_degree060 <- ad.test(data$degree_6.0)

adTest_numNeurites005 <- ad.test(data$numNeurites_0.5)
adTest_numNeurites010 <- ad.test(data$numNeurites_1.0)
adTest_numNeurites015 <- ad.test(data$numNeurites_1.5)
adTest_numNeurites020 <- ad.test(data$numNeurites_2.0)
adTest_numNeurites030 <- ad.test(data$numNeurites_3.0)
adTest_numNeurites040 <- ad.test(data$numNeurites_4.0)
adTest_numNeurites060 <- ad.test(data$numNeurites_6.0)

adTest_totCPs005 <- ad.test(data$totCPs_0.5)
adTest_totCPs010 <- ad.test(data$totCPs_1.0)
adTest_totCPs015 <- ad.test(data$totCPs_1.5)
adTest_totCPs020 <- ad.test(data$totCPs_2.0)
adTest_totCPs030 <- ad.test(data$totCPs_3.0)
adTest_totCPs040 <- ad.test(data$totCPs_4.0)
adTest_totCPs060 <- ad.test(data$totCPs_6.0)

adTest_avgSegLen005 <- ad.test(data$avgSegLen_0.5)
adTest_avgSegLen010 <- ad.test(data$avgSegLen_1.0)
adTest_avgSegLen015 <- ad.test(data$avgSegLen_1.5)
adTest_avgSegLen020 <- ad.test(data$avgSegLen_2.0)
adTest_avgSegLen030 <- ad.test(data$avgSegLen_3.0)
adTest_avgSegLen040 <- ad.test(data$avgSegLen_4.0)
adTest_avgSegLen060 <- ad.test(data$avgSegLen_6.0)

adTest_totLen005 <- ad.test(data$totLen_0.5)
adTest_totLen010 <- ad.test(data$totLen_1.0)
adTest_totLen015 <- ad.test(data$totLen_1.5)
adTest_totLen020 <- ad.test(data$totLen_2.0)
adTest_totLen030 <- ad.test(data$totLen_3.0)
adTest_totLen040 <- ad.test(data$totLen_4.0)
adTest_totLen060 <- ad.test(data$totLen_6.0)

adTest_avgTort005 <- ad.test(data$avgTort_0.5)
adTest_avgTort010 <- ad.test(data$avgTort_1.0)
adTest_avgTort015 <- ad.test(data$avgTort_1.5)
adTest_avgTort020 <- ad.test(data$avgTort_2.0)
adTest_avgTort030 <- ad.test(data$avgTort_3.0)
adTest_avgTort040 <- ad.test(data$avgTort_4.0)
adTest_avgTort060 <- ad.test(data$avgTort_6.0)

#anderson darling values are equivalent to Minitab results

#Kruskal-Wallis Test
#http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r#multiple-pairwise-comparison-between-groups
kw_avgAbsRelTurnAngle <- kruskal.test(avgAbsRelTurnAngle~div,data=my_data)
kw_avgSegLen <- kruskal.test(avgSegLen~div,data=my_data)
kw_totLen <- kruskal.test(totLen~div,data=my_data)
kw_avgTort <- kruskal.test(avgTort~div,data=my_data)

kw_degree <- kruskal.test(degree~div,data=my_data)
kw_numNeurites <- kruskal.test(numNeurites~div,data=my_data)
kw_totCPs <- kruskal.test(totCPs~div,data=my_data)

# Dunn Test with Bonferroni Correction
#https://stats.stackexchange.com/questions/355990/is-dunn-test-for-two-samples-equivalent-to-wilcox-test
#https://microbiome.github.io/tutorials/post_hoc.html
#https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
#"Compared to the Wilcoxon's test, the Dunn's test takes into account the ranking used by the Kruskal-Wallis test. It also does ties adjustments"
#https://stats.stackexchange.com/questions/126686/how-to-read-the-results-of-dunns-test
#https://www.programmingr.com/kruskal-wallis-rank-test/
#Actually, I think R's kruskal.test() automatically adjusts for ties. When running the KW test in Minitab and in R for degree, the chi-sq is the same for the 'Adjusted for Ties' values

dunn_avgAbsRelTurnAngle <- dunn.test(my_data$avgAbsRelTurnAngle, my_data$div,method = "bonferroni")
dunn_avgSegLen <- dunn.test(my_data$avgSegLen, my_data$div,method = "bonferroni")
dunn_totLen <- dunn.test(my_data$totLen, my_data$div,method = "bonferroni")
dunn_avgTort <- dunn.test(my_data$avgTort, my_data$div,method = "bonferroni")

dunn_degree <- dunn.test(my_data$degree, my_data$div,method = "bonferroni")
dunn_numNeurites <- dunn.test(my_data$numNeurites, my_data$div,method = "bonferroni")
dunn_totCPs <- dunn.test(my_data$totCPs, my_data$div,method = "bonferroni")