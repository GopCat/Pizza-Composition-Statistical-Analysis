# Pizza Composition Classification Analysis
# Author: Aleksandra Ivanova
# Statistical modeling project

# Diskriminacni analyza
install.packages("MVN") # Q-Q plot pro normalitu
install.packages('heplots') #box test pro shodu kovariancnich matic
install.packages('klaR')
install.packages("ggpubr") #grafy
install.packages("MLmetrics")  # konfusni matice 
install.packages("pROC") # ROC krivka
install.packages("MASS")

# Logisticka regrese
install.packages("rcompanion") #pseudo R
install.packages("tidyverse")
install.packages("lmtest")
install.packages("caret")
install.packages("nnet")

# Analyza hlavnich komponent
install.packages("REdaS") #KMO statistiky
install.packages("ggfortify")
install.packages("factorextra")


# Import Datasets -> pizza.csv
PZZ <- pizza
# vylouceni nepotrebnych promennych = ID
PZZ <- PZZ[,-c(1)]
# oddeleni neoznacenych pizz
PZZ_NA <- PZZ[is.na(PZZ[,8]),]
PZZ <- PZZ[!is.na(PZZ[,8]),]

# Mahalanobis - odlehla pozorovani
PZZ_x <- PZZ[,-c(8)] # vylouceni promenne Type
means <- sapply(PZZ_x, mean, na.rm = TRUE)
mahal_pzz <- mahalanobis(PZZ_x, means, cov(PZZ_x), inverted = FALSE)      
summary(mahal_pzz)
plot(mahal_pzz)
plot(density(mahal_pzz, bw = 0.3), main="Squared Mahalanobis distances") ; rug(mahal_pzz)
alfa = 0.05
out_mah <- which(mahal_pzz > qchisq(1-alfa,df=ncol(PZZ_x)))
out_mah

# vylouceni odlehlych pozorovani
PZZ2 <- PZZ[-out_mah,]

# Multikolinearita
cor_PZZ2 <- cor(PZZ2)
library(corrplot)
corrplot(cor_PZZ2, "number")
library(matlib)
det(cor_PZZ2) # determinant matice R
cor_PZZ2_inv<-solve(cor_PZZ2) # inverzni matice R-1
VIF_PZZ2 <- diag(cor_PZZ2_inv) # VIF faktory
VIF_PZZ2 > 5	# VIF faktory vetsi nez 5 -> problematicke

# grafy dvou promennych (x = Class, y = lze menit jednotlive promenne)
plot(PZZ2$Type, PZZ2$voda)
plot(PZZ2$Type, PZZ2$bilkoviny)
plot(PZZ2$Type, PZZ2$tuk)
plot(PZZ2$Type, PZZ2$popeloviny)
plot(PZZ2$Type, PZZ2$sodik)
plot(PZZ2$Type, PZZ2$uhlovodik)
plot(PZZ2$Type, PZZ2$kalorie)

# jednofaktorova analyza rozptylu
PZZ2$Type <- as.factor(PZZ2$Type) # pokud by promenna nebyla kodovana jako faktor, zde se prekoduje
PZZ_anova_v <- aov(voda ~ Type, data=PZZ2)
summary(PZZ_anova_v) # vysledky Analyzy rozptylu - P-value by mela byt zamitnuta (H0 = shoda prumeru ve skupinach -> faktor by nemel vliv a typy vina se nelisi)
PZZ_anova_t <- aov(tuk ~ Type, data=PZZ2)
summary(PZZ_anova_t) #atd.

# test normality matice X
PZZ2_x <- PZZ2[,-c(8)]
library(mvnormtest)
mshapiro.test(t(PZZ2_x)) # statisticky test
library(MVN)
qqplot(qchisq(ppoints(1000), df = ncol(PZZ2_x)), mahal_pzz)
abline(0, 1, col = 'red')
# test normality nevysel, protoze je tam hodne odlehlych pozorovani - kdybychom vsechna odstranili, rozdeleni by bylo normalni

# jednorozmerna normalita
shapiro.test(PZZ2$voda)
shapiro.test(PZZ2$sodik)

# histogramy - normalita
hist(PZZ2$voda,probability=TRUE)
hist(PZZ2$sodik,probability=TRUE) #atd.

# grafy - normalita a korelace
library(klaR)
library(MASS)
library(car)
scatterplotMatrix(PZZ2[1:7])
# na diagonale - rozdeleni nah.velicin
# mimo diagonalu - zavislosti mezi promennymi = korelace by byt nemely, ale jsou



#### Analyza hlavnich komponent ####

PZZ2_all <- as.data.frame(rbind(PZZ2[,],PZZ_NA[,]))
PZZ2_allx <- PZZ2_all[,-c(1,5,7,8)]

# Kaiser_Meyer-Olkin statistika
library(REdaS)
KMOS(PZZ2_allx, use = "complete.obs")

# kovariancní a korelacni matice
cov_pzz<-cov(PZZ2_allx)
cor_pzz<-cor(PZZ2_allx)

# vlastni cisla
vlc_cov_pzz <- eigen(cov_pzz)
vlc_cor_pzz <- eigen(cor_pzz)
vlc_cor_pzz$values

# optimalni pocet hlavnich komponent
m_pzz=ncol(PZZ2_allx) # pocet sloupcu
m_pzz
sum(vlc_cov_pzz$values)/m_pzz	#minimalni "akceptovatelna" hodnota vlastniho cisla pro hlavni komponenty, v pripade analyzy zalozene na kovariancni matici
D_e_pzz = 1-(det(cor_pzz))^(1/m_pzz)  # efektivni rozptyl
k_pzz = m_pzz*(0.8-0.5*D_e_pzz) 
k_pzz # pocet vyznamnych hlavnich komponent

# analyza hlavnich komponent
pca_pzz <- prcomp(PZZ2_allx, center = TRUE,scale = TRUE) #scale=TRUE -> nejsou stejne merne jednotky, nutno znormovat data
summary(pca_pzz)
pca_pzz$rotation #komponentni zateze/korelace (loadings)
pzz_skore <- pca_pzz$x #komponentni skore
pzz_skore
screeplot(pca_pzz) 

# novy soubor s komponentami
PZZ_pca <- as.data.frame(cbind(pzz_skore[,1],PZZ2_all$voda,PZZ2_all$sodik,PZZ2_all$kalorie,PZZ2_all$Type))
library(data.table)
setnames(PZZ_pca, old=c("V1","V2","V3","V4","V5"), new=c("K1","voda","sodik","kalorie","Type"))
head(PZZ_pca)
sapply(PZZ_pca,class)

# rozdeleni souboru na predikcni a testovaci 80:20
PZZ_pca_NA <- PZZ_pca[is.na(PZZ_pca[,5]),]
PZZ_pca <- PZZ_pca[!is.na(PZZ_pca[,5]),]
library(caTools)
set.seed(123)
PZZ_pca_split = sample.split(PZZ_pca$Type,SplitRatio = 0.8)
PZZ_pca_test <- PZZ_pca[PZZ_pca_split, ]
PZZ_pca_pred <- PZZ_pca[!PZZ_pca_split, ]
table(PZZ_pca_test$Type) # kontrola, zda zustali rovnomerne kategorie


#### Diskriminacni analyza ####

## Boxuv test o shode kovariancnich matic (H0:kovariancni matice jsou shodne, H1: kovariancni matice jsou rozdilne)
library(car)
library(heplots)
PZZ_pca_box <- boxM(PZZ_pca_test[, 1:4], PZZ_pca_test$Type) #1:4 dle poctu porovnanych promennych
PZZ_pca_box # H0 zamitnuta, kovariance nejsou stejne (!!! nicmene tento test je citlivy na normalitu)
summary(PZZ_pca_box) # radove jsou vlastni cisla v radcich stejna

## Linearni diskriminacni analyza (LDA) na testovacim souboru
library(MASS)
PZZ_pca.lda.test <- lda(Type ~ ., data=PZZ_pca_test) #linearni disriminanci fukce kde vysvetlovanou promennou je Type
PZZ_pca.lda.test #vypsani vysledku (v datech mame 4 typy a tak pouzijeme 3 diskriminacni funkci)

# klasifikace jednotek v predikcnim souboru 
PZZ_pca.lda.pred <- predict(PZZ_pca.lda.test, PZZ_pca_pred)
PZZ_pca.lda.pred$class

# konfusni matice = kontingencni tabulka - porovnani zarazeni do kategorii (zjisteno x predikovano)
library(MLmetrics)
PZZ_pca_CM_lda <- ConfusionMatrix(PZZ_pca.lda.pred$class,PZZ_pca_pred$Type)
PZZ_pca_CM_lda

## Kvadraticka diskriminacni analyza (QDA) na testovacim souboru
PZZ_pca.qda.test <- qda(Type ~ ., data=PZZ_pca_test) #kvadraticka diskriminanci fce, kde vysvetlovanou promennou je Type
PZZ_pca.qda.test #vypsani vysledku 

# klasifikace jednotek v predikcnim souboru 
PZZ_pca.qda.pred <- predict(PZZ_pca.lda.test, PZZ_pca_pred)
PZZ_pca.qda.pred$class

# konfusni matice
PZZ_pca_CM_qda <- ConfusionMatrix(PZZ_pca.qda.pred$class,PZZ_pca_pred$Type)
PZZ_pca_CM_qda

# ROC krivka
library(pROC)
roccurve_l=roc(PZZ_pca_pred$Type ~ as.numeric(PZZ_pca.lda.pred$class))
plot(roccurve_l, print.auc = TRUE)
auc(roccurve_l)

roccurve_q=roc(PZZ_pca_pred$Type ~ as.numeric(PZZ_pca.qda.pred$class))
plot(roccurve_q, print.auc = TRUE)
auc(roccurve_q)


#### Logisticka regrese ####

# model logistické regrese
library(nnet)
PZZ_pca.lr.test <- multinom(as.factor(Type) ~ as.numeric(K1) + as.numeric(voda) + as.numeric(sodik) + as.numeric(kalorie), data = PZZ_pca_test)
summary(PZZ_pca.lr.test)

# predikce
PZZ_pca.lr.pred <- predict(PZZ_pca.lr.test, PZZ_pca_pred)  
PZZ_pca.lr.pred

# konfusni matice
library(MLmetrics)
PZZ_pca_CM_lr <- ConfusionMatrix(PZZ_pca.lr.pred,PZZ_pca_pred$Type)
PZZ_pca_CM_lr


#### Porovnani predikci na neznamych vzorcich ####

# Linearni diskriminacni analyza (LDA)
PZZ_NA.lda.pred <- predict(PZZ_pca.lda.test, PZZ_pca_NA)
PZZ_NA.lda.pred$class

# Kvadraticka diskriminacni analyza (QDA)
PZZ_NA.qda.pred <- predict(PZZ_pca.qda.test, PZZ_pca_NA)
PZZ_NA.qda.pred$class

# Logisticka regrese
PZZ_NA.lr.pred <- predict(PZZ_pca.lr.test, PZZ_pca_NA)
PZZ_NA.lr.pred

