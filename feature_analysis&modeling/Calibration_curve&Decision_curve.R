library(rms)
library(pROC)
library(rmda)
train <-read.csv("E:/Experiments/KangLi/COVID-19/run/analysis/data/nom_train.csv")
test <-read.csv("E:/Experiments/KangLi/COVID-19/run/analysis/data/nom_test.csv")

############
dd=datadist(train)
options(datadist="dd")
f1 <- lrm(y~ Rad_Score
          #+Age
          #+sex
          #+histological_subtype
          #+EGFR_mutation +targeted_drug
          ,data = train,x = TRUE,y = TRUE)

nom <- nomogram(f1, fun=plogis,fun.at=c(.001, .01, seq(.1,.9, by=.4)), lp=F, funlabel="Risk Index")
plot(nom)
#################
# ROC train
f2 <- glm(y~ Rad_Score
          #+Age
          #+histological_subtype
          ,data = train,family = "binomial")

pre <- predict(f2, type='response')
plot.roc(train$y, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, ci.type="bars", 
         of="thresholds",
         thresholds="best",
         print.thres="best",
         col="blue"
         #,identity=TRUE
         ,legacy.axes=TRUE,
         print.auc.x=ifelse(50,50),
         print.auc.y=ifelse(50,50)
)

#######################

# ROC test
pre1 <- predict(f2,newdata = test)
plot.roc(test$y, pre1,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, ci.type="bars",
         of="thresholds",
         thresholds="best",
         print.thres="best",
         col="blue",legacy.axes=TRUE,
         print.auc.x=ifelse(50,50),
         print.auc.y=ifelse(50,50)
)
#######################
## CC test
rocplot1 <- roc(test$y, pre)
#ci.auc(rocplot1)
f3 <- lrm(test$y ~ pre1,x = TRUE,y = TRUE)
cal <- calibrate(f3,  method = "boot", B = 1000)
plot(cal, xlab = "Nomogram Predicted Survival", ylab = "Actual Survival",main = "Calibration Curve")

#################
# CC train
rocplot <- roc(train$y, pre)
#ci.auc(rocplot1)
f4 <- lrm(train$y ~ pre,x = TRUE,y = TRUE)
cal <- calibrate(f4,  method = "boot", B = 1000)
plot(cal, xlab = "Nomogram Predicted Survival", ylab = "Actual Survival",main = "Calibration Curve")

##########################
# Decision Curve test
Score<- decision_curve(y~ 
                          Rad_Score, data = test, family = binomial(link ='logit'),
                          thresholds= seq(0,1, by = 0.01),
                          confidence.intervals =0.95,study.design = 'case-control',
                          population.prevalence = 0.9)

#Nomogram<- decision_curve(y~ Rad_Score,
#                             #+Age, 
#                             data = test,
#                             family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
#                             confidence.intervals= 0.95,study.design = 'case-control',
#                             population.prevalence= 0.9)

List<- list(Score)
plot_decision_curve(List,curve.names= c('Radiomics'),
                    cost.benefit.axis =FALSE,col = c('blue'),
                    confidence.intervals =FALSE,standardize = FALSE,
                    #legend.position = "none",
                    legend.position = "bottomleft"
)

######################
# Decision Curve train
Score<- decision_curve(y~ 
                         Rad_Score, data = train, family = binomial(link ='logit'),
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals =0.95,study.design = 'case-control',
                       population.prevalence = 0.9)

#Nomogram<- decision_curve(y~ Rad_Score,
#                             #+Age, 
#                             data = test,
#                             family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
#                             confidence.intervals= 0.95,study.design = 'case-control',
#                             population.prevalence= 0.9)

List<- list(Score)
plot_decision_curve(List,curve.names= c('Radiomics'),
                    cost.benefit.axis =FALSE,col = c('blue'),
                    confidence.intervals =FALSE,standardize = FALSE,
                    #legend.position = "none",
                    legend.position = "bottomleft"
)
