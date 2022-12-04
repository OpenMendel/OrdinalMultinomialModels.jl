library(MASS)

# based on the example from "Moderna and applied statistics with R", page 204


# A proportional-odds model
house.plr.logistic <- polr(Sat ~ Infl + Type + Cont, data = housing, weights = Freq, method='logistic')
write.csv(fitted(house.plr.logistic),'fitted.logistic.csv',row.names=F,col.names=F)

house.plr.probit <- polr(Sat ~ Infl + Type + Cont, data = housing, weights = Freq, method='probit')
write.csv(fitted(house.plr.probit),'fitted.probit.csv',row.names=F,col.names=F)

# complimentary log-logistic
house.plr.cloglog <- polr(Sat ~ Infl + Type + Cont, data = housing, weights = Freq, method='cloglog')
write.csv(fitted(house.plr.cloglog),'fitted.cloglog.csv',row.names=F,col.names=F)

