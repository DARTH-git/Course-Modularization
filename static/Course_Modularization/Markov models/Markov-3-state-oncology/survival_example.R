scale        <- 222 
shape        <- 1.5
hr_C1        <- 0.9
n_t          <- 1000
cycle_length <- 28


S_BSC <- pweibull(seq(1, n_t, cycle_length), shape ,scale, lower.tail = F)
d_BSC <- dweibull(seq(1, n_t, cycle_length), shape ,scale)
h_BSC <- hweibull(seq(1, n_t, cycle_length), shape ,scale)
H_BSC <- Hweibull(seq(1, n_t, cycle_length), shape ,scale)

tp <- 1 - exp(- diff(H_BSC))

tp_C1A <- 1 - exp(- hr_C1 * diff(H_BSC) )



data1 <- cbind(rlnorm(100,1,2
                      ),1,1)
data2 <- cbind(rlnorm(100,2,2),1,2)
data <- as.data.frame(rbind(data1,data2))

colnames(data) <- c("time", "cens", "trt")

test.fit <- flexsurvreg(Surv(time,cens)~trt, data = data , dist = "lognormal")
summary(test.fit)
# Trans Prob
# convert into rate
# apply hazard
# convert back 





scalePH    <- scale ^ (-shape)
scalePH_C1 <- scalePH * hr_C1


S_C1PH <- pweibullPH(seq(1, n_t, cycle_length), shape ,scalePH_C1, lower.tail = F)
d_C1PH <- dweibullPH(seq(1, n_t, cycle_length), shape ,scalePH_C1)
h_C1PH <- hweibullPH(seq(1, n_t, cycle_length), shape ,scalePH_C1)
H_C1PH <- HweibullPH(seq(1, n_t, cycle_length), shape ,scalePH_C1)


tp_C1B <- 1 - exp(- diff(H_BSCPH))

plot( tp_C1B, type = 'l')
lines(tp_C1A, col = 2)

