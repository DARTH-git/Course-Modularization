probs <- c(0.714,0.0619,.2241, 0, 0.5728, 0.4272,0,0,1)

m_P <- matrix(probs, byrow=T, nrow =3)
m_e_P <- eigen(m_P)
m_P_2 <- m_e_P$vectors %*% diag(m_e_P$value )^(1/12) %*% solve(m_e_P$vectors)

m_M <- matrix(0, nrow = 14,ncol =3 )

r_P <- 1-exp(-m_P/12)
m_P_m <- -log(1-r_P)
for (i in 1:3){
  m_P_m[i,i] <- 1- (sum(m_P_m[i,]) - m_P_m[i,i]) 
}
m_M[1,1] <- 1
for (t in 1:12){
   m_M[t + 1,] <- m_M[t,] %*% (m_P_2)
}



