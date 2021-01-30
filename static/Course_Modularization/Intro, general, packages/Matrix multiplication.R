
# Matrix algebra
m_A <- m_B <- matrix(c(0.75, 0, 0, 0.2, 0.85, 0, 0.05, 0.15, 1), nrow = 3, ncol = 3)

## Matrix addition
m_A + m_B

## Matrix substration
m_A - m_B

## Element-wise Matrix multipicatcation

## matrix multipliaion 
m_A * m_B

## matrix devision - less relevant in decision making 
m_A / m_B  # note our example code has some zero values. Therefore you get NaN. 

m_D <- matrix(c(0.75, 0.07, 0.18, 0.10, 0.85, 0.05, 0.05, 0.15, 0.8), nrow = 3, ncol = 3)
m_E <- matrix(c(0.85, 0.05, 0.10, 0.04, 0.95, 0.01, 0.15, 0.35, 0.5), nrow = 3, ncol = 3)
m_D / m_E

## Matrix multipicatcation
m_A %*% m_B

v_A <- c(1, 0, 0) # make vector A

v_A %*% m_B



