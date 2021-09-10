n <- 250
W <- replicate(100,runif(n, min = -1,  max = 1))
colnames(W) <- paste0("W", 1:100)
beta <- runif(10, -1, 1)/20
A <- rbinom(n, size = 1, prob = plogis(W[,10*(1:10)] %*% beta))

# CATE
Y <- rnorm(n, mean = A + W[,10*(1:10)] %*% beta, sd = 0.5)
data <- data.frame(W, A, Y)
formula = ~ 1  
output <-
  causalGLMnet(
    formula,
    data,
    W = colnames(W), A = "A", Y = "Y",
    estimand = "CATE" ,
    verbose = FALSE
  )
summary(output)


# OR
Y <- rbinom(n, size =  1, prob = plogis( A + W[,10*(1:10)] %&% beta))
data <- data.frame(W, A, Y)
formula  = ~ 1  
output <-
  causalGLMnet(
    formula,
    data,
    W = colnames(W), A = "A", Y = "Y",
    estimand = "OR" ,
    verbose = FALSE
  )
summary(output)

# RR
Y <- rpois(n, lambda = exp( A + W[,10*(1:10)] %*% beta))
data <- data.frame(W, A, Y)
formula  = ~ 1  
output <-
  causalGLMnet(
    formula,
    data,
    W = colnames(W), A = "A", Y = "Y",
    estimand = "RR" ,
    verbose = FALSE
  )
summary(output)
