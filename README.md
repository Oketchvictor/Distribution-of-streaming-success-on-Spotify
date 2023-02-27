# Distribution-of-streaming-success-on-Spotify
This is to help those who are learning data analysis on raw data using R software
a)Simulation study
```{r}
# Set parameters for the simulation
n_list <- c(25, 50, 100, 500) # Sample sizes
num_simulations <- 2000 # Number of simulations for each sample size
true_mean <- c(0, 0, 0, 0) # True mean vector
true_cov <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3) # True covariance matrix
true_df <- 5 # True degrees of freedom parameter
print(true_mean)
print(true_cov)
```

```{r}
set.seed(123)
# Simulation parameters
n_sim <- 2000
n_vec <- c(25, 50, 100, 500)
Sigma <- matrix(c(1,1/3,1/3,2), nrow=2)
Sigma
```

Simulation
```{r}
set.seed(123) # for reproducibility
n_vec <- c(25, 50, 100, 500)
n_sims <- 2000
mu <- c(1, 1)
Sigma <- matrix(c(1, 1/3, 1/3, 2), nrow = 2)
nu <- 4
theta_hat <- matrix(nrow = n_sims * length(n_vec), ncol = 3)
theta_se <- matrix(nrow = n_sims * length(n_vec), ncol = 3)
for (i in 1:length(n_vec)) {
  n <- n_vec[i]
  for (j in 1:n_sims) {
    X <- matrix(rt(n*2, df = nu), nrow = n)
    S <- cov(X)
    df <- n - 1
    Sigma_hat <- S * df / nu
    theta_hat[(i-1)*n_sims+j, ] <- c(Sigma_hat[1,1], Sigma_hat[1,2], Sigma_hat[2,2])
    se_Sigma_hat <- sqrt(2/n) * Sigma_hat / sqrt(df)
    theta_se[(i-1)*n_sims+j, ] <- c(se_Sigma_hat[1,1], se_Sigma_hat[1,2], se_Sigma_hat[2,2])
  }
}
print(Sigma_hat)
print(theta_hat[(i-1)*n_sims+j, ])
print(theta_se[(i-1)*n_sims+j, ])
print(se_Sigma_hat)
```

```{r}
set.seed(123)
n <- c(25, 50, 100, 500)  # sample sizes
B <- 2000  # number of datasets
mu <- c(1, 1)
Sigma <- matrix(c(1, 1/3, 1/3, 2), nrow = 2)
nu <- 4
# Inverse transform method
Z <- matrix(rnorm(max(n) * B * 2), nrow = max(n) * B, ncol = 2)
X <- matrix(rchisq(max(n) * B, df = nu), nrow = max(n) * B, ncol = 1)
# Reshape data into a list of datasets
data_list <- lapply(n, function(nn) array(X[1:(B*nn), ], dim = c(B, nn, 2)))
```

```{r}
library(mvtnorm)
# Function to compute the log-likelihood
loglik <- function(theta, data) {
  mu <- theta[1:2]
  Sigma <- matrix(theta[3:6], nrow = 2)
  nu <- exp(theta[7])
  d <- 2
  ll <- sum(log(dmvt(data, mu = mu, sigma = Sigma, df = nu)))
  ll <- ll - 0.5 * nu * sum(log(diag(chol(Sigma))))
  ll <- ll - d/2 * log(nu) - d/2 * log(2*pi)
  return(ll)
}
# Function to compute the gradient of the log-likelihood
grad_loglik <- function(theta, data) {
  mu <- theta[1:2]
  Sigma <- matrix(theta[3:6], nrow = 2)
  nu <- exp(theta[7])
  d <- ncol(data)
  # Compute the gradient of the log-likelihood
  g1 <- rep(0, 2)
  g2 <- matrix(0, nrow = 2, ncol = 2)
  g3 <- 0
  for (i in 1:nrow(data)) {
    x <- data[i, ]
    t1 <- dt(x, df = nu, ncp = 0, log = TRUE) - log(nu) - log(nu + t(x - mu) %*% solve(Sigma) %*% (x - mu)) - d/2 * log(2*pi)
    t2 <- (nu + d) / (nu + t(x - mu) %*% solve(Sigma) %*% (x - mu))
    t3 <- 0.5 * solve(Sigma) %*% (x - mu) %*% t(x - mu) %*% solve(Sigma)
    g1 <- g1 + t2 * t3 %*% (x - mu)
    g2 <- g2 + t2 * (t3 - solve(Sigma))
    g3 <- g3 + t1
  }
}
```


##
```{r}
library(mvtnorm) # for generating multivariate normal samples
set.seed(123) # for reproducibility
# parameters of the multivariate t-distribution
mu <- c(1,1)
Sigma <- matrix(c(1,1/3,1/3,2), nrow=2)
nu <- 4
# function to generate samples from multivariate t-distribution
generate_t_samples <- function(n) {
  z <- rmvnorm(n, mean=mu, sigma=Sigma)
  w <- rchisq(n, df=nu)
  x <- t(t(z)/sqrt(w/nu))
  return(x)
}
# generate 2000 samples of size 25, 50, 100, and 500
n_vals <- c(25,50,100,500)
num_sims <- 2000
t_samples <- lapply(n_vals, function(n) replicate(num_sims, generate_t_samples(n), simplify = FALSE))
```

```{r}
# negative log-likelihood function for multivariate t-distribution
neg_loglik_t <- function(x, mu, Sigma, nu) {
  d <- 2
  n <- nrow(x)
  const <- -lgamma(nu/2 + d/2) + lgamma(nu/2) + d/2*log(nu*pi) + 0.5*log(det(Sigma))
  term1 <- -(nu + d)/2 * sum(log(1 + (1/nu)*Mahalanobis(x, mu, Sigma, inverted = TRUE)))
  return(const + term1)
}
```

```{r, error=FALSE}
# function to estimate parameters using maximum likelihood
fit_t_params <- function(x) {
  start_vals <- list(mu = colMeans(x), Sigma = cov(x), nu = ncol(x))
  mle_obj <- mle(neg_loglik_t, start=start_vals, x = x, method = "L-BFGS-B",
                 lower = list(mu = rep(-Inf, 2), Sigma = diag(2)*1e-6, nu = 2+2),
                 upper = list(mu = rep(Inf, 2), Sigma = Inf, nu = Inf))
  return(mle_obj)
}
# Set TRUE/FALSE
debug(fit_t_params)
```

```{r}
# function to estimate ML parameters, SEs, and conduct hypothesis test
estimate_params <- function(x, true_params, alpha = 0.05) {
  # estimate ML parameters
  mle_obj <- fit_t_params(x)
  est_params <- coef(mle_obj)
  # compute standard errors
  vcov <- vcov(mle_obj)
  se_params <- sqrt(diag(vcov))
  
  # hypothesis test for Sigma
  df <- ncol(x)
  test_stat <- df * log(det(true_params$Sigma)) - df * log(det(est_params$Sigma)) -
    (df * (df+1))/2 * sum(log(true_params$Sigma)) + (df * (df+1))/2 * sum(log(est_params$Sigma))
  p_val <- pchisq(test_stat, df*(df+1)/2, lower.tail = FALSE)
  reject <- p_val < alpha
  
  return(list(est_params = est_params, se_params = se_params, 
              reject = reject, p_val = p_val))
}
```

b) To begin, we can import the “streaming_popul_data.csv” file into R and load the necessary libraries for our analysis.
To estimate the parameters of the Multivariate T distribution, we first need to load the data and install and load the mvtnorm package in R:
```{r}
# create separate datasets for each genre
rap_data <- data$rap
pop_data <- data$pop
metal_data <- data$metal
rock_data <- data$rock
```

```{r}
Next, we can create a matrix of the data for each genre:
library(mvtnorm)
# Create matrices of the data for each genre
rap_data <- as.matrix(data$rap)
pop_data <- as.matrix(data$pop)
metal_data <- as.matrix(data$metal)
rock_data <- as.matrix(data$rock)
# Calculate the mean vector and covariance matrix for the data
mu <- c(mean(rap_data), mean(pop_data), mean(metal_data), mean(rock_data))
mu
sigma <- cov(cbind(rap_data, pop_data, metal_data, rock_data))
sigma
```

Finally, we can use the rmvt() function from the mvtnorm package to generate samples from the Multivariate T distribution, with the estimated mean vector and covariance matrix. We can also calculate the sample means and variances for each genre from the samples, and use these to make inferences about the streaming success of each genre.
```{r}
# Generate 1000 samples from the Multivariate T distribution
set.seed(123) # Set a seed for reproducibility
samples <- rmvt(n = 1000, sigma = sigma, df = 5, delta = mu)
# Calculate the sample means and variances for each genre
sample_means <- apply(samples, 2, mean)
sample_vars <- apply(samples, 2, var)
# Print the sample means and standard errors for each genre
cat("Rap: mean =", round(sample_means[1], 2), "SE =", round(sqrt(sample_vars[1]/1000), 2), "\n")
cat("Pop: mean =", round(sample_means[2], 2), "SE =", round(sqrt(sample_vars[2]/1000), 2), "\n")
cat("Metal: mean =", round(sample_means[3], 2), "SE =", round(sqrt(sample_vars[3]/1000), 2), "\n")
cat("Rock: mean =", round(sample_means[4], 2), "SE =", round(sqrt(sample_vars[4]/1000), 2), "\n")
```

c)To test the hypothesis that the expected popularity is the same in every genre, we can use a one-way ANOVA (Analysis of Variance) test with a significance level of 0.05.
Here’s the R code to perform the test:
```{r}
# Create a data frame with the relevant columns
genres_data <- data.frame(
  rap = data$rap,
  pop = data$pop,
  metal = data$metal,
  rock = data$rock
)
# Perform a ANOVA test 
anova_result <- anova(lm(data$...1 ~ data$rap+data$metal+data$rock+data$pop, data = data))
print(anova_result)
```

```{r}
To test the hypothesis that the expected popularity is the same in every genre, we can use an analysis of variance (ANOVA) test in R. Here is the code:
# Extract the columns for each genre
rap_data <- data$rap
pop_data <- data$pop
metal_data <- data$metal
rock_data <- data$rock
# Perform ANOVA test
result <- aov(c(rap_data, pop_data, metal_data, rock_data) ~ rep(c("rap", "pop", "metal", "rock"), c(length(rap_data), length(pop_data), length(metal_data), length(rock_data))))
summary(result)
```


d)Population mean and hypothesis testing
```{r}
# Extract the pop data
pop_data <- data$pop
# Calculate the sample mean and standard deviation
pop_mean <- mean(pop_data)
pop_sd <- sd(pop_data)
# Set the null hypothesis mean
null_mean <- 40
# Calculate the t-value and the p-value
t_value <- (pop_mean - null_mean) / (pop_sd / sqrt(length(pop_data)))
p_value <- pt(t_value, df = length(pop_data) - 1, lower.tail = FALSE) * 2
```

```{r}
# Print the results
cat("t-value:", t_value, "\n")
cat("p-value:", p_value, "\n")
# Check if the null hypothesis is rejected or not
if (p_value < 0.05) {
  cat("Reject the null hypothesis. The expected popularity for pop songs is not equal to 40.\n")
} else {
  cat("Fail to reject the null hypothesis. The expected popularity for pop songs is equal to 40.\n")
}
```
