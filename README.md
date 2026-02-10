## Example

Below is a minimal example using the included toy dataset.

```r
library(dplyr)
library(BARTSIMP)

############################################
# Load toy dataset included with package
############################################
data("toy_data")

############################################
# Prepare model inputs
############################################

# Predictor matrix
x_train <- toy_data %>%
  select(x1, x2, x3, x4, x5) %>%
  data.frame()

# Response
y_train <- toy_data$y

# Spatial coordinates
s_train <- toy_data %>% select(s1, s2)

# Replicate counts per unique spatial location
# (robust to row ordering)
size <- toy_data %>%
  distinct(s1, s2, size) %>%
  pull(size)

############################################
# Fit model
############################################
fit <- bartsimp(
  x.train = x_train,
  y.train = y_train,
  x.test = x_train,
  s1 = s_train$s1,
  s2 = s_train$s2,
  size = size,
  ntree = 10, # number of trees
  ndpost = 10, # number of posterior samples
  nwarmup = 10, # number of warmup samples 
  seed = 1
)


# spatial hyperparameters (spatial std, scale parameter and nugget std)
fit$sigmams
fit$kappas
fit$sigma

# fitted values
fit$yhat.train
```
