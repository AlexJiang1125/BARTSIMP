#' @export
library(INLA)
tunehyperpriors = function(x_train,
                           y_train,
                           s1,
                           s2) {
  s_train <- data.frame(s1 = s1, s2 = s2)
  s_train_unique <-
    s_train[!duplicated(s_train),] %>% as.data.frame()
  colnames(s_train_unique) <- c("s1_train", "s2_train")
  dat_train <- as.data.frame(cbind(y_train, x_train))
  # crude linear model
  lm_mod <- lm(y_train ~ ., data = dat_train)
  # crude overestimation of spatial and nonspatial parameters
  sigest <- sigma(lm_mod)
  sigmam_0 <- sigest * 2.5

  sigquant <- 0.9
  nu <- 3
  qchi <- qchisq(1.0 - sigquant, nu)
  lambda <- (sigest * qchi) / nu


  buff <- 0
  # get extreme spatial points in the datasets
  xmin <- min(s_train_unique$s1_train) - buff
  xmax <- max(s_train_unique$s1_train) + buff
  ymin <- min(s_train_unique$s2_train) - buff
  ymax <- max(s_train_unique$s2_train) + buff
  avg_lth <- (ymax - ymin + xmax - xmin) * 0.5

  rho_0 <- sqrt((ymax - ymin)^2 + (xmax - xmin)^2)/5
  alpha_1 <- alpha_2 <- 0.05

  ### spatial domain (xmin, xmax) x (ymin, ymax)
  mesh1 <- inla.mesh.2d(
    loc.domain = cbind(c(xmin, xmax, xmax, xmin), c(ymin, ymin, ymax, ymax)),
    max.edge = c(.1, .2) * avg_lth,
    offset = c(.1, .2) * avg_lth
  )

  loc_unique <- s_train_unique
  loc <- s_train

  mesh2 <- inla.mesh.2d(
    loc.domain = data.frame(loc_unique) ,
    max.edge = c(.1, .2) * avg_lth,
    cutoff = 0.05 * avg_lth
  )
  myspde <- inla.spde2.pcmatern(
    mesh = mesh2,
    prior.range = c(rho_0, alpha_2),
    prior.sigma = c(sigmam_0, alpha_1)
  )
  fyA <- inla.spde.make.A(mesh = mesh2, loc = as.matrix(loc))
  mystack <-
    inla.stack(
      data = list(resp = y_train),
      A = list(fyA, 1),
      effects = list(i = 1:myspde$n.spde,
                     df = x_train)
    )

  fml <- paste(colnames(x_train) , collapse = "+")
  fml <- paste0("resp ~ 0 + ", fml, "+ f(i, model=myspde)")
  fml <- as.formula(fml)
  result <- inla(
    fml,
    data = inla.stack.data(mystack),
    control.family = list(hyper = list(prec = list(
      prior = "loggamma",
      param = c(nu / 2, lambda *
                  nu / 2)
    ))),
    control.predictor = list(A = inla.stack.A(mystack))
  )
  mod_res <- result$summary.hyperpar
  hyperpars <- mod_res[, c("0.5quant")]
  sigest_new <- sqrt(1 / hyperpars[1])
  range_new <- hyperpars[2]
  sig_m_new <- hyperpars[3]
  print(paste0(
    "residual std: ",
    round(sigest_new, 2),
    " range: ",
    round(range_new, 2),
    " matern std: ",
    round(sig_m_new, 2)
  ))
  iters <- 5
  qchi = qchisq(1.0 - sigquant, nu)
  lambda_new = (sigest_new * sigest_new * qchi) / nu
  myspde <- inla.spde2.pcmatern(
    mesh = mesh2,
    prior.range = c(range_new / 2.5, alpha_2),
    prior.sigma = c(sig_m_new * 2.5, alpha_1)
  )

  for (iter in 1:iters) {
    result <- inla(
      fml,
      data = inla.stack.data(mystack),
      control.family = list(hyper = list(prec = list(
        prior = "loggamma",
        param = c(nu / 2, lambda_new *
                    nu / 2)
      ))),
      control.predictor = list(A = inla.stack.A(mystack))
    )
    mod_res <- result$summary.hyperpar
    hyperpars <- mod_res[, c("0.5quant")]
    sigest_new <- sqrt(1 / hyperpars[1])
    range_new <- hyperpars[2]
    sig_m_new <- hyperpars[3]
    qchi = qchisq(1.0 - sigquant, nu)
    lambda_new = (sigest_new * sigest_new * qchi) / nu
    print(paste0(
      "residual std: ",
      round(sigest_new, 2),
      " range: ",
      round(range_new, 2),
      " matern std: ",
      round(sig_m_new, 2)
    ))
    myspde <- inla.spde2.pcmatern(
      mesh = mesh2,
      prior.range = c(range_new / 2.5, alpha_2),
      prior.sigma = c(sig_m_new * 2.5, alpha_1)
    )

  }
  res <- list(
    sigest_new = sigest_new,
    range_new = range_new,
    sig_m_new = sig_m_new,
    sigmam_0 = sigmam_0,
    rho_0 = rho_0
  )
  return(res)
}
