phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")
data_simple <- read.table(
  "https://paul-buerkner.github.io/data/data_simple.txt", 
  header = TRUE
)
head(data_simple)
A <- ape::vcv.phylo(phylo)

data_repeat <- read.table(
  "https://paul-buerkner.github.io/data/data_repeat.txt", 
  header = TRUE
)
data_repeat$spec_mean_cf <- 
  with(data_repeat, sapply(split(cofactor, phylo), mean)[phylo])
head(data_repeat)

model_repeat1 <- brm(
  phen ~ spec_mean_cf + (1|gr(phylo, cov = A)) + (1|species), 
  data = data_repeat, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2, 
  iter = 4000, warmup = 1000
)

------------------------------------------------
set.seed(5)
n = 10
n_condition = 5
ABC =
  tibble(
    condition = rep(c("A","B","C","D","E"), n),
    response = rnorm(n * 5, c(0,1,2,1,-1), 0.5)
  )
m = brm(
  response ~ (1|condition), 
  data = ABC, 
  prior = c(
    prior(normal(0, 1), class = Intercept),
    prior(student_t(3, 0, 1), class = sd),
    prior(student_t(3, 0, 1), class = sigma)
  ),
  control = list(adapt_delta = .99),
  file = "G:\\My Drive\\Philippine rodents\\Chrotomyini_brms\\sandboxModels\\tidy-brms_m.rds" # cache model (can be removed)  
)
get_variables(m)

m %>%
  spread_draws(r_condition[condition,term]) %>%
  head(10)

