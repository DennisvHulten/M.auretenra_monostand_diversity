# Data handling
library(tidyverse)

# Mixed models
library(lme4)
library(lmerTest)

# Model diagnostics
library(performance)

# Post-hoc tests
library(emmeans)
library(tidyr)

#### Full Surveyed Area #######
df <- read.table(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/Monostand_metrics_combined.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Convert to factors
df$Site <- factor(df$Site)
df$Environment <- factor(df$Environment)
rugosity_lookup <- tibble(
  Site = c("Playa Kalki","Playa Kalki","Playa Kalki",
           "Snakebay","Snakebay","Snakebay",
           "Sea Aquarium","Sea Aquarium"),
  Environment = c("Slope","Crest","Terrace",
                  "Slope","Crest","Terrace",
                  "Slope","Crest"),
  rugosity = c(1.785,1.784,1.401,
               1.703,1.479,1.282,
               1.791,1.718)
)

df <- df %>%
  left_join(rugosity_lookup, by = c("Site", "Environment"))
summary_df <- df %>%
  filter(area_2d >= 3.14) %>%
  group_by(Site, Environment) %>%
  summarise(
    area_2d_max   = max(area_2d, na.rm = TRUE),
    area_2d_avg   = mean(area_2d, na.rm = TRUE),
    perimeter_avg = mean(perimeter, na.rm = TRUE),
    perimeter_max = max(perimeter, na.rm = TRUE),
    .groups = "drop"
  )


df_merged <- df %>%
  group_by(across(-c(Environment, rugosity))) %>%
  summarise(
    Environment = paste(sort(unique(Environment)), collapse = ". "),
    rugosity = mean(rugosity, na.rm = TRUE),
    .groups = "drop"
  )

df2 <- df_merged %>%
  mutate(size_class = case_when(
    area_2d < 1               ~ "individual",
    area_2d >= 1 & area_2d <= 3.14 ~ "patch",
    area_2d > 3.14               ~ "stand"
  ))
counts <- df2 %>%
  group_by(Site, size_class, Environment) %>%
  summarise(n = n(), .groups = "drop")
counts_wide <- counts %>%
  pivot_wider(
    names_from = size_class,
    values_from = n,
    values_fill = 0
  )

counts_wide
area_by_class <- df2 %>%
  group_by(Site, size_class) %>%
  summarise(
    total_area_m2 = sum(area_2d, na.rm = TRUE),
    n_shapes = n(),
    .groups = "drop"
  )
area_by_class

area_by_environment <- df2 %>%
  group_by(Site, Environment) %>%
  summarise(
    total_area_m2 = sum(area_2d, na.rm = TRUE),
    n_shapes = n(),
    .groups = "drop"
  )
area_by_environment

area_by_class_by_environment <- df2 %>%
  group_by(Site, Environment, size_class) %>%
  summarise(
    total_area_m2 = sum(area_2d, na.rm = TRUE),
    n_shapes = n(),
    .groups = "drop"
  )

area_by_class_by_environment

df <- df2 %>%
  filter(area_2d >= 3.14)

# Quick sanity check
str(df)
summary(df)
# Each Site × Environment should have exactly one rugosity value
df %>%
  group_by(Site, Environment) %>%
  summarise(n = n(), rugosity = unique(rugosity))

# Area distributions
ggplot(df2, aes(Environment, area_2d, fill = Environment)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Monostand area (m²)")

# Perimeter distributions
ggplot(df2, aes(Environment, perimeter, fill = Environment)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Monostand perimeter (m)")

m_area <- lmer(
  log(area_2d) ~ Environment + rugosity + (1 | Site),
  data = df2
)

summary(m_area)
anova(m_area)
check_model(m_area)
performance::check_collinearity(m_area)
emmeans(m_area, pairwise ~ Environment)

model_gamma <- glm(area_2d ~ Site,
                   family = Gamma(link = "log"),
                   data = df2)

summary(model_gamma)
anova(model_gamma, test = "Chisq")

model_env_rug <- glm(area_2d ~ Environment + rugosity,
                     family = Gamma(link = "log"),
                     data = df2)

summary(model_env_rug)
anova(model_env_rug, test = "Chisq")

###high correlation between env and rugo
# Only Environment
m_area_env <- lmer(log(area_2d) ~ Environment + (1 | Site), data = df)
summary(m_area_env)
anova(m_area_env)
check_model(m_area_env)
performance::check_collinearity(m_area_env)

# Only rugosity
m_area_rug <- lmer(log(area_2d) ~ rugosity + (1 | Site), data = df)
summary(m_area_rug)
anova(m_area_rug)
check_model(m_area_rug)
performance::check_collinearity(m_area_rug)

m_area <- lmer(
  log(area_2d) ~ Environment + rugosity + (1 | Site/Environment),
  data = df
)

summary(m_area)
anova(m_area)
check_model(m_area)
performance::check_collinearity(m_area)

m_perim <- lmer(
  log(perimeter) ~ Environment + rugosity + (1 | Site),
  data = df
)

summary(m_perim)
anova(m_perim)
check_model(m_perim)

emm_area <- emmeans(m_area, pairwise ~ Environment)
emm_area$contrasts

emm_perim <- emmeans(m_perim, pairwise ~ Environment)
emm_perim$contrasts

df_counts <- df %>%
  group_by(Site, Environment) %>%
  summarise(n_monostands = n(), .groups = "drop")

m_count <- glmer(
  n_monostands ~ Environment + rugosity + (1 | Site),
  data = df_counts,
  family = poisson
)

summary(m_count)
check_overdispersion(m_count)

library(glmmTMB)

m_count_nb <- glmmTMB(
  n_monostands ~ Environment + rugosity + (1 | Site),
  data = df_counts,
  family = nbinom2
)

summary(m_count_nb)

# Area vs rugosity
ggplot(df, aes(rugosity, area_2d, color = Environment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10() +
  theme_minimal()

# Perimeter vs rugosity
ggplot(df, aes(rugosity, perimeter, color = Environment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10() +
  theme_minimal()


#### Cropped Area Only #########
df <- read.table(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/Monostand_metrics_combined_clipped.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Convert to factors
df$Site <- factor(df$Site)
df$Environment <- factor(df$Environment)
rugosity_lookup <- tibble(
  Site = c("Playa Kalki","Playa Kalki","Playa Kalki",
           "Snakebay","Snakebay","Snakebay",
           "Sea Aquarium","Sea Aquarium"),
  Environment = c("Slope","Crest","Terrace",
                  "Slope","Crest","Terrace",
                  "Slope","Crest"),
  rugosity = c(1.785,1.784,1.401,
               1.703,1.479,1.282,
               1.791,1.718)
)

df <- df %>%
  left_join(rugosity_lookup, by = c("Site", "Environment"))
summary_df <- df %>%
  filter(area_2d >= 3.14) %>%
  group_by(Site, Environment) %>%
  summarise(
    area_2d_max   = max(area_2d, na.rm = TRUE),
    area_2d_avg   = mean(area_2d, na.rm = TRUE),
    perimeter_avg = mean(perimeter, na.rm = TRUE),
    perimeter_max = max(perimeter, na.rm = TRUE),
    .groups = "drop"
  )


df_merged <- df %>%
  group_by(across(-c(Environment, rugosity))) %>%
  summarise(
    Environment = paste(sort(unique(Environment)), collapse = ". "),
    rugosity = mean(rugosity, na.rm = TRUE),
    .groups = "drop"
  )

df2 <- df_merged %>%
  mutate(size_class = case_when(
    area_2d < 1               ~ "individual",
    area_2d >= 1 & area_2d <= 3.14 ~ "patch",
    area_2d > 3.14               ~ "stand"
  ))
counts <- df2 %>%
  group_by(Site, size_class, Environment) %>%
  summarise(n = n(), .groups = "drop")
counts_wide <- counts %>%
  pivot_wider(
    names_from = size_class,
    values_from = n,
    values_fill = 0
  )

counts_wide
area_by_class_cropped <- df2 %>%
  group_by(Site, size_class) %>%
  summarise(
    total_area_m2 = sum(area_2d, na.rm = TRUE),
    n_shapes = n(),
    .groups = "drop"
  )
area_by_class_cropped

area_by_environment_cropped <- df2 %>%
  group_by(Site, Environment) %>%
  summarise(
    total_area_m2 = sum(area_2d, na.rm = TRUE),
    n_shapes = n(),
    .groups = "drop"
  )
area_by_environment_cropped


df <- df2 %>%
  filter(area_2d >= 3.14)

# Quick sanity check
str(df)
summary(df)
# Each Site × Environment should have exactly one rugosity value
df %>%
  group_by(Site, Environment) %>%
  summarise(n = n(), rugosity = unique(rugosity))

# Area distributions
ggplot(df2, aes(Environment, area_2d, fill = Environment)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Monostand area (m²)")

# Perimeter distributions
ggplot(df2, aes(Environment, perimeter, fill = Environment)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Monostand perimeter (m)")

m_area <- lmer(
  log(area_2d) ~ Environment + rugosity + (1 | Site),
  data = df2
)

summary(m_area)
anova(m_area)
check_model(m_area)
performance::check_collinearity(m_area)
emmeans(m_area, pairwise ~ Environment)

model_gamma <- glm(area_2d ~ Site,
                   family = Gamma(link = "log"),
                   data = df2)

summary(model_gamma)
anova(model_gamma, test = "Chisq")

model_env_rug <- glm(area_2d ~ Environment + rugosity,
                     family = Gamma(link = "log"),
                     data = df2)

summary(model_env_rug)
anova(model_env_rug, test = "Chisq")

###high correlation between env and rugo
# Only Environment
m_area_env <- lmer(log(area_2d) ~ Environment + (1 | Site), data = df)
summary(m_area_env)
anova(m_area_env)
check_model(m_area_env)
performance::check_collinearity(m_area_env)

# Only rugosity
m_area_rug <- lmer(log(area_2d) ~ rugosity + (1 | Site), data = df)
summary(m_area_rug)
anova(m_area_rug)
check_model(m_area_rug)
performance::check_collinearity(m_area_rug)

m_area <- lmer(
  log(area_2d) ~ Environment + rugosity + (1 | Site/Environment),
  data = df
)

summary(m_area)
anova(m_area)
check_model(m_area)
performance::check_collinearity(m_area)

m_perim <- lmer(
  log(perimeter) ~ Environment + rugosity + (1 | Site),
  data = df
)

summary(m_perim)
anova(m_perim)
check_model(m_perim)

emm_area <- emmeans(m_area, pairwise ~ Environment)
emm_area$contrasts

emm_perim <- emmeans(m_perim, pairwise ~ Environment)
emm_perim$contrasts

df_counts <- df %>%
  group_by(Site, Environment) %>%
  summarise(n_monostands = n(), .groups = "drop")

m_count <- glmer(
  n_monostands ~ Environment + rugosity + (1 | Site),
  data = df_counts,
  family = poisson
)

summary(m_count)
check_overdispersion(m_count)

library(glmmTMB)

m_count_nb <- glmmTMB(
  n_monostands ~ Environment + rugosity + (1 | Site),
  data = df_counts,
  family = nbinom2
)

summary(m_count_nb)

# Area vs rugosity
ggplot(df, aes(rugosity, area_2d, color = Environment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10() +
  theme_minimal()

# Perimeter vs rugosity
ggplot(df, aes(rugosity, perimeter, color = Environment)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  scale_y_log10() +
  theme_minimal()


