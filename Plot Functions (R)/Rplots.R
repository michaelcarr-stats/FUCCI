library(readxl)
library(tidyverse)
library(ggpubr)

#Experiment Plots - CellTracking
setwd("Z:/R/ExperimentRuns/CellTrajectory")
theta_regressed <- read_excel("SMCABC_Samples_20.xlsx", sheet = "theta_regressed")
summaries <- read_excel("SMCABC_Samples_20.xlsx", sheet = "Samples") %>% select(contains("sx"))

Nred <- 566
Nyellow <- 111
Ngreen <- 166
RedDistance <- 105.38626
YellowDistance <- 39.95925
GreenDistance <- 100.36051

theta_regressed %>% 
  as_tibble() %>% 
  gather(Parameter, Value) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Rr", "Ry", "Rg","Mr", "My", "Mg"))) %>% 
  group_by(Parameter) %>% 
  summarise(
    mu = signif(mean(Value),3),
    sd = signif(sd(Value),3),
    "2.5%" = signif(quantile(Value, 0.025),3),
    "50%" = signif(quantile(Value, 0.5),3),
    "97.5%" = signif(quantile(Value, 0.975),3),
    CV = signif(sd/mu,3)
            )
  
p1 <- 
  theta_regressed %>% 
  as_tibble() %>%
  gather(Parameter,Value) %>% 
  mutate(Type = case_when(
    Parameter == "Rr" | Parameter == "Ry" | Parameter == "Rg" ~ "Transition",
    Parameter == "Mr" | Parameter == "My" | Parameter == "Mg" ~ "Motility" 
  )) %>% 
  mutate(Type = factor(Type, levels = c("Transition", "Motility"))) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Rr", "Ry", "Rg","Mr", "My", "Mg"))) %>% 
  filter(Type == "Transition") %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  facet_wrap(~Type, scale = "free") +
  scale_fill_manual(values = c("Rr" = "red", "Ry" = "gold", "Rg" = "green", "Mr" = "red", "My" = "gold", "Mg" = "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g])))+ 
  scale_linetype(labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(x = expression("h"^"-1"), title = "(a)")

p2 <- 
  theta_regressed %>% 
  as_tibble() %>%
  gather(Parameter,Value) %>% 
  mutate(Type = case_when(
    Parameter == "Rr" | Parameter == "Ry" | Parameter == "Rg" ~ "Transition",
    Parameter == "Mr" | Parameter == "My" | Parameter == "Mg" ~ "Motility" 
  )) %>% 
  mutate(Type = factor(Type, levels = c("Transition", "Motility"))) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Rr", "Ry", "Rg","Mr", "My", "Mg"))) %>% 
  filter(Type == "Motility") %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  facet_wrap(~Type, scale = "free") +
  scale_fill_manual(values = c("Rr" = "red", "Ry" = "gold", "Rg" = "green", "Mr" = "red", "My" = "gold", "Mg" = "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g])))+ 
  scale_linetype(labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(x = expression("h"^"-1"), title = "(b)")

p3 <- 
  summaries %>% 
  gather(Parameter, Value) %>% 
  mutate(Type = case_when(
    Parameter == "sx1" | Parameter == "sx2" | Parameter == "sx3" ~ "Transition",
    Parameter == "sx4" | Parameter == "sx5" | Parameter == "sx6" ~ "Motility" 
  )) %>% 
  mutate(Type = factor(Type, levels = c("Transition", "Motility"))) %>% 
  mutate(TrueVal = case_when(
    Parameter == "sx1" ~ Nred,
    Parameter == "sx2" ~ Nyellow,
    Parameter == "sx3" ~ Ngreen,
    Parameter == "sx4" ~ RedDistance,
    Parameter == "sx5" ~ YellowDistance,
    Parameter == "sx6" ~ GreenDistance
  )) %>% 
  filter(Type == "Transition") %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  facet_wrap(~Type, scale = "free") + 
  geom_vline(aes(xintercept = TrueVal, color = Parameter), linetype = "dashed") +
  scale_fill_manual(values = c("sx1" = "red", "sx2" = "gold", "sx3" = "green", "sx4" = "red", "sx5" = "gold", "sx6" = "green"), labels = c(expression(S[1]),expression(S[2]),expression(S[3])))+ 
  scale_color_manual(values = c("sx1" = "red", "sx2" = "goldenrod3", "sx3" = "green", "sx4" = "red", "sx5" = "goldenrod3", "sx6" = "green"), labels = c(expression(S[1]),expression(S[2]),expression(S[3])))+ 
  scale_linetype(labels = c(expression(S[1]),expression(S[2]),expression(S[3]))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(x = "cells", title = "(c)")

p4 <- 
  summaries %>% 
  gather(Parameter, Value) %>% 
  mutate(Type = case_when(
    Parameter == "sx1" | Parameter == "sx2" | Parameter == "sx3" ~ "Transition",
    Parameter == "sx4" | Parameter == "sx5" | Parameter == "sx6" ~ "Motility" 
  )) %>% 
  mutate(Type = factor(Type, levels = c("Transition", "Motility"))) %>% 
  mutate(TrueVal = case_when(
    Parameter == "sx1" ~ Nred,
    Parameter == "sx2" ~ Nyellow,
    Parameter == "sx3" ~ Ngreen,
    Parameter == "sx4" ~ RedDistance,
    Parameter == "sx5" ~ YellowDistance,
    Parameter == "sx6" ~ GreenDistance
  )) %>% 
  filter(Type == "Motility") %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  facet_wrap(~Type, scale = "free") + 
  geom_vline(aes(xintercept = TrueVal, color = Parameter), linetype = "dashed") +
  scale_fill_manual(values = c("sx1" = "red", "sx2" = "gold", "sx3" = "green", "sx4" = "red", "sx5" = "gold", "sx6" = "green"), labels = c(expression(S[4]),expression(S[5]),expression(S[6])))+ 
  scale_color_manual(values = c("sx1" = "red", "sx2" = "goldenrod3", "sx3" = "green", "sx4" = "red", "sx5" = "goldenrod3", "sx6" = "green"), labels = c(expression(S[4]),expression(S[5]),expression(S[6])))+ 
  scale_linetype(labels = c(expression(S[4]),expression(S[5]),expression(S[6]))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(x = "cells", title = "(c)")

ggarrange(p1,p2,p3,p4, nrow = 1)
ggsave("ExperimentResults_CellTracking.pdf", path = "Z:/R", width = 9, height = 4, units = "in")

#Experiment Plots - SpatioTemporal
setwd("Z:/R/ExperimentRuns/SpatioTemporal")
theta_regressed <- read_excel("SMCABC_Samples_20.xlsx", sheet = "theta_regressed")
summaries <- read_excel("SMCABC_Samples_20.xlsx", sheet = "Samples") %>% select(contains("sx"))

p1 <- 
  theta_regressed %>% 
  as_tibble() %>%
  gather(Parameter,Value) %>% 
  mutate(Type = case_when(
    Parameter == "Rr" | Parameter == "Ry" | Parameter == "Rg" ~ "Transition",
    Parameter == "Mr" | Parameter == "My" | Parameter == "Mg" ~ "Motility" 
  )) %>% 
  mutate(Type = factor(Type, levels = c("Transition", "Motility"))) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Rr", "Ry", "Rg","Mr", "My", "Mg"))) %>% 
  filter(Type == "Transition") %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  facet_wrap(~Type, scale = "free") +
  scale_fill_manual(values = c("Rr" = "red", "Ry" = "gold", "Rg" = "green", "Mr" = "red", "My" = "gold", "Mg" = "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g])))+ 
  scale_linetype(labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(x = expression("h"^"-1"), title = "(a)")

p2 <- 
  theta_regressed %>% 
  as_tibble() %>%
  gather(Parameter,Value) %>% 
  mutate(Type = case_when(
    Parameter == "Rr" | Parameter == "Ry" | Parameter == "Rg" ~ "Transition",
    Parameter == "Mr" | Parameter == "My" | Parameter == "Mg" ~ "Motility" 
  )) %>% 
  mutate(Type = factor(Type, levels = c("Transition", "Motility"))) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Rr", "Ry", "Rg","Mr", "My", "Mg"))) %>% 
  filter(Type == "Motility") %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  facet_wrap(~Type, scale = "free") +
  scale_fill_manual(values = c("Rr" = "red", "Ry" = "gold", "Rg" = "green", "Mr" = "red", "My" = "gold", "Mg" = "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g])))+ 
  scale_linetype(labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(x = expression("h"^"-1"), title = "(b)")

ggarrange(p1,p2, nrow = 1)
ggsave("ExperimentResults_SpatioTemporal.pdf", path = "Z:/R", width = 4.5, height = 4, units = "in")

#SixParameters Plots
setwd("Z:/R/SyntheticRuns/SixParameters")
Sim1 <- read_excel("Sim1/SMCABC_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 1)
Sim2 <- read_excel("Sim2/SMCABC_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 2)
Sim3 <- read_excel("Sim3/SMCABC_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 3)
Sim4 <- read_excel("Sim4/SMCABC_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 4)

p1 <-
  bind_rows(Sim1,Sim2,Sim3,Sim4) %>% 
  select(c(1:3,ncol(Sim1))) %>% 
  gather(Parameter, Value, -Simulation) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Rr", "Ry", "Rg"))) %>% 
  mutate(Simulation = case_when(Simulation == 1 ~ "(a)", Simulation == 2 ~ "(b)", Simulation == 3 ~ "(c)", Simulation == 4 ~ "(d)")) %>% 
  ggplot(aes(x = Value, fill = Parameter)) + 
  geom_density(alpha = .3) + 
  facet_wrap(~Simulation, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g])))+ 
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  )
  
p2 <- 
  bind_rows(Sim1,Sim2,Sim3,Sim4) %>% 
  select(c(4:6,ncol(Sim1))) %>% 
  gather(Parameter, Value, -Simulation) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Mr", "My", "Mg"))) %>%
  mutate(Simulation = case_when(Simulation == 1 ~ "(e)", Simulation == 2 ~ "(f)", Simulation == 3 ~ "(g)", Simulation == 4 ~ "(h)")) %>% 
  ggplot(aes(x = Value, fill = Parameter)) + 
  geom_density(alpha = .3) + 
  facet_wrap(~Simulation, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +  
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  )
ggarrange(p1,p2, nrow = 2)

#Proliferation Plots
setwd("Z:/R/SyntheticRuns/Proliferation")
Sim1 <- read_excel("Sim1/SMCABC_3params_Samples.xlsx", sheet = "Samples_regressed") %>% mutate(Simulation = 1)
Sim2 <- read_excel("Sim2/SMCABC_3params_Samples.xlsx", sheet = "Samples_regressed") %>% mutate(Simulation = 2)
Sim3 <- read_excel("Sim3/SMCABC_3params_Samples.xlsx", sheet = "Samples_regressed") %>% mutate(Simulation = 3)
Sim4 <- read_excel("Sim4/SMCABC_3params_Samples.xlsx", sheet = "Samples_regressed") %>% mutate(Simulation = 4)

SixParameterData <-
  bind_rows(Sim1,Sim2,Sim3,Sim4) %>% 
  gather(Parameter, Value, -Simulation) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Rr", "Ry", "Rg","Mr", "My", "Mg"))) %>% 
  mutate(TrueValue = case_when(
    Simulation == 1 & Parameter == "Mr" ~ 4,
    Simulation == 1 & Parameter == "My" ~ 4,
    Simulation == 1 & Parameter == "Mg" ~ 4,
    Simulation == 2 & Parameter == "Mr" ~ 2,
    Simulation == 2 & Parameter == "My" ~ 5,
    Simulation == 2 & Parameter == "Mg" ~ 8,
    Simulation == 3 & Parameter == "Mr" ~ 8,
    Simulation == 3 & Parameter == "My" ~ 2,
    Simulation == 3 & Parameter == "Mg" ~ 5,
    Simulation == 4 & Parameter == "Mr" ~ 5,
    Simulation == 4 & Parameter == "My" ~ 8,
    Simulation == 4 & Parameter == "Mg" ~ 2,
    Parameter == "Rr" ~ 0.04,
    Parameter == "Ry" ~ 0.17,
    Parameter == "Rg" ~ 0.08
  )) %>% 
  mutate(Parameter_type = case_when(
    Parameter %in% c("Rr", "Ry", "Rg") ~ "Transition",
    Parameter %in% c("Mr", "My", "Mg") ~ "Motility",
  )) %>% 
  mutate(Parameter_type = factor(Parameter_type, levels = c("Transition", "Motility"))) %>% 
  mutate(facets = case_when(
    Simulation == 1 & Parameter_type == "Transition" ~ "(a)",
    Simulation == 2 & Parameter_type == "Transition" ~ "(b)",
    Simulation == 3 & Parameter_type == "Transition" ~ "(c)",
    Simulation == 4 & Parameter_type == "Transition" ~ "(d)",
    Simulation == 1 & Parameter_type == "Motility" ~ "(e)",
    Simulation == 2 & Parameter_type == "Motility" ~ "(f)",
    Simulation == 3 & Parameter_type == "Motility" ~ "(g)",
    Simulation == 4 & Parameter_type == "Motility" ~ "(h)"
  )) 

p1 <- 
  SixParameterData %>% 
  filter(Parameter %in% c("Rr", "Ry", "Rg")) %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  geom_vline(aes(xintercept = TrueValue, color = Parameter), linetype='dashed') + 
  facet_wrap(~facets, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  )

p2 <- 
  SixParameterData %>% 
  filter(Parameter %in% c("Mr", "My", "Mg")) %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) + 
  geom_vline(aes(xintercept = TrueValue, color = Parameter), linetype='dashed') + 
  facet_wrap(~facets, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  )

ggarrange(p1,p2, nrow = 2)
ggsave("SixParameters_SIM.pdf", path = "Z:/R/SyntheticRuns", width = 9, height = 8, units = "in")


#Motility Plots
setwd("Z:/R/SyntheticRuns/Motility")
Sim1_10 <- read_excel("Sim1/SMCABC_3params_Samples_10.xlsx", sheet = "Samples") %>% mutate(Simulation = 1, ntrack = 10)
Sim1_20 <- read_excel("Sim1/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 1, ntrack = 20)
Sim1_30 <- read_excel("Sim1/SMCABC_3params_Samples_30.xlsx", sheet = "Samples") %>% mutate(Simulation = 1, ntrack = 30)
Sim1_40 <- read_excel("Sim1/SMCABC_3params_Samples_40.xlsx", sheet = "Samples") %>% mutate(Simulation = 1, ntrack = 40)
Sim1_50 <- read_excel("Sim1/SMCABC_3params_Samples_50.xlsx", sheet = "Samples") %>% mutate(Simulation = 1, ntrack = 50)
Sim1 <- bind_rows(Sim1_10, Sim1_20, Sim1_30, Sim1_40, Sim1_50)

Sim2_10 <- read_excel("Sim2/SMCABC_3params_Samples_10.xlsx", sheet = "Samples") %>% mutate(Simulation = 2, ntrack = 10)
Sim2_20 <- read_excel("Sim2/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 2, ntrack = 20)
Sim2_30 <- read_excel("Sim2/SMCABC_3params_Samples_30.xlsx", sheet = "Samples") %>% mutate(Simulation = 2, ntrack = 30)
Sim2_40 <- read_excel("Sim2/SMCABC_3params_Samples_40.xlsx", sheet = "Samples") %>% mutate(Simulation = 2, ntrack = 40)
Sim2_50 <- read_excel("Sim2/SMCABC_3params_Samples_50.xlsx", sheet = "Samples") %>% mutate(Simulation = 2, ntrack = 50)
Sim2 <- bind_rows(Sim2_10, Sim2_20, Sim2_30, Sim2_40, Sim2_50)

Sim3_10 <- read_excel("Sim3/SMCABC_3params_Samples_10.xlsx", sheet = "Samples") %>% mutate(Simulation = 3, ntrack = 10)
Sim3_20 <- read_excel("Sim3/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 3, ntrack = 20)
Sim3_30 <- read_excel("Sim3/SMCABC_3params_Samples_30.xlsx", sheet = "Samples") %>% mutate(Simulation = 3, ntrack = 30)
Sim3_40 <- read_excel("Sim3/SMCABC_3params_Samples_40.xlsx", sheet = "Samples") %>% mutate(Simulation = 3, ntrack = 40)
Sim3_50 <- read_excel("Sim3/SMCABC_3params_Samples_50.xlsx", sheet = "Samples") %>% mutate(Simulation = 3, ntrack = 50)
Sim3 <- bind_rows(Sim3_10, Sim3_20, Sim3_30, Sim3_40, Sim3_50)

Sim4_10 <- read_excel("Sim4/SMCABC_3params_Samples_10.xlsx", sheet = "Samples") %>% mutate(Simulation = 4, ntrack = 10)
Sim4_20 <- read_excel("Sim4/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 4, ntrack = 20)
Sim4_30 <- read_excel("Sim4/SMCABC_3params_Samples_30.xlsx", sheet = "Samples") %>% mutate(Simulation = 4, ntrack = 30)
Sim4_40 <- read_excel("Sim4/SMCABC_3params_Samples_40.xlsx", sheet = "Samples") %>% mutate(Simulation = 4, ntrack = 40)
Sim4_50 <- read_excel("Sim4/SMCABC_3params_Samples_50.xlsx", sheet = "Samples") %>% mutate(Simulation = 4, ntrack = 50)
Sim4 <- bind_rows(Sim4_10, Sim4_20, Sim4_30, Sim4_40, Sim4_50)

MotilityData <-
  bind_rows(Sim1,Sim2,Sim3,Sim4) %>% 
  gather(Parameter, Value, -c(Simulation, ntrack)) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Mr", "My", "Mg"))) %>% 
  mutate(TrueValue = case_when(
    Simulation == 1 & Parameter == "Mr" ~ 4,
    Simulation == 1 & Parameter == "My" ~ 4,
    Simulation == 1 & Parameter == "Mg" ~ 4,
    Simulation == 2 & Parameter == "Mr" ~ 2,
    Simulation == 2 & Parameter == "My" ~ 5,
    Simulation == 2 & Parameter == "Mg" ~ 8,
    Simulation == 3 & Parameter == "Mr" ~ 8,
    Simulation == 3 & Parameter == "My" ~ 2,
    Simulation == 3 & Parameter == "Mg" ~ 5,
    Simulation == 4 & Parameter == "Mr" ~ 5,
    Simulation == 4 & Parameter == "My" ~ 8,
    Simulation == 4 & Parameter == "Mg" ~ 2,
  )) %>% 
  mutate(Parameter_type = case_when(
    Parameter %in% c("Mr", "My", "Mg") ~ "Motility",
  )) %>% 
  mutate(Parameter_type = factor(Parameter_type, levels = c("Transition", "Motility"))) %>% 
  mutate(facets = case_when(
    Simulation == 1 & Parameter_type == "Motility" ~ "(a)",
    Simulation == 2 & Parameter_type == "Motility" ~ "(b)",
    Simulation == 3 & Parameter_type == "Motility" ~ "(c)",
    Simulation == 4 & Parameter_type == "Motility" ~ "(d)"
  )) 


MotilityData %>% 
  filter(Parameter %in% c("Mr", "My", "Mg")) %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) + 
  geom_vline(aes(xintercept = TrueValue, color = Parameter), linetype='dashed') + 
  facet_wrap(~facets, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  )
ggsave("Motility_Sim.pdf", path = "Z:/R/SyntheticRuns", width = 9, height = 4, units = "in")


#CellTracking vs CellDensity Plots
setwd("Z:/R/SyntheticRuns/Motility")
Sim1_CellTracking <- read_excel("Sim1/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 1, ntrack = 20, CellTracking = T)
Sim2_CellTracking <- read_excel("Sim2/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 2, ntrack = 20, CellTracking = T)
Sim3_CellTracking <- read_excel("Sim3/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 3, ntrack = 20, CellTracking = T)
Sim4_CellTracking <- read_excel("Sim4/SMCABC_3params_Samples_20.xlsx", sheet = "Samples") %>% mutate(Simulation = 4, ntrack = 20, CellTracking = T)

setwd("Z:/R/SyntheticRuns/Motility/SpatioTemporal")
Sim1_CellDensity <- read_excel("Sim1/SMCABC_3params_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 1, ntrack = NA, CellTracking = F)
Sim2_CellDensity <- read_excel("Sim2/SMCABC_3params_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 2, ntrack = NA, CellTracking = F)
Sim3_CellDensity <- read_excel("Sim3/SMCABC_3params_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 3, ntrack = NA, CellTracking = F)
Sim4_CellDensity <- read_excel("Sim4/SMCABC_3params_Samples.xlsx", sheet = "Samples") %>% mutate(Simulation = 4, ntrack = NA, CellTracking = F)

MotilityData <-
  bind_rows(Sim1_CellTracking,Sim2_CellTracking,Sim3_CellTracking,Sim4_CellTracking,Sim1_CellDensity,Sim2_CellDensity,Sim3_CellDensity,Sim4_CellDensity) %>% 
  gather(Parameter, Value, -c(Simulation, ntrack, CellTracking)) %>% 
  mutate(Parameter = factor(Parameter, levels = c("Mr", "My", "Mg"))) %>% 
  mutate(TrueValue = case_when(
    Simulation == 1 & Parameter == "Mr" ~ 4,
    Simulation == 1 & Parameter == "My" ~ 4,
    Simulation == 1 & Parameter == "Mg" ~ 4,
    Simulation == 2 & Parameter == "Mr" ~ 2,
    Simulation == 2 & Parameter == "My" ~ 5,
    Simulation == 2 & Parameter == "Mg" ~ 8,
    Simulation == 3 & Parameter == "Mr" ~ 8,
    Simulation == 3 & Parameter == "My" ~ 2,
    Simulation == 3 & Parameter == "Mg" ~ 5,
    Simulation == 4 & Parameter == "Mr" ~ 5,
    Simulation == 4 & Parameter == "My" ~ 8,
    Simulation == 4 & Parameter == "Mg" ~ 2,
  )) %>% 
  mutate(Parameter_type = case_when(
    Parameter %in% c("Mr", "My", "Mg") ~ "Motility",
  )) %>% 
  mutate(Parameter_type = factor(Parameter_type, levels = c("Transition", "Motility"))) %>% 
  mutate(facets = case_when(
    Simulation == 1 & Parameter_type == "Motility" ~ "(a)",
    Simulation == 2 & Parameter_type == "Motility" ~ "(b)",
    Simulation == 3 & Parameter_type == "Motility" ~ "(c)",
    Simulation == 4 & Parameter_type == "Motility" ~ "(d)"
  )) 

MotilityData %>% 
  filter(Parameter %in% c("Mr", "My", "Mg")) %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) + 
  geom_vline(aes(xintercept = TrueValue, color = Parameter), linetype='dashed') + 
  facet_wrap(~facets, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c(expression(M[r]),expression(M[y]),expression(M[g]))) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  )
ggsave("Motility_Sim.pdf", path = "Z:/R/SyntheticRuns", width = 9, height = 4, units = "in")
