library(readxl)
library(tidyverse)
library(ggpubr)

#set working directory to the directory of this R file

#load data
Experiment_CellTracking <- read_excel("../Data/SMCABC returned data/SMCABC_DATA.xlsx", sheet = "Experiment_CellTracking")
Experiment_CellDensity <- read_excel("../Data/SMCABC returned data/SMCABC_DATA.xlsx", sheet = "Experiment_CellDensity")
ProliferationData <- read_excel("../Data/SMCABC returned data/SMCABC_DATA.xlsx", sheet = "ProliferationData")
MotilityData <- read_excel("../Data/SMCABC returned data/SMCABC_DATA.xlsx", sheet = "MotilityData")
SixParameterData <- read_excel("../Data/SMCABC returned data/SMCABC_DATA.xlsx", sheet = "SixParameterData")
CellTracking <- read_excel("../Data/DataProcessing/FUCCI_processed.xlsx", sheet = "CellTracking")


#observed summary statistics
Nred <- 566
Nyellow <- 111
Ngreen <- 166
RedDistance <- 105.38626
YellowDistance <- 39.95925
GreenDistance <- 100.36051


#Experiment Plots - Cell Tracking

Experiment_CellTracking %>% 
  as_tibble() %>% 
  select(1:6) %>% 
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
  Experiment_CellTracking %>% 
  as_tibble() %>%
  select(1:6) %>% 
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
  Experiment_CellTracking %>% 
  as_tibble() %>%
  select(1:6) %>% 
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
  Experiment_CellTracking %>% 
  select(contains("sx")) %>% 
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
  Experiment_CellTracking %>% 
  select(contains("sx")) %>% 
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
ggsave("ExperimentResults_CellTracking.pdf", width = 9, height = 4, units = "in")

#Experiment Plots - SpatioTemporal
p1 <- 
  Experiment_CellDensity %>% 
  as_tibble() %>%
  select(1:6) %>% 
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
  Experiment_CellDensity %>% 
  as_tibble() %>%
  select(1:6) %>% 
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
ggsave("ExperimentResults_SpatioTemporal.pdf", width = 4.5, height = 4, units = "in")

#SixParameters Plots - CellTracking

SixParameterData_CellTracking <-
  SixParameterData %>% 
  filter(MotilityData == "CellTracking") %>% 
  select(c(1:6,"Simulation")) %>% 
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
  SixParameterData_CellTracking %>%
  filter(Parameter_type == "Transition") %>% 
  ggplot(aes(x = Value, fill = Parameter)) + 
  geom_density(alpha = .3) + 
  facet_wrap(~facets, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g])))+ 
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  ) + 
  xlim(0,0.4)

p2 <- 
  SixParameterData_CellTracking %>%
  filter(Parameter_type == "Motility")  %>% 
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

ggsave("SixParameters_SIM_CellTracking.pdf", width = 9, height = 8, units = "in")

#SixParameters Plots - CellDensity

SixParameterData_CellDensity <-
  SixParameterData %>% 
  filter(MotilityData == "CellDensity") %>% 
  select(c(1:6,"Simulation")) %>% 
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
  SixParameterData_CellDensity %>%
  filter(Parameter_type == "Transition") %>% 
  ggplot(aes(x = Value, fill = Parameter)) + 
  geom_density(alpha = .3) + 
  facet_wrap(~facets, nrow = 1) + 
  scale_fill_manual(values=c("red", "gold", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g]))) + 
  scale_color_manual(values=c("red", "darkgoldenrod3", "green"), labels = c(expression(R[r]),expression(R[y]),expression(R[g])))+ 
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text.x = element_text(angle = 0, hjust = 0),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  ) + 
  xlim(0,0.25)

p2 <- 
  SixParameterData %>%
  filter(Parameter_type == "Motility")  %>% 
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

ggsave("SixParameters_SIM_CellDensity.pdf", width = 9, height = 8, units = "in")

#Proliferation Plots
p1 <- 
  ProliferationData %>% 
  gather(Parameter,Value,-Simulation) %>% 
  mutate(TrueValue = case_when(
    Simulation == 1 & Parameter == "Rr" ~ 0.04,
    Simulation == 1 & Parameter == "Ry" ~ 0.17,
    Simulation == 1 & Parameter == "Rg" ~ 0.08,
    Simulation == 2 & Parameter == "Rr" ~ 0.25,
    Simulation == 2 & Parameter == "Ry" ~ 0.15,
    Simulation == 2 & Parameter == "Rg" ~ 0.22,
    Simulation == 3 & Parameter == "Rr" ~ 0.12,
    Simulation == 3 & Parameter == "Ry" ~ 0.07,
    Simulation == 3 & Parameter == "Rg" ~ 0.03,
    Simulation == 4 & Parameter == "Rr" ~ 0.3,
    Simulation == 4 & Parameter == "Ry" ~ 0.36,
    Simulation == 4 & Parameter == "Rg" ~ 0.28
    )
  ) %>% 
  ggplot(aes(x = Value, fill = Parameter, linetype = Parameter)) + 
  geom_density(alpha = 0.5, size = 0.75) +
  geom_vline(aes(xintercept = TrueValue, color = Parameter), linetype='dashed') + 
  facet_wrap(~Simulation, nrow = 1) + 
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

ggsave("Proliferation_SIM.pdf", width = 5, height = 5, units = "in")


#Motility Plots 
MotilityData <-
  MotilityData %>% 
  gather(Parameter, Value, -c(Simulation, ntrack, MotilityData)) %>% 
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
    Simulation == 4 & Parameter == "Mg" ~ 2 
  )) %>% 
  mutate(facets = case_when(
    Simulation == 1 ~ "(a)",
    Simulation == 2 ~ "(b)",
    Simulation == 3 ~ "(c)",
    Simulation == 4 ~ "(d)"
  )) 


MotilityData %>% 
  filter(MotilityData == "CellTracking") %>% 
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

ggsave("Motility_SIM_CellTracking.pdf", width = 9, height = 4, units = "in")

MotilityData %>% 
  filter(MotilityData == "CellDensity") %>% 
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

ggsave("Motility_SIM_CellDensity.pdf", width = 9, height = 4, units = "in")

## Cell Trajectories ## - 8x8

CellTracking %>% 
  mutate(firstID = ifelse(frame == 1, ntrack, NA)) %>%
  mutate(
    color = case_when(color == 1 ~ "G1", color == 2 ~ "eS", color == 3 ~ "S/G2/M"),
    color = factor(color, levels = c("G1", "eS", "S/G2/M"))
    ) %>% 
  mutate(facets = case_when(
    ntrack == 1 ~ "(a)",
    ntrack == 2 ~ "(b)",
    ntrack == 3 ~ "(c)",
    ntrack == 4 ~ "(d)",
    ntrack == 5 ~ "(e)",
    ntrack == 6 ~ "(f)",
    ntrack == 7 ~ "(g)",
    ntrack == 8 ~ "(h)",
    ntrack == 9 ~ "(i)",
    ntrack == 10 ~ "(j)",
    ntrack == 11 ~ "(k)",
    ntrack == 12 ~ "(l)",
    ntrack == 13 ~ "(m)",
    ntrack == 14 ~ "(n)",
    ntrack == 15 ~ "(o)",
    ntrack == 16 ~ "(p)",
    ntrack == 17 ~ "(q)",
    ntrack == 18 ~ "(r)",
    ntrack == 19 ~ "(s)",
    ntrack == 20 ~ "(t)")) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(color), shape = as.factor(color))) +
  geom_path(color = "black") +
  facet_wrap(~facets, scale = "free", nrow = 5) +
  labs(x = expression(paste("x ",mu,"m")), y = expression(paste("y ",mu,"m"))) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0, hjust = 0)
  ) +
  scale_colour_manual(values = c("red", "gold", "green")) 

ggsave("CellTrajectories.pdf", width = 8, height = 8, units = "in")

#cell trajectories (sub sample)
CellTracking %>% 
  mutate(firstID = ifelse(frame == 1, ntrack, NA)) %>% 
  mutate(
    color = case_when(color == 1 ~ "G1", color == 2 ~ "eS", color == 3 ~ "S/G2/M"),
    color = factor(color, levels = c("G1", "eS", "S/G2/M"))
  ) %>%  
  filter(ntrack %in% c(1,12,11,15,2,9,16,18,3,17,5,14,19)) %>% 
  ggplot(aes(x = x, y = y, group = ntrack)) +
  geom_point(aes(color = as.factor(color), shape = as.factor(color))) +
  geom_path(color = "black") +
  labs(x = "", y = "") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_colour_manual(values = c("red", "gold", "green"))

ggsave("CellTrajectories_subsample.pdf", width = 5, height = 5, units = "in")
