library(ggplot2)
library(dplyr)
library(readxl)
library(writexl)
library(magick)
library(RSAGA)

#set working directory to the directory of this R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

N <- 110 #number of images

##-- Convert Pixels to microns --##
#Pixels: 960x720
#Microns: 1309.09x1745.35
x_scale <- 1309.09/720
y_scale <- 1745.35/960


img_df2 <- NULL 
for (i in 1:N) {
  
  if (i < 10) {
    num <- paste0("00",i)
  } else if (i < 100) {
    num <- paste0("0",i)
  } else {
    num <- i
  }
  
  #read image
  img <- image_read(paste0("../Raw Data/scene00",num,".png")) #read image from parent folder (Data)
  img_df <- grid.to.xyz(as.matrix(as.raster(img)))
  img_df <- img_df %>% filter(z != "#000000ff") #remove black background
  
  #calculate rgb decimal codes
  tmp <- as.matrix(t(col2rgb(as.character(img_df$z))))
  img_df <- cbind(img_df,tmp)
  
  img_df <- img_df %>% 
    mutate(color = case_when(
    red > 100 & green <= 100 ~ 1, #red
    green > 100 & red <= 100 ~ 3, #green
    red > 100 & green > 100 ~ 2, #yellow
    red < 50 & green < 50 & blue < 50 ~ 0, #black
  )) %>% 
    mutate(color = case_when( 
      is.na(color) & red > green ~ 1, #red
      is.na(color) & red < green ~ 3, #green
      TRUE ~ color
    )) %>% 
    filter(color != 0) %>% #remove black pixels
    mutate(frame = i) 
  
  img_df2 <- bind_rows(img_df2,img_df)
}

FUCCI_Data <- img_df2
InitPos <- read_excel("../Raw Data/FUCCI_DATA_ImageJ.xlsx", sheet = "InitPos")  #read ImageJ data from parent folder (Data)
FinalPos <- read_excel("../Raw Data/FUCCI_DATA_ImageJ.xlsx", sheet = "FinalPos")  #read ImageJ data from parent folder (Data)
CellTracking <- read_excel("../Raw Data/FUCCI_DATA_ImageJ.xlsx", sheet = "CellTracking")  #read ImageJ data from parent folder (Data)

InitPos_color <- 
  InitPos %>% 
  as_tibble() %>% 
  mutate(frame = 1) %>% 
  mutate(x = round(x), y = round(y)) %>% 
  left_join(.,FUCCI_Data, by = c("x", "y", "frame")) %>% 
  select(x,y, color) %>% 
  filter(!is.na(color)) #remove unidentified cells

FinalPos_color <- 
  FinalPos %>% 
  as_tibble() %>% 
  mutate(frame = 110) %>% 
  mutate(x = round(x), y = round(y)) %>% 
  left_join(.,FUCCI_Data, by = c("x", "y", "frame")) %>% 
  select(x,y, color) %>% 
  filter(!is.na(color)) #remove unidentified cells 

CellTracking_color <- 
  CellTracking %>% 
  as_tibble() %>% 
  mutate(x = round(x), y = round(y)) %>% 
  left_join(.,FUCCI_Data, by = c("x", "y", "frame")) %>% 
  select(x,y, color, ntrack, frame) 

## cleaning data
for (i in 1:nrow(CellTracking_color)){
  if (CellTracking_color[i,"frame"] != 1){
    if (is.na(CellTracking_color[i,"color"])){
      CellTracking_color[i,"color"] <- CellTracking_color[i-1,"color"]
    }
    if (CellTracking_color[i,"color"] < CellTracking_color[i-1,"color"]) {
      CellTracking_color[i,"color"] <- CellTracking_color[i-1,"color"]
    }
  }
}

InitPos_color <- InitPos_color %>% mutate(x = x*x_scale, y = y*y_scale)
FinalPos_color <- FinalPos_color %>% mutate(x = x*x_scale, y = y*y_scale)
CellTracking_color <- CellTracking_color %>% mutate(x = x*x_scale, y = y*y_scale)

write_xlsx(list(InitPos = InitPos_color, FinalPos = FinalPos_color, CellTracking = CellTracking_color), "FUCCI_processed.xlsx")
