#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 

Rcpp::List GenerateSummaryStatData(bool &CellTracking, Rcpp::List &SummaryStatsData, int &ntrack, double &T_record, double &t, int &Nred, int &Nyellow, int &Ngreen, ivec &RowPosCell, ivec &ColPosCell, ivec &CellSelected, ivec &CellSelectedStart, ivec &MigPosition, int &rowIndex, int &columnIndex, bool &migFailed, int &transID, double &delta, int &simuIndex, int &simuNum, imat &domain, mat &domain_x) {
  
  if (CellTracking) {
  
    if (SummaryStatsData.length() == 0) { //initialise list
      //summary stats
      ivec Nred_ss(simuNum); //count of red cells at T_record
      ivec Nyellow_ss(simuNum); //count of yellow cells at T_record
      ivec Ngreen_ss(simuNum); //count of green cells at T_record
      Rcpp::NumericMatrix red_distance(simuNum,ntrack); //distance traveled in red phase of tracked cells
      Rcpp::NumericMatrix yellow_distance(simuNum,ntrack); //distance traveled in yellow phase of tracked cells
      Rcpp::NumericMatrix green_distance(simuNum,ntrack); //distance traveled in green phase of tracked cells
      
      //additional helper variables
      Rcpp::NumericMatrix red_time(simuNum,ntrack); //time spent in red phase for tracked cells
      Rcpp::NumericMatrix yellow_time(simuNum,ntrack); //time spent in yellow phase for tracked cells
      Rcpp::NumericMatrix green_time(simuNum,ntrack); //time spent in green phase for tracked cells
      mat phasestart_time(simuNum,ntrack); //time which previous phases ended - initially 0
      imat TrackingCompleted(simuNum,ntrack,fill::zeros); //Indicates whether cell cycle has be completed
      
      SummaryStatsData = Rcpp::List::create(
        Rcpp::Named("Nred_ss") = Nred_ss,
        Rcpp::Named("Nyellow_ss") = Nyellow_ss,
        Rcpp::Named("Ngreen_ss") = Ngreen_ss, 
        Rcpp::Named("red_distance") = red_distance,
        Rcpp::Named("yellow_distance") = yellow_distance,
        Rcpp::Named("green_distance") = green_distance,
        Rcpp::Named("red_time") = red_time,
        Rcpp::Named("yellow_time") = yellow_time,
        Rcpp::Named("green_time") = green_time,
        Rcpp::Named("phasestart_time") = phasestart_time,
        Rcpp::Named("TrackingCompleted") = TrackingCompleted
      );
    } 
    
    //summary stats
    ivec Nred_ss = SummaryStatsData["Nred_ss"];
    ivec Nyellow_ss = SummaryStatsData["Nyellow_ss"];
    ivec Ngreen_ss = SummaryStatsData["Ngreen_ss"];
    Rcpp::NumericMatrix red_distance = SummaryStatsData["red_distance"];
    Rcpp::NumericMatrix yellow_distance = SummaryStatsData["yellow_distance"];
    Rcpp::NumericMatrix green_distance = SummaryStatsData["green_distance"];
    
    //additional helper variables
    Rcpp::NumericMatrix red_time = SummaryStatsData["red_time"];
    Rcpp::NumericMatrix yellow_time = SummaryStatsData["yellow_time"];
    Rcpp::NumericMatrix green_time = SummaryStatsData["green_time"];
    mat phasestart_time = SummaryStatsData["phasestart_time"];
    imat TrackingCompleted = SummaryStatsData["TrackingCompleted"];
    
    //Record Summary Statistics
    if (t >= T_record - 1) {
      Nred_ss(simuIndex - 1) = Nred;
      Nyellow_ss(simuIndex - 1) = Nyellow;
      Ngreen_ss(simuIndex - 1) = Ngreen;
    }
  
    int cellID = -1; //ID of tracked cells 0:ntrack-1
    if (!migFailed || transID != 0) {
      for (int i = 0; i < ntrack; i++) {    
        if (RowPosCell(i) == rowIndex && ColPosCell(i) == columnIndex) { //determine which tracked cell moved/transitioned
          cellID = i;
          break;
        }
      }
      if (cellID != -1) {
        if (transID != 0) { //record phase change
  
          if (TrackingCompleted(simuIndex - 1, cellID) == 0) {
            if (transID == 1) {
              green_time(simuIndex - 1, cellID) = t - phasestart_time(simuIndex - 1, cellID); //calculate phase duration 
              if (red_time(simuIndex - 1, cellID) == 0) { //ensure not overiding first phase
                red_time(simuIndex - 1, cellID) = T_record - t; //assign next phase duration based on termination time - to be overriden
                phasestart_time(simuIndex - 1, cellID) = t;
                CellSelected(cellID) = transID;
              }
            } else if (transID == 2) {
              red_time(simuIndex - 1, cellID) = t - phasestart_time(simuIndex - 1, cellID); //calculate phase duration 
              if (yellow_time(simuIndex - 1, cellID) == 0) { //ensure not overiding first phase
                yellow_time(simuIndex - 1, cellID) = T_record - t; //assign next phase duration based on termination time - to be overriden
                phasestart_time(simuIndex - 1, cellID) = t;
                CellSelected(cellID) = transID;
              }
            } else if (transID == 3) {
              yellow_time(simuIndex - 1, cellID) = t - phasestart_time(simuIndex - 1, cellID); //calculate phase duration 
              if (green_time(simuIndex - 1, cellID) == 0) { //ensure not overiding first phase
                green_time(simuIndex - 1, cellID) = T_record - t; //assign next phase duration based on termination time - to be overriden
                phasestart_time(simuIndex - 1, cellID) = t;
                CellSelected(cellID) = transID;
              }
            }
          }
  
          if (CellSelected(cellID) == CellSelectedStart(cellID)){
            TrackingCompleted(simuIndex - 1, cellID) = 1;
          }
  
        } else if (!migFailed) { //record distance traveled
          if (TrackingCompleted(simuIndex - 1, cellID) == 0){
            if (CellSelected(cellID) == 1) { //if cell is in red phase
              red_distance(simuIndex - 1, cellID) += delta;
            } else if (CellSelected(cellID) == 2) { //if cell is in yellow phase
              yellow_distance(simuIndex - 1, cellID) += delta;
            } else if (CellSelected(cellID) == 3) { //if cell is in green phase
              green_distance(simuIndex - 1, cellID) += delta;
            }
          }
        }
      }
    }
    
    //summary stats
    SummaryStatsData["Nred_ss"] = Nred_ss;
    SummaryStatsData["Nyellow_ss"] = Nyellow_ss;
    SummaryStatsData["Ngreen_ss"] = Ngreen_ss;
    SummaryStatsData["red_distance"] = red_distance;
    SummaryStatsData["yellow_distance"] = yellow_distance;
    SummaryStatsData["green_distance"] = green_distance;
    
    //additional helper variables
    SummaryStatsData["red_time"] = red_time;
    SummaryStatsData["yellow_time"] = yellow_time;
    SummaryStatsData["green_time"] = green_time;
    SummaryStatsData["phasestart_time"] = phasestart_time;
    SummaryStatsData["TrackingCompleted"] = TrackingCompleted;
  
  } else {

    //Record Summary Statistics
    if (t >= T_record - 1) {
      
      if (SummaryStatsData.length() == 0) {
      
        ivec Nred_ss(simuNum); //count of red cells at T_record
        ivec Nyellow_ss(simuNum); //count of yellow cells at T_record
        ivec Ngreen_ss(simuNum); //count of green cells at T_record
        
        Rcpp::List red_position(simuNum); //list with x positions of all red cels at T_record
        Rcpp::List yellow_position(simuNum); //list with x positions of all red cels at T_record
        Rcpp::List green_position(simuNum); //list with x positions of all red cels at T_record
        
        SummaryStatsData = Rcpp::List::create(
          Rcpp::Named("Nred_ss") = Nred_ss,
          Rcpp::Named("Nyellow_ss") = Nyellow_ss,
          Rcpp::Named("Ngreen_ss") = Ngreen_ss, 
          Rcpp::Named("red_position") = red_position,
          Rcpp::Named("yellow_position") = yellow_position,
          Rcpp::Named("green_position") = green_position
        );
      
      }
      
      ivec Nred_ss = SummaryStatsData["Nred_ss"];
      ivec Nyellow_ss = SummaryStatsData["Nyellow_ss"];
      ivec Ngreen_ss = SummaryStatsData["Ngreen_ss"];
      
      Rcpp::List red_position  = SummaryStatsData["red_position"];
      Rcpp::List yellow_position = SummaryStatsData["yellow_position"];
      Rcpp::List green_position = SummaryStatsData["green_position"]; 
      
      Nred_ss(simuIndex - 1) = Nred;
      Nyellow_ss(simuIndex - 1) = Nyellow;
      Ngreen_ss(simuIndex - 1) = Ngreen;
      
      red_position(simuIndex - 1) = (vec) domain_x.elem(find(domain == 1));
      yellow_position(simuIndex - 1) = (vec)  domain_x.elem(find(domain == 2));
      green_position(simuIndex - 1) = (vec) domain_x.elem(find(domain == 3));
        
      SummaryStatsData = Rcpp::List::create(
        Rcpp::Named("Nred_ss") = Nred_ss,
        Rcpp::Named("Nyellow_ss") = Nyellow_ss,
        Rcpp::Named("Ngreen_ss") = Ngreen_ss, 
        Rcpp::Named("red_position") = red_position,
        Rcpp::Named("yellow_position") = yellow_position,
        Rcpp::Named("green_position") = green_position
      );
      
    }
  }
  
  return SummaryStatsData;
}
