#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


//Functions
#include "Migration1.h"
#include "Migration2.h"
#include "GenerateSummaryStatData.h"

//[[Rcpp::export]]
Rcpp::List Main_Simulate(arma::vec theta, Rcpp::List SetupVars, double T_record, bool CellTracking) {
  
  //////////////////////////////////////////////
  //		      	Variables/Parmaters				    //
  //////////////////////////////////////////////						
  
  //Input variables
  double Kry = theta(0); // transition rate of red to yellow per time unit
  double Kyg = theta(1); // transition rate of yellow to green per time unit
  double Kgr = theta(2); // transition rate of green to red per time unit
  double Pr = theta(3); // migration rate of red per time unit
  double Py = theta(4); // migration rate of yellow per time unit
  double Pg = theta(5); // migration rate of green per time unit
  
  int ntrack = Rcpp::as<int>(SetupVars["ntrack"]);
  int simuNum = Rcpp::as<int>(SetupVars["simuNum"]);
  double delta = Rcpp::as<double>(SetupVars["delta"]);
  double t = Rcpp::as<double>(SetupVars["t"]);
  double tstop = Rcpp::as<double>(SetupVars["tstop"]);
  double Xmax = Rcpp::as<double>(SetupVars["Xmax"]);
  double Ymax = Rcpp::as<double>(SetupVars["Ymax"]);
  int columnNum = Rcpp::as<int>(SetupVars["columnNum"]);
  int rowNum = Rcpp::as<int>(SetupVars["rowNum"]);
  int BC = Rcpp::as<int>(SetupVars["BC"]);
  
  imat domain = Rcpp::as<imat>(SetupVars["domain"]);
  mat domain_x = Rcpp::as<mat>(SetupVars["domain_x"]);
  mat domain_y = Rcpp::as<mat>(SetupVars["domain_y"]);
  ivec NStartCount = Rcpp::as<ivec>(SetupVars["NStartCount"]);
  ivec RowPosCell = Rcpp::as<ivec>(SetupVars["RowPosCell"]);
  ivec ColPosCell = Rcpp::as<ivec>(SetupVars["ColPosCell"]);
  ivec CellSelected = Rcpp::as<ivec>(SetupVars["CellSelected"]);
  ivec CellSelectedStart = CellSelected;
  
  int Nred = (int) NStartCount(0);
  int Nyellow = (int) NStartCount(1);
  int Ngreen = (int) NStartCount(2);
  
  //local variables
  double ar; // red movement
  double ay; // yellow movement
  double ag; // green movement
  double tr; // red transition
  double ty; // yellow transition
  double tg; // green transition
  double a0; // total propensity function
  double tau; // time step
  bool ExitSimStatus = 1; // indicate if simulation was terminated early  
  int simuIndex = 1; //The simulation index starting from 1
  double R; //probability determining proliferation or motility event
  
  int rowIndex; //row index of selected cell
  int columnIndex; //column index of selected cell
  ivec MigPosition(2); // migration row and column position
  uvec MatchedIndicies; // vectot of found indicies (singular)
  int Index; // vector of indicies in the vector MatchedIndicies
  
  //output objects
  Rcpp::List L; 
  
  Rcpp::List SummaryStatsData; //inititialse empty list
  
  //////////////////////////////////////////////
  //				      Main Framework				      //
  //////////////////////////////////////////////	
  
  while (simuIndex <= simuNum) {
    
    bool migFailed = 1; //indicator if chosen site is not vacant
    int transID = 0; //ID of phase transitioning to: green->red - 1, red->yellow - 2, yellow->green - 3 
    
    while (t < tstop) {
      
      ar = Pr * Nred; // red movement
      ay = Py * Nyellow; // yellow movement
      ag = Pg * Ngreen; // green movement
      tr = Kry * Nred; // red transition
      ty = Kyg * Nyellow; // yellow transition
      tg = Kgr * Ngreen; // green transition
      a0 = ar + ay + ag + tr + ty + tg; // total propensity function
      
      R = R::runif(0,1);

      //red Migration  
      if (R <= ar/a0) {
        MatchedIndicies = find(domain == 1);
        Index = Rcpp::as<int>(Rcpp::sample(MatchedIndicies.n_elem,1,false))-1; // chosen cell row index
        rowIndex = MatchedIndicies(Index) % rowNum;
        columnIndex = floor(MatchedIndicies(Index)/rowNum);
        if ((rowNum - rowIndex) % 2 == 0) { //If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
          MigPosition = Migration1(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 0; //Remove agent at prevous site
            domain(MigPosition(0), MigPosition(1)) = 1; //Move agent to the new site
            migFailed = 0;
          }
        } else { //If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
          MigPosition = Migration2(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 0; //Remove agent at prevous site
            domain(MigPosition(0), MigPosition(1)) = 1; //Move agent to the new site
            migFailed = 0;
          }
        }
        
        //Yellow Migration
      } else if (R > ar/a0 && R <= (ay + ar)/a0) {
        MatchedIndicies = find(domain == 2);
        Index = Rcpp::as<int>(Rcpp::sample(MatchedIndicies.n_elem,1,false))-1; // chosen cell row index
        rowIndex = MatchedIndicies(Index) % rowNum;
        columnIndex = floor(MatchedIndicies(Index)/rowNum);
        if ((rowNum - rowIndex) % 2 == 0) { //If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
          MigPosition = Migration1(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 0; //Remove agent at prevous site
            domain(MigPosition(0), MigPosition(1)) = 2; //Move agent to the new site
            migFailed = 0;
          }
        } else { //If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
          MigPosition = Migration2(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 0; //Remove agent at prevous site
            domain(MigPosition(0), MigPosition(1)) = 2; //Move agent to the new site
            migFailed = 0;
          }
        }
        
        //Green Migration
      } else if (R > (ar + ay)/a0 && R <= (ar + ay + ag)/a0) {
        MatchedIndicies = find(domain == 3);
        Index = Rcpp::as<int>(Rcpp::sample(MatchedIndicies.n_elem,1,false))-1; // chosen cell row index
        rowIndex = MatchedIndicies(Index) % rowNum;
        columnIndex = floor(MatchedIndicies(Index)/rowNum);
        if ((rowNum - rowIndex) % 2 == 0) { //If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
          MigPosition = Migration1(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 0; //Remove agent at previous site
            domain(MigPosition(0), MigPosition(1)) = 3; //Move agent to the new site
            migFailed = 0;
          }
        } else { //If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
          MigPosition = Migration2(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 0; //Remove agent at previous site
            domain(MigPosition(0), MigPosition(1)) = 3; //Move agent to the new site
            migFailed = 0;
          }
        }
        
        //Transition Red to Yellow
      } else if (R > (ar + ay + ag)/a0 && R <= (ar + ay + ag + tr)/a0) {
        MatchedIndicies = find(domain == 1);
        Index = Rcpp::as<int>(Rcpp::sample(MatchedIndicies.n_elem,1,false))-1; // chosen cell row index
        rowIndex = MatchedIndicies(Index) % rowNum;
        columnIndex = floor(MatchedIndicies(Index)/rowNum);
        domain(rowIndex, columnIndex) = 2; //change domain entry to a yellow agent
        Nyellow += 1; //Update yellow population
        Nred -= 1; //Update Red population
        transID = 2;
        
        
        //Transition Yellow to Green
      } else if (R > (ar + ay + ag + tr)/a0 && R <= (ar + ay + ag + tr + ty)/a0) {
        MatchedIndicies = find(domain == 2);
        Index = Rcpp::as<int>(Rcpp::sample(MatchedIndicies.n_elem,1,false)) - 1; // chosen cell index
        rowIndex = MatchedIndicies(Index) % rowNum;
        columnIndex = floor(MatchedIndicies(Index)/rowNum);
        domain(rowIndex, columnIndex) = 3; //change domain entry to a yellow agent
        Ngreen += 1; //Update green population
        Nyellow -= 1; //Update yellow population
        transID = 3;
        
        
        //Transition Green to Red
      } else if (R > (ar + ay + ag + tr + ty)/ a0 && R <= (ar + ay + ag + tr + ty + tg)/a0) {
        MatchedIndicies = find(domain == 3);
        Index = Rcpp::as<int>(Rcpp::sample(MatchedIndicies.n_elem,1,false))-1; // chosen cell row index
        rowIndex = MatchedIndicies(Index) % rowNum;
        columnIndex = floor(MatchedIndicies(Index)/rowNum);
        if ((rowNum - rowIndex) % 2 == 0) { //If the selected agent is at bottom/(bottom-2)/(bottom-4)/... rows
          MigPosition = Migration1(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 1; 
            domain(MigPosition(0), MigPosition(1)) = 1; //Move agent to the new site
            Nred += 2; //update Red cell population
            Ngreen -= 1; //update Green cell population
            transID = 1;
          }
        }
        else { //If the selected agent is at (bottom-1)/(bottom-3)/(bottom-5)/... rows
          MigPosition = Migration2(rowIndex, columnIndex, rowNum, columnNum, BC);
          if (domain(MigPosition(0), MigPosition(1)) == 0 && domain_x(MigPosition(0), MigPosition(1)) <= Xmax && domain_y(MigPosition(0), MigPosition(1)) <= Ymax) {
            domain(rowIndex, columnIndex) = 1; 
            domain(MigPosition(0), MigPosition(1)) = 1; //Move agent to the new site
            Nred += 2; //update Red cell population
            Ngreen -= 1; //update Green cell population
            transID = 1;
          }
        }
      }
      
      SummaryStatsData = GenerateSummaryStatData(CellTracking, SummaryStatsData, ntrack, T_record, t, Nred, Nyellow, Ngreen, RowPosCell, ColPosCell, CellSelected, CellSelectedStart, rowIndex, columnIndex, migFailed, transID, delta, simuIndex, simuNum, domain, domain_x);
      
      if (!migFailed && (int) SummaryStatsData["cellID"] != -1) {
        RowPosCell((int) SummaryStatsData["cellID"]) = MigPosition(0); //update location of tracked cells
        ColPosCell((int) SummaryStatsData["cellID"]) = MigPosition(1); //update location of tracked cells
      } else if (transID != 0 && (int) SummaryStatsData["cellID"] != -1) {
        CellSelected((int) SummaryStatsData["cellID"]) = transID //update cell ID
      }
      
      //reset transition/migration indicators
      transID = 0;
      migFailed = 1;
      
      
      //	calculate time step
      tau = (1 / a0) * log(1 / R::runif(0,1));
      t += tau;
    }
    
    
    if (simuNum != simuIndex) {
      
      //Re-initialse simulation
      t = Rcpp::as<double>(SetupVars["t"]);
      domain = Rcpp::as<imat>(SetupVars["domain"]);
      domain_x = Rcpp::as<mat>(SetupVars["domain_x"]);
      domain_y = Rcpp::as<mat>(SetupVars["domain_y"]);
      NStartCount = Rcpp::as<ivec>(SetupVars["NStartCount"]);
      RowPosCell = Rcpp::as<ivec>(SetupVars["RowPosCell"]);
      ColPosCell = Rcpp::as<ivec>(SetupVars["ColPosCell"]);
      CellSelected = Rcpp::as<ivec>(SetupVars["CellSelected"]);
      Nred = (int) NStartCount(0);
      Nyellow = (int) NStartCount(1);
      Ngreen = (int) NStartCount(2);
      
    }
    
    simuIndex += 1;
  }
  
  ExitSimStatus = 0; //indicate simulation finished successfully
  
  if (CellTracking) {
    return L = Rcpp::List::create(Rcpp::Named("ExitSimStatus") = Rcpp::wrap(ExitSimStatus), Rcpp::Named("Nred") = SummaryStatsData["Nred_ss"], Rcpp::Named("Nyellow") = SummaryStatsData["Nyellow_ss"], Rcpp::Named("Ngreen") = SummaryStatsData["Ngreen_ss"], Rcpp::Named("red_distance") = SummaryStatsData["red_distance"], Rcpp::Named("yellow_distance") = SummaryStatsData["yellow_distance"], Rcpp::Named("green_distance") = SummaryStatsData["green_distance"], Rcpp::Named("red_time") = SummaryStatsData["red_time"], Rcpp::Named("yellow_time") = SummaryStatsData["yellow_time"], Rcpp::Named("green_time") = SummaryStatsData["green_time"]);
  } else { //spatio-temporal
    return L = Rcpp::List::create(Rcpp::Named("ExitSimStatus") = Rcpp::wrap(ExitSimStatus), Rcpp::Named("Nred") = SummaryStatsData["Nred_ss"], Rcpp::Named("Nyellow") = SummaryStatsData["Nyellow_ss"], Rcpp::Named("Ngreen") = SummaryStatsData["Ngreen_ss"], Rcpp::Named("red_position") = SummaryStatsData["red_position"], Rcpp::Named("yellow_position") = SummaryStatsData["yellow_position"], Rcpp::Named("green_position") = SummaryStatsData["green_position"]);
  }
}
