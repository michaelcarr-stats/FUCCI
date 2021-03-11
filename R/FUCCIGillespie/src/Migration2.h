#pragma once

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

ivec Migration2(int& rowIndex, int& columnIndex, int& rowNum, int& columnNum, int& BC) {
    int MigRow = rowIndex;
    int MigColumn = columnIndex;
    double P = R::runif(0,1);
    
        /*
    Target Sites arround lattice site L.
       ___  ___                  _______
      / 3 \/ 5 \                | 3 | 5 |
      \___/\___/__           ___|___|___|
    / 1 \/ L \/ 2 \    =>   | 1 | L | 2 |
    \___/\___/\___/         |___|___|___|
      / 4 \/ 6 \                | 4 | 6 |
      \___/\___/                |___|___|

    */

    if (BC == 1) {
        if (P < 1.0 / 6.0) { //Target site 1
            if (columnIndex == 0) {
                MigColumn = columnNum - 1;
            }
            else {
                MigColumn = columnIndex - 1;
            }
        }
        else if (P >= 1.0 / 6.0 && P <= 2.0 / 6.0) { //Target site 2
            if (columnIndex == columnNum - 1) {
                MigColumn = 0;
            }
            else {
                MigColumn = columnIndex + 1;
            }
        }
        else if (P >= 2.0 / 6.0 && P <= 3.0 / 6.0) { //Target site 3
            if (rowIndex == 0) {
                MigRow = rowNum - 1;
            }
            else {
                MigRow = rowIndex - 1;
            }
        }
        else if (P >= 3.0 / 6.0 && P <= 4.0 / 6.0) { //Target site 4
            if (rowIndex == rowNum - 1) {
                MigRow = 0;
            }
            else {
                MigRow = rowIndex + 1;
            }
        }
        else if (P >= 4.0 / 6.0 && P <= 5.0 / 6.0) { //Target site 5
            if (rowIndex == 0) {
                MigRow = rowNum - 1;
            }
            else {
                MigRow = rowIndex - 1;
            }
            if (columnIndex == columnNum - 1) {
                MigColumn = 0;
            }
            else {
                MigColumn = columnIndex + 1;
            }
        }
        else if (P >= 5.0 / 6.0) { //Target site 6
            if (rowIndex == rowNum - 1) {
                MigRow = 0;
            }
            else {
                MigRow = rowIndex + 1;
            }
            if (columnIndex == 0) {
                MigColumn = columnNum - 1;
            }
            else {
                MigColumn = columnIndex + 1;
            }
        }
    }
    else if (BC == 2) {
        if (P < 1.0 / 6.0) { //Target site 1
            if (columnIndex == 0) {
                MigColumn = 0;
            }
            else {
                MigColumn = columnIndex - 1;
            }
        }
        else if (P >= 1.0 / 6.0 && P <= 2.0 / 6.0) { //Target site 2
            if (columnIndex == columnNum - 1) {
                MigColumn = columnNum - 1;
            }
            else {
                MigColumn = columnIndex + 1;
            }
        }
        else if (P >= 2.0 / 6.0 && P <= 3.0 / 6.0) { //Target site 3
            if (rowIndex == 0) {
                MigRow = 0;
            }
            else {
                MigRow = rowIndex - 1;
            }
        }
        else if (P >= 3.0 / 6.0 && P <= 4.0 / 6.0) { //Target site 4
            if (rowIndex == rowNum - 1) {
                MigRow = rowNum - 1;
            }
            else {
                MigRow = rowIndex + 1;
            }
        }
        else if (P >= 4.0 / 6.0 && P <= 5.0 / 6.0) { //Target site 5
            if (rowIndex == 0) {
                MigRow = 0;
            }
            else {
                MigRow = rowIndex - 1;
            }
            if (columnIndex == columnNum - 1) {
                MigColumn = columnNum - 1;
            }
            else {
                MigColumn = columnIndex + 1;
            }
        }
        else if (P >= 5.0 / 6.0) { //Target site 6
            if (rowIndex == rowNum - 1) {
                MigRow = rowNum - 1;
            }
            else {
                MigRow = rowIndex + 1;
            }
            if (columnIndex == columnNum - 1) {
                MigColumn = columnNum - 1;
            }
            else {
                MigColumn = columnIndex + 1;
            }
        }
    }
    ivec MigPosition = {MigRow, MigColumn};
    
    return MigPosition;
}
