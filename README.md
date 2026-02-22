# CONTROL-SURVEY-PROJECT-GROUP-14
semi-graphical fix for computation of final co-ordinate
# Surveying Computation Suite

A comprehensive Python toolkit for surveying computations including coordinate adjustment, error analysis, and network adjustment using rigorous least squares methods.

## Overview

This suite implements classical surveying computation methods with modern numerical techniques. It handles coordinate determination from angular and distance observations, performs rigorous error analysis, and provides professional visualization of results.

## Features

### 1. Provisional Coordinate Computation
- **Intersection Method**: Compute unknown point from two known stations with bearings
- **Resection Method**: Compute position from angles observed at unknown point
- **Combined Angles**: Use observed angles α and β at known stations
- **Direct Bearings**: Input observed bearings directly

### 2. Cut Computation (Coordinate Adjustment)
- **Original Method**: Classical coordinate cuts using cotangent and tangent formulas
- **Perpendicular Projection**: Shortest distance to bearing line (alternative method)
- **Multiple Stations**: Handle redundant observations from 3+ stations
- **Error Detection**: Identify and flag inconsistent observations

### 3. Least Squares Adjustment
- **Rigorous Adjustment**: Full least squares solution for overdetermined systems
- **Weight Matrix**: Distance-based weighting (optional)
- **Residual Analysis**: Compute and display observation residuals
- **Precision Estimates**: Standard deviations of adjusted coordinates

### 4. Visualization & Reporting
- **Cartesian Plots**: Cut lines with intersection point
- 
- **Residual Plots**: Analysis of observation quality
- **Professional Tables**: Formatted output matching field book standards



