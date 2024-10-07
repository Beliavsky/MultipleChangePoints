# MultipleChangePoints
Simulate a Gaussian time series with several changes in the mean (changepoints) and identify them

Compile the programs with
```
gfortran kind.f90 random.f90 change_point_util.f90 split_data.f90 xxsplit_data.f90
gfortran kind.f90 random.f90 change_point_util.f90 split_data.f90 xoptimize_change_points.f90
gfortran change_point_util.f90 xneighbors.f90
```
