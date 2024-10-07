# MultipleChangePoints
Simulate a Gaussian time series with several changes in the mean (changepoints) and identify them

Compile the programs with
```
gfortran kind.f90 random.f90 change_point_util.f90 split_data.f90 xxsplit_data.f90
gfortran kind.f90 random.f90 change_point_util.f90 split_data.f90 xoptimize_change_points.f90
gfortran change_point_util.f90 xneighbors.f90
```
The output of running `xxsplit_data.f90` is below. The Estimated Change Points are usually close to
the True Change Points.
```
 #obs:        1000
 True Change Points:         100         200         700
      Segment Means: 5.0000 6.0000 7.0000 8.0000

 Estimated Change Points:         100         197         697
 Total SSE:   1036.1501559502469     

 Estimated Change Points:          89         194         706
 Total SSE:   1018.7717060194260     

 Estimated Change Points:         101         199         695
 Total SSE:   1008.4290460992340     

 Estimated Change Points:          94         200         706
 Total SSE:   1027.5405334996076     

 Estimated Change Points:          99         195         699
 Total SSE:   910.42981169375798
```    
