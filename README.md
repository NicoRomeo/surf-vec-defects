# surf-vec-defects
Locates topological defects in tangent vector fields to a closed surface.

The code has been developed for usage on vector fields living on spheres, sampled on a (quasi-)regular point grid. However, the code should be able to work on closed surfaces with no overhangs, eg ellipsoids or weakly perturbed spheres.


the function DefectTrackSurf takes 2 mandatory arguments and a third optional:
- the array of points r of size (# of points)x3
- the array of sampled vectors v of size (# of points)x3x(# of time points)
- the distance dr between nearest neighbors points. If not provided, dr defaults to sin(2pi/sqrt(# grid points)), the approximate spacing between points on a regular grid on the sphere.


The algorithm is detailed [here](https://drive.google.com/file/d/1s_5GXzcphbt8bVqCTsxOaLMuaihzjxvJ/view?usp=sharing).
