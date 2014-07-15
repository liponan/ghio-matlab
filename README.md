ghio-matlab
===========

Guided hybrid input-and-output (GHIO), HIO and shrink-wrap Matlab functions

# Introduction 

This Matlab function set consists of some useful functions for phase retrival and image reconstruction that widely used in the community of X-ray coherent diffraction imaging (XCDI). 

# Function usage

## hio2d.m
Developed by Fienup [1], well known and widely used algorithm for XCDI resarch.Supports oversampling smoothness develped by Miao's group [2].

### Syntax
```
function R = hio2d(Fabs, S, n)
function R = hio2d(Fabs, S, n, ukwn, alpha)
```

## ghio.m
Developed by Chien-Chun Chen *et al*. when he was in Institute of Physics in Academia Sinica [3].

### Syntax
```
function R = ghio2d(Fabs, S, n, gen, rep, checker, alpha)
function [R, G, efs] = ghio2d(Fabs, S, n, gen, rep, checker, alpha)
```


## shrinkwrap.m
Developed by S. Marchesini *et al*. [4].

### Syntax
```
gshrinkwrap(Fabs, n1, uknwn, gen, n2);
gshrinkwrap(Fabs, n1, ukwn, gen, n2, alpha);
gshrinkwrap(Fabs, n1, ukwn, gen, n2, alpha, sigma, cutoff1, cutoff2);
```

### Description


## gshrinkwrap.m

### Syntax
```
gshrinkwrap(Fabs, n1, uknwn, gen, n2, rep);
gshrinkwrap(Fabs, n1, ukwn, gen, n2, rep, alpha);
gshrinkwrap(Fabs, n1, ukwn, gen, n2, rep, alpha, sigma, cutoff1, cutoff2);
```


# References

[1] J.R. Fienup, Appl. Opt. **21**, 2758 (1982).

[2] J.A. Rodriguez *et al*., J. Appl. Cryst. **46** (2013).

[3] C.-C. Chen *et al*., Phys. Rev. B **59**, 064113 (2007).

[4] S. Marchesini *et al*., Phys. Rev. B **68**, 140101 (2003).
