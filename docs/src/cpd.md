# CPD model

The acronym "CPD" stands for "Chemical Percolation and Devolatilization".  To [cite the original CPD authors](http://www.et.byu.edu/~tom/cpd/cpd.html):

> The chemical percolation devolatilization (CPD) model was developed to model
> coal devolatilization based on characteristics of the chemical structure of
> the parent coal. [It] describes the devolatilization behavior of rapidly
> heated coal based on the chemical structure of the parent coal.

The model predicts gas, tar, metaplast and char yield of coal pyrolysis over the course of the process.  It is based on percolation statistics theory, where the coal is modeled as a matrix of connected aromatic 'clusters' or 'sites' connected by bridges that may break up, react and re-link during pyrolysis, thereby creating finite fragments that will end up as released volatiles, metaplast, or re-linked back into the original material.

In the literature, the model has also been used to model biomass devolatilization, by considering biomass a sum of its main components, namely cellulose, hemicellulose and lignin. See e.g. [Fletcher et al. 2012](https://pubs.acs.org/doi/abs/10.1021/ef300574n) or [Lewis and Fletcher, 2013](https://pubs.acs.org/doi/10.1021/ef3018783).

## Implementation

The `cpd()` function below aims to provide a basic Julia implementation of the CPD model as described in the following early papers:

- [Chemical model of coal devolatilization using percolation lattice statistics](https://pubs.acs.org/doi/pdf/10.1021/ef00014a011)
  David M. Grant, Ronald J. Pugmire, Thomas H. Fletcher, and Alan R. Kerstein (1988)
- [Chemical percolation model for devolatilization. 2. Temperature and heating rate effects on product yields](https://pubs.acs.org/doi/pdf/10.1021/ef00019a010)
  Thomas H. Fletcher, Alan R. Kerstein, Ronald J. Pugmire, and David M. Grant (1990)
- [Chemical percolation model for devolatilization. 3. Direct use of carbon-13 NMR data to predict effects of coal type](https://pubs.acs.org/doi/10.1021/ef00034a011)
  Thomas H. Fletcher, Alan R. Kerstein, Ronald J. Pugmire, Mark S. Solum, and David M. Grant (1992)

In addition, the `cpd()` function proposes an experimental, modified metaplast model to ensure mass conservation when integrating with the spatial model.

## Online resources

The original creators of the CPD model have made a number of resources [available online](http://www.et.byu.edu/~tom/devolatilization/CPD%20model.html).  These web pages are provided and maintained by Thomas Fletcher, who co-authored the original papers cited above.   The provided material also includes Fortran and Matlab code for the original CPD model.

We have significantly benefitted from these resources during the development of `CPDSpatial`, but do not claim any affiliation with the  authors of the original papers.

# Source code documentation
```@docs
ReactionRateParams
MaterialParams
cpd(AEσb::ReactionRateParams, AEσg::ReactionRateParams,
AEσρ::ReactionRateParams, mpar::MaterialParams,
duration::Float64,Tfun::Function; metaplast_model, Pfun::Function, num_tar_bins, max_tstep)

```
