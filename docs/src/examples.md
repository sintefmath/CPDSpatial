# Examples

These examples are intended to demonstrate the functionality and use of `CPDSpatial`.

## Testing the CPD model without a spatial model

This example is provide to demonstrate the use of the basic `cpd()` algorithm
provided by `CPDSpatial`.  This algorithm aims to reproduce the results of the
[original FORTRAN/Matlab computer code](http://www.et.byu.edu/~tom/cpd/cpdcodes.html)
issued by the original inventors/researchers.  Not all the features have been
included, such as nitrogen release calculations and light gas composition.  On
the other hand, a modified, experimental metaplast model has been proposed in
addition to the original metaplast model.

```@docs
    cpd_benchmarking
```

## Testing the CPD model embedded in a spatial model using JutulDarcy

### Simple case
### Composite biomaterial case
