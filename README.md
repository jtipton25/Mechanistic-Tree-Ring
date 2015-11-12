# Mechanistic-Tree-Ring

This is a repository for code associated with the paper: "Reconstruction of late Holocene climate based on tree growth and mechanistic hierarchical models."


### Installation of `R` packages 

Users must install the two tarball packages myFunctions_1.0.tar.gz and bayesTreeRing_1.0.tar.gz to load the necessary c++ helper functions for code to run. Once downloaded, move to the download directory and use the commands

```
R CMD INSTALL myFunctions_1.0.tar.gz
```
and
```
R CMD INSTALL bayesTreeRing_1.0.tar.gz
```
To install the needed R packages.

### Climate Data
Climate data is accessed from the `.RData` file using the command
```
data("hudsonValleyData")
```

The temperature data used for the Hudson Valley reconstruction are `Temp.avg.dat` and the Precipitation data are `Precip.dat`.

