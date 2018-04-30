Bootstrapping and hazard ratio thresholds
================
Dominic Pearce, The Institute of Genetics and Molecular Medicine, The University of Edinburgh
2018-04-30



 

Bootstrapping and hazard ratio thresholds
-----------------------------------------

 

### Libraries

``` r
library(survivALL)
library(Biobase)
library(knitr)
```

 

To determine and ensure reliable prognostic association as a measure of significance, *survivALL* can perform a non-parametric bootstrapping procedure. In short we calculate, for each point-of-separation a distribution of expected hazard ratios (HRs), against which we're able to compare our observed HRs as part of our analysis.

To achieve this, we randomly sample our survival data with replacement and then calculate survival statistics for all points-of-separation, exactly as we would for a biomarker under investigation. By repeating this procedure 1,000s or 10,000s of times, we produce our distribution of *expected* hazard ratios.

 

``` r
data(nki_subset)

#bootstrapping data should be in the format of 1 repeat per column
bs_mtx <- matrix(nrow = ncol(nki_subset), ncol = 20)

system.time(
            for(i in 1:ncol(bs_mtx)){
                bs_mtx[, i] <- allHR(measure = sample(1:ncol(nki_subset), 
                                                      replace = TRUE),
                                     srv = pData(nki_subset),
                                     time = "t.dmfs",
                                     event = "e.dmfs")
            }
)
```

user system elapsed 25.547 0.124 25.995

``` r

kable(bs_mtx[1:20, 1:5])
```

|            |            |            |           |            |
|-----------:|-----------:|-----------:|----------:|-----------:|
|  -1.4869615|          NA|          NA|         NA|          NA|
|  -1.3504229|          NA|          NA|         NA|  -1.7082819|
|  -1.5213504|   0.1817535|  -0.1495991|         NA|          NA|
|  -1.0592437|   0.5781454|  -1.0484108|         NA|          NA|
|  -1.2792451|  -0.3664483|  -0.4756618|         NA|  -1.4188765|
|  -0.9183075|  -0.7742038|  -1.0205127|         NA|  -1.0457754|
|  -0.5600026|  -0.3516894|  -0.6952990|         NA|  -0.7446602|
|  -0.6625810|  -0.0846293|  -0.4491413|         NA|  -0.4585176|
|  -0.4402430|  -0.3920619|  -0.1236995|         NA|  -0.7981985|
|  -0.6104234|  -0.5738814|   0.1220281|         NA|  -0.5706505|
|  -0.8201632|  -0.7471677|  -0.2475666|         NA|  -0.3688480|
|  -0.6410602|  -0.5235132|  -0.0946439|         NA|  -0.6812726|
|  -0.5232169|  -0.6256363|   0.1241489|         NA|  -0.5315006|
|  -0.4092515|  -0.8223939|   0.2379108|         NA|  -0.3805137|
|  -0.3006547|  -0.7062154|  -0.0470965|         NA|  -0.2579031|
|  -0.2065570|  -0.6527314|  -0.2364833|         NA|  -0.1301706|
|  -0.0653427|  -0.5516527|  -0.0712919|         NA|  -0.0248992|
|  -0.2597032|  -0.3988798|   0.0287874|         NA|   0.1212919|
|  -0.4079458|  -0.4995051|   0.1177846|  0.9829079|   0.2787528|
|  -0.4867013|  -0.6063988|   0.1966554|  1.0553475|   0.4116118|

 

Having calculated our bootstrapped data we then simply hand the matrix to either the `survivALL()` or `plotALL()` functions (using the `bs_dfr =` argument) to handle the subsequent significance calculations. It should be noted that bootstrapping up to 10,000x can be a long process requiring an investment of time.
