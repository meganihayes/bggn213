Class07
================
Megan Hayes
January 30, 2019

Functions revisted
------------------

Load (i.e. **source**) our rescale() function from last class.

``` r
source("http://tinyurl.com/rescale-R")
```

    ## Warning in file(filename, "r", encoding = encoding): "internal" method
    ## cannot handle https redirection to: 'https://bioboot.github.io/bggn213_f17/
    ## class-material/rescale.R'

    ## Warning in file(filename, "r", encoding = encoding): "internal" method
    ## failed, so trying "libcurl"

Test this function

``` r
rescale(1:5)
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

``` r
#rescale(c(1:5, "string"))
```

We want to make this function more robust

``` r
#rescale2(c(1:5, "string"))
```

``` r
is.numeric(1:5)
```

    ## [1] TRUE

``` r
is.numeric("string")
```

    ## [1] FALSE

``` r
is.numeric(c(1:5, "string"))
```

    ## [1] FALSE

``` r
!is.numeric(c(1:5, "string"))
```

    ## [1] TRUE

``` r
!is.numeric(1:5)
```

    ## [1] FALSE

``` r
x <- c( 1, 2, NA, 3, NA)
y<-c(NA,3,NA,3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum( is.na(x) & is.na(y) )
```

    ## [1] 1

``` r
both.na <- function(x, y) {
  sum( is.na(x) & is.na(y) )
}
```

``` r
both.na(x,y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c(1, NA, NA)
y2 <- c(1, NA, NA, NA)

both.na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
both_na3(x,y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3
