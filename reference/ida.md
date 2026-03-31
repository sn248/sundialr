# ida

IDA solver to solve stiff DAEs

## Usage

``` r
ida(
  time_vector,
  IC,
  IRes,
  input_function,
  Parameters,
  reltolerance = 1e-04,
  abstolerance = 1e-04
)
```

## Arguments

- time_vector:

  time vector

- IC:

  Initial Value of y

- IRes:

  Inital Value of ydot

- input_function:

  Right Hand Side function of DAEs

- Parameters:

  Parameters input to ODEs

- reltolerance:

  Relative Tolerance (a scalar, default value = 1e-04)

- abstolerance:

  Absolute Tolerance (a scalar or vector with length equal to ydot,
  default = 1e-04)

## Value

A data frame. First column is the time-vector, the other columns are
values of y in order they are provided.
