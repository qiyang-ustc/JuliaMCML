## Julia MCML

### Path Setting: 

The second line in "mcml.jl": 

modify  `abspath("D:\\IOPHY\\jmc\\")` to your path 

### Global Settings: 

in "global.jl"

```julia
Random.seed!(12345)
const SpinValueType = Int8
const layers = 3
const dim = 3 
const Nsamp = 1000
const Nblck = 1024
const CUT_ENERGY = 100  # the value in WangLandau for cut off low degeneracy state
 #----Global Variables

 #number of A,B,Vectors, and normailzation coefficient
 # Auto----
const NA = layers*dim^2 
const NB = layers*dim
const NV = 2^dim
const NC = 1/(NV*dim)
#----
const Jcp = collect(Float64,0:1:1) 
#If you use Wang Landau. You will not need this line

```

### Target Setting:

Set **target** by modifying

```julia
A[:,1] = [1 -1 -1;-1 1 -1;1 1 1]
A[:,2] = [1 1 -1 ;-1 1  1;1 -1 1]
A[:,3] = [1 -1 -1; 1 1 -1;1 1 -1]
B[:,1] = [-1;-1;1]
B[:,2] = [1 ; 1;1]
B[:,3] = [1 ;-1;-1]
C[:,1] = [-0.5;0.5,-0.5]
C[:,2] = [0.5;0.5;0.5]
C[:,3] = [-0.5;0.5;-0.5]
```



