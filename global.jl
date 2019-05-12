Random.seed!(12345)
const SpinValueType = Int8
const layers = 3
const dim = 3
const Nsamp = 100000
const Nblck = 1024
const CUT_ENERGY = 100  # the value in WangLandau for cut off low degeneracy state
 #----Global Variables

 #number of A,B,Vectors, and normailzation coefficient
const NA = layers*dim^2 
const NB = layers*dim
const NV = 2^dim
const NC = 1/(NV*dim)
const Jcp = collect(Float64,0:1:1)

#--------------The set up for version 1.0
# Random.seed!(12345)
# const SpinValueType = Int8
# const layers = 3
# const dim = 3
# const Nsamp = 1000
# const Nblck = 1024
# const CUT_ENERGY = 100  # the value in WangLandau for cut off low degeneracy state
#  #----Global Variables

#  #number of A,B,Vectors, and normailzation coefficient
# const NA = layers*dim^2 
# const NB = layers*dim
# const NV = 2^dim
# const NC = 1/(NV*dim)
# const Jcp = collect(Float64,0:1:1)
