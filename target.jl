#Don't modify the following 3 lines
A = ones(SpinValueType,dim,dim,layers)
B = ones(SpinValueType,dim,layers)
C = 0.5*ones(SpinValueType,dim,layers)
# you could modify A,B,C here just like this
#A[:,1] = [1 2 3;4 5 6;7 8 9]
#B[:,1] = [1;2;3]

A[:,:,1] = [1 -1 -1;-1 1 -1;-1 -1 1]
A[:,:,2] = [1 1 -1 ;-1 -1 -1;1 1 1]
A[:,:,3] = [1 -1 -1; -1 1 -1;-1 1 -1]
B[:,1] = [1 ;-1;1]
B[:,2] = [1 ; 1;-1]
B[:,3] = [1 ; 1;-1]
C[:,1] = [-0.5;-0.5;0.5]
C[:,2] = [0.5;-0.5;-0.5]
C[:,3] = [-0.5;-0.5;0.5]

# #---------Model C---------
# A[:,:,1] = [1 -1 -1;-1 1 -1;1 1 1]
# A[:,:,2] = [1 1 -1 ;-1 1  1;1 -1 1]
# A[:,:,3] = [1 -1 -1; 1 1 -1;1 1 -1]
# B[:,1] = [-1;-1;1]
# B[:,2] = [1 ; 1;1]
# B[:,3] = [1 ;-1;-1]
# C[:,1] = [-0.5;0.5;-0.5]
# C[:,2] = [0.5;0.5;0.5]
# C[:,3] = [-0.5;0.5;-0.5]
#------------------------------------
#Don't modify the following lines
A = SpinValueType.(A)
B = SpinValueType.(B)
C = Float64.(C)