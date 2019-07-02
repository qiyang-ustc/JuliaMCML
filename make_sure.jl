#this is result for model - A (energy spectrum)
#----FOR MODEL A-----
t = [0 133488;0.5 12075120;1 14436612;1.5 297648738; 2 1140637138;2.5 7371173628;3 5334431624; 3.5 9161971938; 4 5908502738;  4.5 9980021530;  5 5099541576;  5.5 7212254910;  6 3815706232;  6.5 5599512690;  7 2359841892;  7.5 2622590122;  8 1097687126;  8.5 1077445458;  9 280104776;  9.5 248642436;  10 61531918;  10.5 21086502; 11 1412928;11.5 1076832;12 8784]
#----FOR MODEL A-----
tsum = 68719476736
energy = t[:,1]/12
count = t[:,2]/tsum
beta = 5.0
# print("Ps=",sum(count),'\n')
# print("E =",sum(energy.*exp.(-beta.*energy).*count)/sum(exp.(-beta.*energy).*count),'\n')
# print("E2=",sum(energy.*energy.*exp.(-beta.*energy).*count)/sum(exp.(-beta.*energy).*count))

#----FOR MODEL B-----
#count = [0.00000918,0.000010645,0.000153906,0.000249316,0.001397949,0.003188965,0.012978711,0.013878125,0.052058008,0.043127148,0.07578916,0.051986035,0.106209277,0.071886426,0.1358375,0.103907813,0.187310352,0.060552539,0.053367969,0.015180762,0.00825166,0.001361133,0.001168945,0.000078516,0.000059961]
print("Ps=",sum(count),'\n')
print("B =",beta,'\n')
print("E =",sum(energy.*exp.(-beta.*energy).*count)/sum(exp.(-beta.*energy).*count),'\n')
print("E2=",sum(energy.*energy.*exp.(-beta.*energy).*count)/sum(exp.(-beta.*energy).*count))
#----FOR MODEL B-----
