import csv
number_of_data = 12
n=3
seed=123

data_file_names = "n=4#1234-temp.dat"
write_file_names = "n=4#1234-temp.csv"

csvfile = open(write_file_names,'w',newline='')
datafile = open(data_file_names,'r',errors='ignore')
temp = datafile.readline()

result_writer = csv.writer(csvfile)
result_writer.writerow(["n="+"{}".format(n),"seed="+"{}".format(seed)])
result_writer.writerow(["beta","E","dE","E^2","dE^2","A1","dA1","A2","dA2","A3","dA3","A1^2","dA1^2","A2^2","dA2^2","A3^2","dA3^2","I01","I02","I03","I12","I13","I23"])


for i in range(number_of_data):
    temp = datafile.readline().split(' ')
    temp = temp[0].split('=')
    beta = temp[1]

    temp = datafile.readline().split('\t')
    E = temp[1]
    dE = temp[2]

    temp = datafile.readline().split('\t')
    E2 = temp[1]
    dE2 = temp[2]

    temp = datafile.readline().split('\t')
    A1 = temp[1]
    dA1 = temp[2]
    temp = datafile.readline().split('\t')
    A2 = temp[1]
    dA2 = temp[2]
    temp = datafile.readline().split('\t')
    A3 = temp[1]
    dA3 = temp[2]

    temp = datafile.readline().split('\t')
    sA1 = temp[1]
    dsA1 = temp[2]
    temp = datafile.readline().split('\t')
    sA2 = temp[1]
    dsA2 = temp[2]
    temp = datafile.readline().split('\t')
    sA3 = temp[1]
    dsA3 = temp[2]

    informaition = datafile.readline().split('\t')[0:6]
    result_writer.writerow([beta,E,dE,E2,dE2,A1,dA1,A2,dA2,A3,dA3,sA1,dsA1,sA2,dsA2,sA3,dsA3]+informaition)
    temp = datafile.readline()