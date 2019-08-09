import csv
number_of_data = 7
n=3
for seed in ["1234567"]:
    data_file_names = "n={}#{}.dat".format(n,seed)
    write_file_names = "n={}#{}.csv".format(n,seed)

    csvfile = open(write_file_names,'w',newline='')
    datafile = open(data_file_names,'r',errors='ignore')
    # temp = datafile.readline()

    result_writer = csv.writer(csvfile)
    result_writer.writerow(["n="+"{}".format(n),"seed="+"{}".format(seed)])
    result_writer.writerow(["beta","E","dE","AutoE","E^2","dE^2","AutoE2","A1","dA1","AutoA1","A2","dA2","AutoA2","A3","dA3","AutoA3","A1^2","dA1^2","AutoA1^2","A2^2","dA2^2","AutoA2^2","A3^2","dA3^2","AutoA3^2","I01","I02","I03","I12","I13","I23"])


    for i in range(number_of_data):
        temp = datafile.readline().split(' ')
        temp = temp[0].split('=')
        beta = temp[1]

        temp = datafile.readline().split('\t')
        E = temp[1]
        dE = temp[2]
        AutoE = float(temp[3])

        temp = datafile.readline().split('\t')
        E2 = temp[1]
        dE2 = temp[2]
        AutoE2 = float(temp[3])

        temp = datafile.readline().split('\t')
        A1 = temp[1]
        dA1 = temp[2]
        AutoA1 = float(temp[3])
        temp = datafile.readline().split('\t')
        A2 = temp[1]
        dA2 = temp[2]
        AutoA2 = float(temp[3])
        temp = datafile.readline().split('\t')
        A3 = temp[1]
        dA3 = temp[2]
        AutoA3 = float(temp[3])
        temp = datafile.readline().split('\t')
        sA1 = temp[1]
        dsA1 = temp[2]
        AutosA1 = float(temp[3])
        temp = datafile.readline().split('\t')
        sA2 = temp[1]
        dsA2 = temp[2]
        AutosA2 = float(temp[3])
        temp = datafile.readline().split('\t')
        sA3 = temp[1]
        dsA3 = temp[2]
        AutosA3 = float(temp[3])
        informaition = datafile.readline().split('\t')[0:6]
        result_writer.writerow([beta,E,dE,AutoE,E2,dE2,AutoE2,A1,dA1,AutoA1,A2,dA2,AutoA2,A3,dA3,AutoA3,sA1,dsA1,AutosA1,sA2,dsA2,AutosA2,sA3,dsA3,AutosA3]+informaition)
        temp = datafile.readline()