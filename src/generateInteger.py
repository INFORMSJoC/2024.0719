from scipy.stats import ortho_group
import random
import json
import sys
import numpy as np
from util import decompose_sysmtic_matrix


def generate():
    N = int(sys.argv[1])
    M = 2
    # positive_number = int(sys.argv[2])

    data_bound = 20
    for data_index in range(20):
        Q = [[0 for j in range(N)] for i in range(N)]
        for i in range(N):
            for j in range(i, N):
                Q[i][j] = random.randrange(0,data_bound)
        for i in range(N):
            for j in range(i):
                Q[i][j] = Q[j][i]
        # aaa,bbb = decompose_sysmtic_matrix(N, Q)
        # for row in range(len(aaa)):
        #     for col in range(len(aaa[0])):
        #         aaa[row][col] = round(aaa[row][col],2)
        # for row in range(len(bbb)):
        #     for col in range(len(bbb[0])):
        #         bbb[row][col] = round(bbb[row][col],2)
        # # print(aaa)
        # # print(bbb)
        # print(np.dot(aaa,np.transpose(aaa)) - np.dot(bbb,np.transpose(bbb)))
        # print(Q)

        c = [0 for i in range(N)]
        for i in range(N):
            c[i] = random.randrange(-data_bound,data_bound)

        b = [0 for i in range(M)]
        for i in range(M):
            b[i] = random.randrange(-data_bound,data_bound) + 50

        A = []
        C = []

        posi_evalue = []
        for m in range(M):
            atemp = [[0 for j in range(N)] for i in range(N)]
            for i in range(N):
                for j in range(i, N):
                    atemp[i][j] = random.randrange(-data_bound,data_bound)
            for i in range(N):
                for j in range(i):
                    atemp[i][j] = atemp[j][i]

            aaa, _ = decompose_sysmtic_matrix(N, atemp)
            posi_evalue.append(len(aaa[0]))

            # for row in range(len(aaa)):
            #     for col in range(len(aaa[0])):
            #         aaa[row][col] = round(aaa[row][col],2)
            # for row in range(len(bbb)):
            #     for col in range(len(bbb[0])):
            #         bbb[row][col] = round(bbb[row][col],2)
            # # print(aaa)
            # # print(bbb)
            # print(np.dot(aaa,np.transpose(aaa)) - np.dot(bbb,np.transpose(bbb)))
            # print(atemp)
            A.append(atemp)

            ctemp = [0 for i in range(N)]
            for i in range(N):
                ctemp[i] = random.randrange(-data_bound,data_bound)
            C.append(ctemp)


        data_dict = {
            "N" : N,
            "M" : M,
            "Q" : Q,
            "c" : c,
            "A" : A,
            "Y" : [],
            "Z" : [],
            "C" : C,
            "b" : b,
            "u" : posi_evalue
        }
        with open("./QCQPData2/" + str(N) + "_" +str(data_index) + ".json","w") as f:
            json.dump(data_dict,f)

if __name__ == "__main__":
    generate()