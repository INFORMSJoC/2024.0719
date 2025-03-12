import matplotlib.pyplot as plt
import math

def print_info(K, epsilon, theta, v):
    number_contraints = 0
    number_variables = 2*(2**theta) - 1
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            number_contraints += 4*v[l-1] + 6
            number_variables += 2*v[l-1] + 2
    print(f'Polyhedral Approximation (#var and #cons): ({number_variables}, {number_contraints})')

    number_variables = 2**theta # layer 0: y 
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            number_variables += 2
            for j in range(1, v[l-1]+1):
                number_variables += 1
    number_variables += 1 # t

    number_contraints = 0
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            if l == 1:
                number_contraints += 4
            else:
                number_contraints += 2

    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            number_contraints += 1
            for j in range(1, v[l-1]+1):
                number_contraints += 2
    number_contraints += 1
    print(f'Simplified Polyhedral Approximation (#var and #cons): ({number_variables}, {number_contraints})')

    layers = 0
    layers_count = [] # 每一 成的数目 
    current_layer_count = 2**theta
    while current_layer_count != 1:
        layers += 1
        current_layer_count = current_layer_count // 2 + current_layer_count%2
        layers_count.append(current_layer_count)

    number_contraints = 0
    # 第一层
    for i in range(layers_count[0]):
        # j 就是index_mul  v是multip_number
        for index_mul in range(0, 2**(v[0]-1)):
            number_contraints += 1


    for lay in range(1, layers):
        for i in range(layers_count[lay]):
            # j 就是index_mul  v是multip_number
            for index_mul in range(0, 2**(v[lay-1]-1)):
                number_contraints += 1
    number_variables = 2*2**theta - 1
    print(f'Our Polyhedral Approximation (#var and #cons): ({number_variables}, {number_contraints})')
 
# print(y)
# plt.plot( x, y, 's-', color = 'r', label="ATT-RLSTM")  
# plt.xlabel("v")  
# plt.ylabel("accuracy") 
# plt.legend(loc = "best") 
# plt.show() 2 4 6 

# v = [2,3]
# layer_accuracy = [1/math.cos( math.pi/(2**(item+1)) ) for item in v]
# total_accuracy = 1
# for i in layer_accuracy:
#     total_accuracy = total_accuracy * i 
# print(total_accuracy)
 

K = 200
epsilon = 0.001
theta = math.ceil(math.log2(K))
print("*******************************************")
################################################## 0. 策略
v = []
xxx = list(range(1, 21))
yyy = [ 1/math.cos( math.pi/(2**(item+1)) ) for item in xxx ] 
kkk = [1, 0.618, 0.544, 0.519, 0.509, 0.504, 0.502, 0.501, 0.5, 0.5, 0.5]
total_accuracy = 1
for layer in range(theta):
    for i in range(len(yyy)):
        if yyy[i] <= 1 + kkk[theta-1]**(1+layer) * epsilon:
            v.append(i+1)
            total_accuracy = total_accuracy * yyy[i]
            break

print(f'total_accuracy: {total_accuracy}, and v: {v}.') 
print_info(K, epsilon, theta, v)

# print("===========================================")
# ################################################## 1. 贪婪策略
# v = []
# xxx = list(range(1, 21))
# yyy = [ 1/math.cos( math.pi/(2**(item+1)) ) for item in xxx ] 
# total_accuracy = 1
# for layer in range(theta):
#     for i in range(len(yyy)):
#         if total_accuracy * yyy[i] < 1+epsilon:
#             v.append(i+1)
#             total_accuracy = total_accuracy * yyy[i]
#             break

# print(f'total_accuracy: {total_accuracy}, and v: {v}.') 
# print_info(K, epsilon, theta, v)

    


# print("===========================================")
# ################################################## 2. same v strategy
# # v = [math.ceil( math.log2( math.pi/math.acos( 1/ ((1+epsilon)**(1/theta)) ) ) )  for i in range(1, theta+1)]
# for v_value in range(1, 15):
#     layer_accuracy = [1/math.cos( math.pi/(2**(v_value+1)) ) for i in range(1, theta+1)]
#     total_accuracy = 1
#     for i in layer_accuracy:
#         total_accuracy = total_accuracy * i 
#     if total_accuracy <= 1+epsilon:
#         break
# v = [v_value  for i in range(1, theta+1)]
# print(f'total_accuracy: {total_accuracy}, and v: {v}.') 
# print_info(K, epsilon, theta, v)

# print("===========================================")
# ################################################## 3. benta strategy

# for o1 in range(1, 201):
#     v = [math.floor((o1/100)*i*math.log(2/epsilon))  for i in range(1, theta+1)]
#     layer_accuracy = [1/math.cos( math.pi/(2**(item+1)) ) for item in v]
#     total_accuracy = 1
#     for i in layer_accuracy:
#         total_accuracy = total_accuracy * i 
#     if total_accuracy <= 1+epsilon:
#         break
# v = [math.floor((o1/100)*i*math.log(2/epsilon))  for i in range(1, theta+1)]
# print(f'total_accuracy: {total_accuracy}, and v: {v}.') 
# # print_info(K, epsilon, theta, v)


# number_contraints = 0
# number_variables = 2*(2**theta) - 1
# for l in range(1, theta+1):
#     for i in range(1, 2**(theta-l)+1):
#         number_contraints += 4*v[l-1] + 6
#         number_variables += 2*v[l-1] + 2
# print(f'Polyhedral Approximation (#var and #cons): ({number_variables}, {number_contraints})')

# number_variables = 2**theta # layer 0: y 
# for l in range(1, theta+1):
#     for i in range(1, 2**(theta-l)+1):
#         number_variables += 2
#         for j in range(1, v[l-1]+1):
#             number_variables += 1
# number_variables += 1 # t

# number_contraints = 0
# for l in range(1, theta+1):
#     for i in range(1, 2**(theta-l)+1):
#         if l == 1:
#             number_contraints += 4
#         else:
#             number_contraints += 2

# for l in range(1, theta+1):
#     for i in range(1, 2**(theta-l)+1):
#         number_contraints += 1
#         for j in range(1, v[l-1]+1):
#             number_contraints += 2
# number_contraints += 1
# print(f'Simplified Polyhedral Approximation (#var and #cons): ({number_variables}, {number_contraints})')
 








 