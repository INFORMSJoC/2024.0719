import math
import numpy as np

K = 100
epsilon = 0.1
theta = math.ceil(math.log2(K))

# 贪婪算法
v = []
xxx = list(range(1, 21))
yyy = [ 1/math.cos( math.pi/(2**(item+1)) ) for item in xxx ] 
total_accuracy = 1
for layer in range(theta):
    for i in range(len(yyy)):
        if total_accuracy * yyy[i] < 1+epsilon:
            v.append(i+1)
            total_accuracy = total_accuracy * yyy[i]
            break

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

print(number_contraints)
 

 