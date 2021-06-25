li = []

def f(li, i):
    li.append(i)

N = 10000000

import time
tic = time.perf_counter()
li = []
for i in range(N):
    if i % 2:
        li.append(1)
    elif i % 3:
        li.append(1)
toc = time.perf_counter()
print(f'one loop took {toc - tic:0.4f} seconds')

import time
tic = time.perf_counter()
li = []
for i in range(N):
    if i % 2:
        li.append(1)
for i in range(N):
    if i % 3 and not i % 2:
        li.append(1)
toc = time.perf_counter()
print(f'two loops took {toc - tic:0.4f} seconds')
import time
tic = time.perf_counter()
li = [1 for i in range(N) if i % 2 and i % 3]
toc = time.perf_counter()
print(f' list comp took {toc - tic:0.4f} seconds')
# # import time
# # tic = time.perf_counter()

# # [f(li, i) for i in range(N)]

# # toc = time.perf_counter()
# # print(f'list comp took {toc - tic:0.4f} seconds')

# # import time
# # tic = time.perf_counter()

# # li = [0] * N
# # for i in range(N):
# #     li[i] = i


# # toc = time.perf_counter()
# # print(f'for loop took {toc - tic:0.4f} seconds')

# setN = set(range(N))
# import time
# tic = time.perf_counter()

# li = set()
# for i in setN:
#     li.add(i)


# toc = time.perf_counter()
# print(f'add set took {toc - tic:0.4f} seconds')

# listN = list(range(N))
# import time
# tic = time.perf_counter()

# li = []
# for i in listN:
#     li.append(i)


# toc = time.perf_counter()
# print(f'append list took {toc - tic:0.4f} seconds')

# import time
# tic = time.perf_counter()

# li = [0] * N
# for i in listN:
#     li[i] = i

# toc = time.perf_counter()
# print(f'pre list took {toc - tic:0.4f} seconds')
