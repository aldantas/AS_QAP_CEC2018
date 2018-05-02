import os
import sys
import pickle
import random
import time
import multiprocessing as multip

instance_file_name = sys.argv[1]
local_optimum_dict = {}
pseudo_global_optimum_dict = {}
best_cost = float('inf')
number_of_tests = 10

def run_ils():
    seed = random.randint(0, 10000)
    query = './bin/ils -f ../qap_instances/{} -s {}'.format(instance_file_name,
                                                        str(seed))
    ils_result = os.popen(query).read()
    ils_result = ils_result.split('\n')
    solution = tuple(map(int, ils_result[0].split()))
    cost = int(ils_result[1])
    return solution, cost

def store_results(results):
    global solution_list, cost_list
    solution, cost = results
    solution_list.append(solution)
    cost_list.append(cost)

start_time = time.time()
for _ in range(100):
    solution_list, cost_list = [], []
    pool = multip.Pool(processes=number_of_tests)
    for _ in range(number_of_tests):
        pool.apply_async(run_ils, callback=store_results)
    pool.close()
    pool.join()
    for solution, cost in zip(solution_list, cost_list):
        if cost < best_cost:
            best_cost = cost
            best_solution = solution
            pseudo_global_optimum_dict = {best_solution: best_cost}
        elif cost == best_cost and solution not in pseudo_global_optimum_dict:
            pseudo_global_optimum_dict[solution] = cost
        local_optimum_dict[solution] = cost
    elapsed_time = time.time() - start_time
    if elapsed_time >= 300:
        break
output_filename = 'QAP_iterated_local_search_' + instance_file_name
with open(output_filename, 'wb') as output_file:
    pickle.dump((local_optimum_dict, pseudo_global_optimum_dict, elapsed_time),
                output_file)
