import os
import sys
import time
import pickle
import random
import multiprocessing as multip

instance_file_name = sys.argv[1]
local_optimum_dict = {}
number_of_tests = 10

def run_two_opt():
    seed = random.randint(0, 10000)
    query = './bin/best_improv -f ../qap_instances/{} -s {}'.format(
        instance_file_name, str(seed))
    two_opt_result = os.popen(query).read()
    two_opt_result = two_opt_result.split('\n')
    solution = tuple(map(int, two_opt_result[0].split()))
    cost = int(two_opt_result[1])
    return solution, cost

def store_results(results):
    global solution_list, cost_list
    solution, cost = results
    solution_list.append(solution)
    cost_list.append(cost)

start_time = time.time()
elapsed_time = 0
while elapsed_time < 300:
    solution_list, cost_list = [], []
    pool = multip.Pool(processes=number_of_tests)
    for _ in range(number_of_tests):
        pool.apply_async(run_two_opt, callback=store_results)
    pool.close()
    pool.join()
    for solution, cost in zip(solution_list, cost_list):
        local_optimum_dict[solution] = cost
        if len(local_optimum_dict) == 5000:
            break
    else:
        elapsed_time = time.time() - start_time
        continue
    break

output_filename = 'QAP_best_improv_local_search' + instance_file_name
with open(output_filename, 'wb') as output_file:
    pickle.dump((local_optimum_dict, elapsed_time), output_file)
