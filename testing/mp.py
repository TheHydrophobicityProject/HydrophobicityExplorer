from multiprocessing import Pool, Process
import time

def f(x):
    time.sleep(1)
    return x*x

if __name__ == '__main__':
    start = time.time()
    with Pool(5) as p:
        print(p.map(f, [1, 2, 3]))
    end = time.time()
    print('Parallel Execution Time: {}'.format(end-start))

    start = time.time()
    results = []
    for n in [1,2,3]:
        results.append(f(n))
    print(results)
    end = time.time()
    print('Sequential Execution Time: {}'.format(end-start))