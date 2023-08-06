from multiprocessing import Pool, Process
import time

class dummy:
    def __init__(self, n, x=0, y=0):
        self.n = n
        self.x = x
        self.y = y

def f(dummy,y = 0):
    time.sleep(1)
    dummy.x = dummy.n**y
    return dummy

if __name__ == '__main__':

    dummies = []
    for i in range(5):
        dummies.append(dummy(i))

    const1 = 3
    const2 = 4
    const3 = 5

    inps=[(dumb, const1) for dumb in dummies]

    start = time.time()
    with Pool(5) as p:
        results = p.starmap(f, inps)
    
    for dumb in results:
        print(dumb.x)
    end = time.time()
    print('Parallel Execution Time: {}'.format(end-start))

    # start = time.time()
    # results = []
    # for n in [1,2,3]:
    #     results.append(f(n))
    # print(results)
    # end = time.time()
    # print('Sequential Execution Time: {}'.format(end-start))

    
# const1 = 3
# const2 = 4
# const3 = 5
# lst = ["1", "2", "3", "4", "5", "6"]

# iterargs = [(flat, const1, const2, const3) for flat in lst]

# print(iterargs)