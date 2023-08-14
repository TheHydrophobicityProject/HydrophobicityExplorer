# from multiprocessing import Pool, Process
from rich.progress import track
from concurrent.futures import ProcessPoolExecutor
import time, os, random
import multiprocessing

from rich import progress

class dummy:
    def __init__(self, n, x=0, y=0):
        self.n = n
        self.x = x
        self.y = y

    def __str__(self):
        return f"{self.n = }, {self.x = }, {self.y = }"
        

def f(progress, task_id):
    n = random.randint(1, 10)
    x = random.randint(1, 10)
    y = random.randint(1, 10)
    dumb = dummy(n, y=y)
    dumb.x = dumb.n**dumb.y
    long_running_fn(progress, task_id)
    return dumb


def long_running_fn(progress, task_id):
    len_of_task = random.randint(3, 20)  # take some random length of time
    for n in range(0, len_of_task):
        time.sleep(.1)  # sleep for a bit to simulate work
        progress[task_id] = {"progress": n + 1, "total": len_of_task}
    # print(progress)

if __name__ == '__main__':

    n_workers = os.cpu_count()

    with progress.Progress(
        "[progress.description]{task.description}",
        progress.BarColumn(),
        "[progress.percentage]{task.percentage:>3.0f}%",
        progress.TimeRemainingColumn(),
        progress.TimeElapsedColumn(),
        refresh_per_second=5,  # bit slower updates
    ) as progress:
        futures = []
        with multiprocessing.Manager() as manager:
            # this is the key - we share some state between our 
            # main process and our worker functions
            _progress = manager.dict()
            overall_progress_task = progress.add_task("[green]All jobs progress:")
            
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                for n in range(0, 20):  # iterate over the jobs we need to run
                    # set visible false so we don't have a lot of bars all at once:
                    task_id = progress.add_task(f"task {n}", visible=False)
                    futures.append(executor.submit(f, _progress, task_id))
                
                # monitor the progress:
                while (n_finished := sum([future.done() for future in futures])) < len(futures):
                    progress.update(
                        overall_progress_task, completed=n_finished, total=len(futures)
                    )
                    for task_id, update_data in _progress.items():
                        latest = update_data["progress"]
                        total = update_data["total"]
                        # update the progress bar for this task:
                        progress.update(
                            task_id,
                            completed=latest,
                            total=total,
                            visible=latest < total,
                        )
                # raise any errors:
                results = [f.result() for f in futures]
                for result in results:
                    print(result)

    # for i in range(5):
    #     dummies.append(dummy(i))

    # const1 = 3
    # const2 = 4
    # const3 = 5

    # inps=[(dumb, const1) for dumb in dummies]

    # start = time.time()
    # with Pool(5) as p:
    #     results = p.starmap(f, inps)
    
    # for dumb in results:
    #     print(dumb.x)
    # end = time.time()
    # print('Parallel Execution Time: {}'.format(end-start))

