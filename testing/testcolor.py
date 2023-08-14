from tqdm import trange, tqdm
from rich.progress import track
import time

# # Define your function or loop
# def some_function():
#     for i in trange(100, desc="Processing", ncols=80, colour='MAGENTA'):
#         # Simulate some processing time
#         time.sleep(0.1)

# # Call your function
# some_function()

lst = range(20)

for i in track(lst, description="Processing...", disable=True):
    print(lst[i])
    time.sleep(.1)  # Simulate work being done


# from rich.progress import Progress

# with Progress() as progress:

#     task1 = progress.add_task("[red]Downloading...", total=1000)
#     task2 = progress.add_task("[green]Processing...", total=1000)
#     task3 = progress.add_task("[cyan]Cooking...", total=1000)

#     while not progress.finished:
#         progress.update(task1, advance=0.5)
#         progress.update(task2, advance=0.3)
#         progress.update(task3, advance=0.9)
#         time.sleep(0.02)