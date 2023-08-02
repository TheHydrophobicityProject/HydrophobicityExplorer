from tqdm import trange, tqdm
import time

# Define your function or loop
def some_function():
    for i in trange(100, desc="Processing", ncols=80, colour='MAGENTA'):
        # Simulate some processing time
        time.sleep(0.1)

# Call your function
some_function()

