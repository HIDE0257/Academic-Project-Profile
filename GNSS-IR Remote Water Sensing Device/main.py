# Import micropython libraries
from micropython import alloc_emergency_exception_buf
from gc import collect

# Asyncio I/O event scheduler
from uasyncio import run

# import the main class that will handle everything
from mainClass import mainClass

# Create an emergency exceiption buffer to store debug messages
alloc_emergency_exception_buf(100)

# Define our main loop, which can run asynchronous code with asyncio
async def main():
    # Construct the main loop
    main = mainClass()
    # Create a task for the main loop
    await main.mainLoop()
    
# Run the main loop
while True:
    run(main())
    # run garbage collection
    collect()