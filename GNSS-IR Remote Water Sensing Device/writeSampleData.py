# Write some sample data for later use
import micropython
# Machine is a standard micropython library
from machine import Pin
import utime as time

# Asyncio I/O event scheduler
import uasyncio as asyncio

async def writeSampleData():
    # Open the IRDP file to place some data in there
    with open("/sd/irdp.csv", "w") as file:
        # write each list into the file
        file.write("1609459269,6,180,23\r\n")
        file.write("1609459270,49,180,43\r\n")
        file.write("1609459271,63,180,38\r\n")
        file.write("1609459272,12,180,41\r\n")
        file.write("1609459273,54,180,38\r\n")
        file.write("1609459274,44,180,36\r\n")
        file.write("1609459275,25,180,43\r\n")
        file.write("1609459276,11,180,32\r\n")
        file.write("1609459277,39,180,43\r\n")
        file.write("1609459278,4,180,37\r\n")
    
    # Read to see if it worked
    with open("/sd/irdp.csv", "r") as file:
        data = file.read()
        print(data)
        
    # wait 1s after done
    await asyncio.sleep(1)
        