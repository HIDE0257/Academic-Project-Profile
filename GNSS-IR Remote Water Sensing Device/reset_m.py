from os import remove
from machine import reset

# Function for persistent resets
def persistent(prn_list, msgQueue):
    # Delete the traceback file
    try:
        remove("/sd/traceback.txt")
    except:
        pass
    
    # Save the current PRN queue to a file to make it persistent
    try:
        with open("/sd/prn.csv", "w") as file:
            # Separate prns and write to file
            for prn in prn_list:
                file.write(str(prn[0]) + ',' + str(prn[1]) + '\n')
    except:
        pass

    # Save current message queue to a file to make it persistent
    try:
        with open("/sd/msg.txt", "w") as file:
            # Separate prns and write to file
            file.write(msgQueue)
    except:
        pass

    # Reset
    reset()

def nonper():
    # call a reset
    reset()