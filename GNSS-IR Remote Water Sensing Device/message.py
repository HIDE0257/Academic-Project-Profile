from log import log_and_print
from gc import collect
from uasyncio import create_task
from sys import print_exception
from utime import time

# Message handling functions
def message_proc(message, reboot_event, updateSample, updateMask):

    try:
        # Enter a while loop based on the length of the message
        while len(message) > 0:
            # Read off the header of the message
            header = int(message[:2],16)
            # Cull the header from the message
            message = message[2:]

            if header == 0:
                # Reboot the system by setting the reboot event
                log_and_print('Reboot command recieved, spawning reboot event')
                # Create the file so we know this was a reboot
                with open("/sd/reboot", "w") as file:
                    file.write("Reboot request recieved at time " + str(time()))
                reboot_event.set()

            elif header == 1:
                log_and_print("New Masks recieved, updating")
                # Update the mask
                aziStart = int(message[0:4], 16)
                aziEnd = int(message[4:8], 16)
                eleStart = int(message[8:10], 16)
                eleEnd = int(message[10:12], 16)
                minHeight = int(message[12:14], 16)
                maxHeight = int(message[14:16], 16)
                # Remove the footer from the message
                message = message[16:]

                # Generate a task to update the maks
                updateMask((aziStart, aziEnd, eleStart, eleEnd, minHeight, maxHeight))

            elif header == 2:
                log_and_print("Downsampling update recieved, updating")
                # Update the sampling rate
                sampleRate = int(message[0:2], 16)
                closeoutTime = int(message[2:6], 16)
                # Remove the footer fromthe message
                message = message[6:]

                # Generate the task to update the sampling information
                create_task(updateSample(closeoutTime, sampleRate))

            else:
                # Other headers are erroneous, so print to console
                log_and_print("Unknown MT Message header, will do nothing")
        del message
    except Exception as e:
        with open('/sd/traceback.txt', 'a') as stdout:
            print_exception(e, stdout)
        log_and_print("Error in message_proc: " + str(e))
        # Problem in message creation, exit
    # Garbage collection
    collect()