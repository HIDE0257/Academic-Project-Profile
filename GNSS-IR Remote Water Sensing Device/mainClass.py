# Machine is a standard micropython library
from machine import Pin, UART, RTC, ADC, I2C
from utime import time, sleep
from math import floor
from random import random
# Pressure sensor class
import sht31
from press_data import take_press_data
# Reflection height
from Reflection_Height import reflection_height as RH
# Iridium functions
from Iridium import irid_buffer, irid_read, irid_SBD, SBDIX_status_proc
# persistent reset
from reset_m import persistent, nonper
# Message processing
from message import message_proc
# Logging fn
from log import log_and_print
# Asyncio I/O event scheduler
import uasyncio as asyncio
# AS_GPS Package from Peter Hinch
from as_GPS import AS_GPS
# Basic OS functions
from os import statvfs, rename, mkdir, remove, listdir
from sys import exit, print_exception
# Garbage collector
from gc import enable, collect, mem_free
# carrier frequency function
from get_cf import get_cf

class mainClass:
    # Initialization function
    def __init__(self):

        # Turn on the LED to indicate entry into the initialization
        # make an instance of the pin class from machine to use our LED
        LEDpin = Pin("LED", Pin.OUT)
        LEDpin.toggle()

        # Enable garbage collection if it isn't already
        enable()

        # Initialize UART connection to ublox NEO-F10N
        # default baud rate: 38400, use UART port 1
        uart_gnss = UART(1, 38400)
        # Use an asyncio stream reader/writer
        gnss_sreader = asyncio.StreamReader(uart_gnss)
        # Instantiate the asynchgps driver
        self.gnss = AS_GPS(gnss_sreader)

        # Initialize the uart object for communiations
        self.uartComms = UART(0, 19200, tx=Pin(0), rx=Pin(1))
        self.uartComms.init(bits=8, parity=None, stop=1)

        # Verify the file system
        vfsStat = statvfs("/sd")
        # Check if the block size of the sd card is as expected
        if(vfsStat[0] != 16384):
            # If we reach here it's reasonable to assume the mount failed
            print_exception("SD Card Mount failed (block size is not that of expected), please check if one is inserted and try again")
            exit()

        dir = listdir("/sd")

        # Azimuth and Elevation Angle Mask Handling
        # Check if there is an azimuth config already existing, open the file
        if "mask.csv" in dir:
            with open("/sd/mask.csv", "r") as file:
                # Attempt to read
                masks = file.readlines()
                # Check if this input is null or not expected
                if len(masks) != 3:
                    print("Writing default masks to file")
                    # The information does not exist, so set to the default values
                    self.aziStart = 0
                    self.aziEnd = 360
                    self.eleStart = 0
                    self.eleEnd = 45
                    self.minHeight = 0
                    self.maxHeight = 5
                    # Write to the expected file location
                    with open("/sd/mask.csv", "w") as file2:
                        file2.write("0,360\n0,45\n0,5")
                else:
                    # Otherwise, pull the information from the files
                    azi = masks[0].split(",")
                    ele = masks[1].split(",")
                    hei = masks[2].split(",")
                    # Apply to the variables
                    self.aziStart = int(azi[0])
                    self.aziEnd = int(azi[1])
                    self.eleStart = int(ele[0])
                    self.eleEnd = int(ele[1])
                    self.minHeight = float(hei[0])
                    self.maxHeight = float(hei[1])
                    del azi, ele, hei
                del masks
        else:
            # Create the file and write the default 
            with open("/sd/mask.csv", "w") as file:
                file.write("0,360\n0,45\n0,5")
            self.aziStart = 0
            self.aziEnd = 360
            self.eleStart = 0
            self.eleEnd = 45
            self.minHeight = 0
            self.maxHeight = 5

        print(self.aziStart)
        print(self.aziEnd)
        print(self.eleStart)
        print(self.eleEnd)
        print(self.minHeight)
        print(self.maxHeight)

        # Downsampling handling
        # Check if there is a onfig already existing, open the file
        if "down.csv" in dir:
            with open("/sd/down.csv", "r") as file:
                # Attempt to read
                samp = file.readlines()
                # Check if this input is null or not expected
                if len(samp) != 1:
                    print("Writing default downsampling to file")
                    # The information does not exist, so set to the default values
                    self.sampleRate = 13
                    self.closeoutTime = 2700
                    # Write to the expected file location
                    with open("/sd/down.csv", "w") as file2:
                        file2.write("13,2700")
                else:
                    # Otherwise, pull the information from the files
                    samp = samp[0].split(",")
                    # Apply to the variables
                    self.sampleRate = int(samp[0])
                    self.closeoutTime = int(samp[1])
                del samp
        else:
            # Create the file and write the default 
            with open("/sd/down.csv", "w") as file:
                file.write("13,2700")
            self.sampleRate = 13
            self.closeoutTime = 2700

        # Preallocation size
        self.preAlloc = floor(self.closeoutTime / (2+self.sampleRate))

        # Create the needed directories
        if not "MO" in dir:
            mkdir("/sd/MO")

        if not "MT" in dir:
            mkdir("/sd/MT")

        if not "hold" in dir:
            mkdir("/sd/hold")

        # Initialize the message queue
        self.msgQueue = ""
        # Check if there is a persistent msg file currently
        if "msg.txt" in dir:
            with open("/sd/msg.txt", "r") as file:
                # Attempt to read
                msg = file.readlines()
                # Add into message queue
                if len(msg) != 0:
                    self.msgQueue += msg[0]
            # Delete the file
            remove("/sd/msg.txt")
            del msg
        
        # Check if this startup was a reboot
        if "reboot" in dir:
            # If it's present, this was a reboot, so add the response to the msg queue
            self.msgQueue += "04"
            # Delete the file
            remove("/sd/reboot")

        # health message booleans
        self.healthBool = False
        self.tempBool = False
        self.procBool = False

        # Lock to only close out a single curve at a time (memory limitation) and sending a measurement
        self.closeout_lock = asyncio.Lock()

        # Event to handle full system reboot
        self.reboot_event = asyncio.Event()

        # Instantiate the pin for the comms module
        self.comms_power_pin = Pin(3, Pin.OUT)

        # Configure the ADC Pin for the on-board temperature sensor
        self.picoTemp = ADC(4)

        # Configure I2C for the sht31
        self.temp_sensor = sht31.SHT31(I2C(0, scl = Pin(21), sda = Pin(20), freq=400000), addr=0x44)

        # Configure I2C for pressure sensor
        self.pressure_sensor = I2C(1, scl = Pin(15), sda = Pin(14), freq=400000)

        # Initialize a list to hold the PRN's currently being searched
        self.prn_list = []
        # Check if there is a persistent prn file currently
        if "prn.csv" in dir:
            with open("/sd/prn.csv", "r") as file:
                # Attempt to read
                prns = file.readlines()
                # Separate into a list of tuples
                for prn in prns:
                    # Split the prn
                    prn = prn.split(',')
                    self.prn_list.append( (prn[0],int(prn[1].replace('\n', ''))) )
                    # Ensure the new closeout time is >= 0
                    closeoutMod = self.closeoutTime - (time() - int(prn[1].replace('\n', '')))
                    # Re-initialize the curve closeouts for these persistent prns
                    asyncio.create_task(self.closeout_curve(prn[0], closeoutMod ) )
            # Delete the persistent prn file
            remove("/sd/prn.csv")
            del prns

        del dir
        # Turn off the LED to inform that setup is complete
        LEDpin.toggle()

    # Mask updater function
    def updateMask(self, new_mask):
        # Ensure the masks are within valid ranges
        if new_mask[0] >= 0 and new_mask[1] <= 360 and new_mask[2] >= 0 and new_mask[3] <= 45 and new_mask[4] >= 0 and new_mask[5] <= 50:
            try:
                # Set the local variable masks to themselves
                self.aziStart = new_mask[0]
                self.aziEnd = new_mask[1]
                self.eleStart = new_mask[2]
                self.eleEnd = new_mask[3]
                self.minHeight = new_mask[4]
                self.maxHeight = new_mask[5]

                # Write to file
                with open("/sd/mask.csv", "w") as file:
                    file.write(str(new_mask[0]) + ',' + str(new_mask[1]) + '\n' + str(new_mask[2]) + ',' + str(new_mask[3]) + '\n' + str(new_mask[4]) + ',' + str(new_mask[5]))

                # Enqueue a message to respond to the mask update
                self.msgQueue += "0500"

                # Set an event for a persistent reset to change the sample rate TODO

            except:
                # Some form of error in mask update, send an error msg
                self.msgQueue += "0502"
        else:
            # Send a bad mask msg
            self.msgQueue += "0501"

        del new_mask
        collect()

    # sampling updater function
    def updateSample(self, new_sampleRate, new_closeoutTime):
        # Compute number of datapoints
        preAlloc = floor(new_closeoutTime / (2+new_sampleRate))
        # Ensure the sampling info is within ranges
        if preAlloc <= 216:
            try:
                # Set the local variable masks to themselves
                self.closeoutTime = new_closeoutTime
                self.sampleRate = new_sampleRate
                self.preAlloc = preAlloc

                # Write to file
                with open("/sd/down.csv", "w") as file:
                    file.write(str(new_sampleRate) + ',' + str(new_closeoutTime))

                # Enqueue a message to respond to the sampling update
                self.msgQueue += "0600"
            except:
                # Some form of error in sampling update, send an error msg
                self.msgQueue += "0602"
        else:
            # Send a bad sampling info msg
            self.msgQueue += "0601"

        del new_sampleRate, new_closeoutTime, preAlloc
        collect()

    # GPS information handler function
    async def gnss_handler(self):
        while True:
            # Wait for the lock
            await self.closeout_lock.acquire()
            # Await recieving GSV message data from the gnss chip
            gsv_dict = await self.gnss.get_satellite_data()

            try:
                # Send to the function that handles saving this data to the sd card
                # the gnss.values
                self.write_gnss(time(), gsv_dict.values())
            except Exception as e:
                with open('/sd/traceback.txt', 'a') as stdout:
                    print_exception(e, stdout)
                log_and_print("Error in gnss_handler: " + str(e))
                # Force a reset as a panic if this is ever an error
                nonper()

            # release the lock
            self.closeout_lock.release()
            collect()
            # Sleep the sample rate time and collect garbage
            await asyncio.sleep(self.sampleRate)

    # Function that closes out and processes an SNR curve
    async def closeout_curve(self, idString, closeoutTimeModif):
        # Wait configured time from being called
        await asyncio.sleep(self.closeoutTime - closeoutTimeModif)
        # Wait to recieve the lock to allow execution
        await self.closeout_lock.acquire()
        collect()
        log_and_print("Closing out curve " + idString + ", current free memory: " + str(mem_free()))
        mem_open_Bool = False
        for i in range(15):
            log_and_print(str(mem_free()))
            if int(mem_free()) < 160000:
                mem_open_Bool = True
                break
            collect()
            await asyncio.sleep(1)
            
        if mem_open_Bool:
            ele_correction_bool = False
            try:
                # Read tempture and himidity from sht31 sensor
                X = self.temp_sensor.get_temp_humi()
                temp_C = X[0]
                hum = X[1]
                press = take_press_data(self.pressure_sensor, 0x5D)
                ele_correction_bool = True
            except Exception as e:
                # Don't run the elevation correction if something is wrong
                log_and_print("Error in closeout_curve environment measurement: " + str(e))
                ele_correction_bool = False
                temp_C = 0
                hum = 0
                press = 0
            # Force perform a garbage collection
            collect()

            try:
                # From here, process the data in the file
                # Initialize the long lists, preallocating them
                times = [None] * self.preAlloc
                eles = [None] * self.preAlloc
                azis = [None] * self.preAlloc
                snrs = [None] * self.preAlloc
                # Extract it
                with open("/sd/" + idString + ".csv", "r") as file:
                    # Read the first line to get the information about this satellitee
                    satInfo = (file.readline()).split(",")
                    sigID = satInfo[0]
                    msgTp = satInfo[2]
                    t_0 = int(satInfo[3])

                    # get the signal type from the function
                    wavelength = get_cf(msgTp,sigID)
                    del msgTp, sigID, satInfo

                    # Read the rest of the ele, azi, and snr data
                    for i in range(self.preAlloc):
                        row = file.readline()
                        # Ensure we have not reached the end of the file or we aren't about to overallocate a list
                        if not row:
                            break
                        rowS = row.split(",")
                        del row
                        # Append to the lists
                        times[i] = int(rowS[0])
                        eles[i] = int(rowS[1])
                        azis[i] = int(rowS[2])
                        snrs[i] = int(rowS[3].replace('\n', ''))

                collect()
                # Once we are done, remove all elements that are 'none' from the lists
                times = [i for i in times if i is not None]
                eles = [i for i in eles if i is not None]
                azis = [i for i in azis if i is not None]
                snrs = [i for i in snrs if i is not None]

                # Garbage collection
                collect()
                
                # Only pass to the function if there is more than 1 datapoint
                if len(times) > 1:
                    # Pass to function to compute h_bar and stuff
                    t, h_bar, e = RH(wavelength, times, eles, azis, snrs, self.minHeight, self.maxHeight, self.aziStart, self.aziEnd, self.eleStart, self.eleEnd, ele_correction_bool, press, temp_C, hum)
                    collect()
                    # If the output is 0,0,0, do nothing
                    if t != 0 and h_bar != 0 and e != 0:
                        # Write into file for h_bar datapoints
                        with open("/sd/hbar.csv", "a") as file:
                            file.write(str(t_0 + int(t)) + "," + str(h_bar) + "," + str(e) + "\n")

            except Exception as e:
                with open('/sd/traceback.txt', 'a') as stdout:
                    print_exception(e, stdout)
                # On any given error, make sure to still release the lock to allow others to continue
                log_and_print("Error in closeout_curve for " + idString + ": " + str(e))
                # add processing error message
                if not self.procBool:
                    self.msgQueue += '0301'
                    self.healthBool = True
                    self.procBool = True
        else:
            log_and_print("Could not clear enogh memory to process curve safely!")

        # Clear memory
        del times, eles, azis, snrs
        collect()

        # Release the lock to allow another task to run
        self.closeout_lock.release()

        # Move the file for safekeeping later
        # List the directory to search
        # dir = listdir("/sd/hold")
        # count = 1
        # while True:
        #     collect()
        #     if idString + "_" + str(count) + ".csv" in dir:
        #         count += 1
        #     else:
        #         # Place the file in
        #         rename("/sd/" + idString + ".csv", "/sd/hold/" + idString + "_" + str(count) + ".csv")
        #         print("Successfully moved file with footer " + str(count))
        #         break

        # del dir, count
        # collect()

        # Once that's done, try to delete the file just in case
        try:
            remove("/sd/" + idString + ".csv")
        except:
            pass

        # Remove the file from the current PRN queue
        self.prn_list = [i for i in self.prn_list if i[0] != idString]
                

    # Function that stores the GNSS information we recieve onto the SD card
    def write_gnss(self, time, gsv_dict):
        
        # Extract the prn numbers of each dict and see if they are currently writing to a curve
        for gsv in gsv_dict:
            collect()
            sigID = gsv[1]
            satPRN = gsv[2]
            
            # Assemble into the standard string
            idString = str(satPRN) + "_" + str(sigID)
            # Check the list of prns being searched
            PRN_info = [i for i in self.prn_list if i[0] == idString]
            if len(PRN_info) != 0:
                start_time = PRN_info[0][1]
                del PRN_info
                
                # If this works, then the file already exists, so APPEND the new data onto it
                with open("/sd/" + idString + ".csv", "a") as file:
                    # Check if any of the values are none, and do nothing if they are
                    ele = gsv[3]
                    azi = gsv[4]
                    snr = gsv[5]
                    if(ele != None and azi != None and snr != None):
                        file.write(str(time - start_time) + "," + str(ele) + "," + str(azi) + "," + str(snr) + "\n")

            # If it's not in there,
            else:
                del PRN_info
                # Check if the PRN list has too many members
                #if len(self.prn_list) <= 12:
                # Initialize a new file to store the data over time (or overwrite an old once once it's done)
                with open("/sd/" + idString + ".csv", "w") as file:
                    ele = gsv[3]
                    azi = gsv[4]
                    snr = gsv[5]
                    msgTp = gsv[0]
                    # The first line will contain the signal id, PRN number, and message type
                    file.write(str(sigID) + "," + str(satPRN) + "," + str(msgTp) + "," + str(time) + "\n")
                    # The following lines will be time, ele, azi, snr
                    if(ele != None and azi != None and snr != None):
                        file.write("0," + str(ele) + "," + str(azi) + "," + str(snr) + "\n")
                # Add the prn number into the list
                self.prn_list.append( (idString, time) )
                # Create the task that closes this file out after 1hr
                asyncio.create_task(self.closeout_curve(idString, 0))
                
        # Clear memory
        del gsv_dict, time

    async def cleanup(self):
        await asyncio.sleep(21900)
        await self.closeout_lock.acquire()
        log_and_print("Beginning 6hr Cleanup")

        # Perform a garbage collection
        collect()

        # Delete measurements older than 3hrs
        data = []
        # Get the reference time to delete everything
        timeToDel = time() - 10800
        # Open the relevant files by checking the range
        if "hbar.csv" in listdir("/sd"):
            with open("/sd/hbar.csv", "r") as file:
                # Start reading all the rows and place them into individual members of a list
                allRows = file.readlines()
                # Check each row individually
                for row in allRows:
                    # Check if the time is within 3 hrs of our current time
                    rowCont = row.split(",")
                    dpTime = int(rowCont[0])
                    if dpTime > timeToDel:
                        data.append(row)

            # Overwrite the file with all of the data measurements we want to keep
            with open("/sd/hbar.csv", "w") as file:
                # Write each row as the data in 'data'
                for dp in data:
                    file.write(dp)
        collect()

        # Move the logfile, MO msges, and MT msges to hold
        dir = listdir("/sd/hold")
        if "MO.csv" in listdir("/sd/MO"):
            # Move into the hold folder       
            count = 1
            while True:
                if "MO_" + str(count) + ".csv" in dir:
                    count += 1
                else:
                    # Place the file in
                    rename("/sd/MO/MO.csv", "/sd/hold/MO_" + str(count) + ".csv")
                    break
        collect()

        if "MT.csv" in listdir("/sd/MT"):
            # Move into the hold folder       
            count = 1
            while True:
                if "MT_" + str(count) + ".csv" in dir:
                    count += 1
                else:
                    # Place the file in
                    rename("/sd/MT/MT.csv", "/sd/hold/MT_" + str(count) + ".csv")
                    break
        collect()

        if "log.txt" in listdir("/sd"):
            # Move into the hold folder       
            count = 1
            while True:
                if "log_" + str(count) + ".txt" in dir:
                    count += 1
                else:
                    # Place the file in
                    rename("/sd/log.txt", "/sd/hold/log_" + str(count) + ".txt")
                    break
        collect()

        # Delete the traceback file
        try:
            remove("/sd/traceback.txt")
        except:
            pass

        # Call a persistent reset once we are done
        collect()
        log_and_print("6hr Cleanup Complete, restarting")
        persistent(self.prn_list, self.msgQueue)
        

    # Function to enqueue a water level measurement
    def queueWaterLevelMessage(self, height, measTime):
        # Generate the water level in the message format
        h_oneten = floor(height)
        h_tenthshundths = floor( (height - h_oneten)*100)
        # Convert to integers
        h_oneten = int(h_oneten)
        h_tenthshundths = int(h_tenthshundths)
        # Convert to hex
        timeHex = hex(measTime)
        onetenHex = hex(h_oneten)
        tenthshundthsHex = hex(h_tenthshundths)
        # Get rid of the 0x footer from all
        timeHex = timeHex[2:]
        onetenHex = onetenHex[2:]
        tenthshundthsHex = tenthshundthsHex[2:]
        # Add on the appropriate amount of leading zeros
        timeHex = ('00000000' + timeHex)[-8:]
        onetenHex = ('00' + onetenHex)[-2:]
        tenthshundthsHex = ('00' + tenthshundthsHex)[-2:]

        # combine along with header into the full message
        msg = "02" + timeHex + onetenHex + tenthshundthsHex

        # Add onto the message queue
        self.msgQueue = self.msgQueue + msg

        del height, measTime


    # Function to generate a datapoint for the last 30 minutes
    async def processLastThirty(self):
        while True:
            # Collect all of the datapoints from the last 30 minutes
            # await for 30 minutes before attempting to do anything
            await asyncio.sleep(1801)

            # Get current time
            timeToRead = time()

            # Create lists for the data
            times = []
            hbars = []
            edots = []
            # Check if there is an hbar file
            if "hbar.csv" in listdir("/sd"): 
                # check the h file for height measurements from the past hour
                with open("/sd/hbar.csv", "r") as file:
                    # split out the rows
                    allRows = file.readlines()

                    # For every row, get the times and datapoints
                    for row in allRows:
                        rowS = row.split(",")
                        # Make sure it is a measurement from within the last hour
                        dptime = int(rowS[0])
                        if (dptime > timeToRead - 3600):
                            # Append to the lists
                            times.append(dptime)
                            hbars.append(float(rowS[1]))
                            edots.append(float(rowS[2]))

                # make sure there is acutally data
                if len(times) != 0:
                    # With all the datapoints, pass to the function to run the dynamic height correction TODO
                    #h = DHC(times, hbars, es, edots)
                    # For now just take an average of each
                    h = sum(hbars) / len(hbars)
                    # Add this value to the message queue to be sent for the time point 30 minutes ago
                    self.queueWaterLevelMessage(h, timeToRead - 1800)
                    log_and_print("processLastThirty: New measurement " + str(h))
            collect()

    # Function that will attempt to send a message every hour
    async def sendMOMsg(self):
        while True:
            # Wait 1 hr to send a message
            await asyncio.sleep(3605)
            #log_and_print("Sending MO Message")

            # Read the message queue, and don't try to send a message if it's empty
            # Check if a health status message is in the message queue
            if self.healthBool and len(self.msgQueue) != 0:
                # Send a message
                asyncio.create_task(self.irid_send(self.msgQueue))
            else:
                # Send the health status message for 'all okay'
                self.msgQueue += "0300"
                asyncio.create_task(self.irid_send(self.msgQueue))
            # Reset all the boolean to false to allow new health status messages
            self.healthBool = False
            self.tempBool = False
        
            # Store the MO message to the SD Card
            with open("/sd/MO/MO.csv", "a") as file:
                file.write(str(time()) + "," + self.msgQueue + "\n")

            # Empty the message queue
            self.msgQueue = ""

    # Async function that blinks the LED as a status indicator
    async def blink(self):
        # make an instance of the pin class from machine to use our LED
        LEDpin = Pin("LED", Pin.OUT)

        # Toggle the LED every one second
        while True:
            LEDpin.toggle()
            await asyncio.sleep(1)


    # Wrapper function for iridium
    async def irid_send(self, msg_to_send):

        # Wait to recieve the lock to allow execution
        await self.closeout_lock.acquire()

        try:

            # Enable the power to the reciever and give it a moment to turn on
            self.comms_power_pin.value(1)
            sleep(1)

            # Buffer the message, retrying 10 times if it doesnt work
            count = 0
            buffer_status = irid_buffer(msg_to_send, self.uartComms)
            while not buffer_status:
                if count > 10:
                    raise ValueError
                # Try again, give up after 10 attempts
                buffer_status = irid_buffer(msg_to_send, self.uartComms)
                count += 1
            del count

            # Queue up a final status task
            final_status = await asyncio.create_task(irid_SBD(self.uartComms))
            #log_and_print(final_status)
            # Status processing
            status = SBDIX_status_proc(final_status)

            # Check the status of the message to decide what to do next
            SBD_bool = False
            if status[0] <= 2:
                log_and_print("Message send success, status " + str(status[0]))
            elif status[0] == 32:
                #log_and_print("Message send failure, status 32: No Network Service. Retrying")
                SBD_bool = True
            else:
                log_and_print("Unaccounted message send status " + str(status[0]))

            # If the retry boolean is true, enter the retry sequence
            retry_count = 0
            wait_time = 0
            while SBD_bool:
                if retry_count < 2:
                    wait_time = random() * 5
                elif retry_count < 4:
                    wait_time = random() * 30
                else:
                    # Power the comms module off and on, then keep trying
                    self.comms_power_pin.value(0)
                    sleep(1)
                    self.comms_power_pin.value(1)
                    sleep(1)
                    wait_time = 60
                    # Re-queue the message after the power off
                    # Buffer the message, retrying 10 times if it doesnt work
                    count = 0
                    buffer_status = irid_buffer(msg_to_send, self.uartComms)
                    while not buffer_status:
                        if count > 10:
                            raise ValueError
                        # Try again, give up after 10 attempts
                        buffer_status = irid_buffer(msg_to_send, self.uartComms)
                        count += 1
                    del count
                    
                # Iniate another SBD Session after the wait time
                await asyncio.sleep(wait_time)
                final_status = await asyncio.create_task(irid_SBD(self.uartComms))
                status = SBDIX_status_proc(final_status)
                # Based on the output, change the boolean to false if the message was successful
                if status[0] <= 2:
                    log_and_print("Message successful after retry " + str(retry_count))
                    SBD_bool = False
                else:
                    retry_count += 1
                    #log_and_print("Message retry failed, attempting again. Attempt number " + str(retry_count))

            # Mailbox check
            MT_bool = False
            if status[2] == 0:
                MT_bool = False
            elif status[2] == 1:
                MT_bool = True
            else:
                MT_bool = False
                #log_and_print("Error in mailbox check")

            # If we have recieved a message, download it
            if MT_bool:
                # Download this MT Message
                message = irid_read(self.uartComms)
                # Convert the byte array into a string
                message = ''.join('{:02x}'.format(x) for x in message)
                # Store MT messsage to the device
                with open("/sd/MT/MT.csv", "a") as file:
                    file.write(str(time()) + "," + message + "\n")

                # Interperet the message and execute it's orders
                message_proc(message, self.reboot_event, self.updateSample, self.updateMask)

            # All is said and done for transmitting. Check the number of messages to read in the MT buffer and read them all
            while status[5] > 0:
                final_status = await asyncio.create_task(irid_SBD(self.uartComms))
                status = SBDIX_status_proc(final_status)
                # Download this MT Message
                message = irid_read(self.uartComms)
                # Convert the byte array into a string
                message = ''.join('{:02x}'.format(x) for x in message)
                # Store MT messsage to the device
                with open("/sd/MT/MT.csv", "a") as file:
                    file.write(str(time()) + "," + message + "\n")

                # Interperet the message and execute it's orders
                message_proc(message, self.reboot_event, self.updateSample, self.updateMask)
        except Exception as e:
            with open('/sd/traceback.txt', 'a') as stdout:
                print_exception(e, stdout)
            log_and_print("Error in irid_send: " + str(e))


        # Disable the power to the comms module
        self.comms_power_pin.value(0)
        sleep(1)

        # Release the lock to allow another task to run
        self.closeout_lock.release()

    async def readIntTemp(self):
        # Wait 15 minutes between measurements
        await asyncio.sleep(900)

        # Read the temperature
        bin_val = self.picoTemp.read_u16() # read sensor voltage bin value
        Voltage = (3.3/65535) * bin_val # determine foltage from bin value
        temp = 23.7 - (Voltage - .724)/.001721 # determine temp from voltage [degC]

        if temp > 55 or temp < -20:
            self.healthBool = True
            # Only add a message if there isnt already a temperature bool
            if not self.tempBool:
                self.msgQueue += "0301"
                self.tempBool = True

        # Create another health status for later on
        asyncio.create_task(self.readIntTemp())

    async def mainLoop(self):

        # Create a task that updates the RTC before doing anything else
        # Wait to recieve a GNSS time
        await self.gnss.data_received(date = True)
        hms = self.gnss.utc
        dmy = self.gnss.date
            
        # Update the real time clock to what we got from the GNSS constellation
        RTC().datetime((2000 + dmy[2], dmy[1], dmy[0], 0, hms[0], hms[1], hms[2], 000))
        
        # Create a global epoch
        self.epoch = time()
        log_and_print("<- The epoch")

        # Make a task to recieve the GNSS data (which will recreate itself)
        asyncio.create_task(self.gnss_handler())

        # Make the task that handles the 6hr cleanup
        asyncio.create_task(self.cleanup())

        # Make the task that blinks the LED
        asyncio.create_task(self.blink())

        # Make a task for the message sender
        asyncio.create_task(self.sendMOMsg())

        # Make a task to process the last 30min of data points
        asyncio.create_task(self.processLastThirty())

        # Task that reads the internal temperature sensor
        asyncio.create_task(self.readIntTemp())
        
        log_and_print("GNSS Handler, Message Sender, and 6hr Cleanup initialized")

        # Start an infinite loop! yipee
        while True:
            # Break the loop if we see the loop breaking event
            await self.reboot_event.wait()
            # Call a non-persistent reset
            nonper()
            break
                
    

