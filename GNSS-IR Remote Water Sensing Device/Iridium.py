# Iridium Functions

import micropython
from utime import sleep, sleep_ms
from uasyncio import sleep

# Sending a message -> Requires input of string of hex and length in bytes
# The message should NOT include the 2 byte checksum
def irid_buffer(message, comms_object):
    length = int(len(message) / 2)
    # AT is the attention command -> Run 2x to avoid errors
    comms_object.write('AT\r')

    # Run a while loop that only breaks when the tx is done
    while not comms_object.txdone():
        sleep_ms(1)

    # Read the buffer for the response
    sleep_ms(50)
    r1 = comms_object.read()
    sleep_ms(50)
    
    # Check if recieved OK
    if r1 != None and 'OK' in r1:
        comms_object.write('AT&K0\r')
        # Wait for transmission to be sent
        while not comms_object.txdone():
            sleep_ms(1)
            # Read the buffer for the response
        sleep_ms(50)
        r2 = comms_object.read()
        sleep_ms(50)
        if r2 != None and 'OK' in r2:
            # Write message length to buffer
            string2send = 'AT+SBDWB=' + str(length) + '\r'
            while not comms_object.txdone():
                sleep_ms(1)
            comms_object.write(string2send)
            while not comms_object.txdone():
                sleep_ms(1)
            sleep_ms(50)
            # Read the buffer for the response
            r3 = comms_object.read()
            sleep_ms(50)
            if r3 != None and 'READY' in r3.decode():
                # Generate checksum (2 LSB) and append to end of message
                checksum = gen_checksum(message)
                hex_mess = message + checksum
                ascii_mess = bytearray.fromhex(hex_mess)
                comms_object.write(ascii_mess)
                # Wait for write to be complete
                while not comms_object.txdone():
                    sleep_ms(1)
                sleep_ms(50)
                SBDWB_stat = comms_object.read()
                sleep_ms(50)

                if SBDWB_stat != None and '0' in SBDWB_stat.decode(): #SBDWB status is 0 -> Successfully wrote to ISU
                    final_status = True
            else:
                final_status = False
        else:
            final_status = False
    else:
        final_status = False

    return(final_status)


#Takes in string of hex characters, and returns checksum string
def gen_checksum(string):
    split_strings = [string[i:i+2] for i in range(0, len(string), 2)]
    length = len(split_strings)

    sum = 0
    for x in range (0, length, 1):
        string_new = '0x' + str(split_strings[x])
        val = int(string_new, 16)
        sum += val

    hex_string = str(hex(sum))
    hex_string = hex_string[2:]
    length_hex = len(hex_string)

    if length_hex < 4:
        zeros_needed = 4 - length_hex
        i = 0
        leading_z_string = ''
        while (i<zeros_needed):
            leading_z_string += '0'
            i += 1
        hex_string = leading_z_string + hex_string

    elif length_hex > 4:
        hex_string = hex_string[length_hex-4:]

    return(hex_string)

async def irid_SBD(comms_object):
    comms_object.write('AT\r')
    while not comms_object.txdone():
        sleep_ms(1)
    sleep_ms(50)
    # Set to short burst data transmission mode
    comms_object.write('AT+SBDIX\r')
    while not comms_object.txdone():
        sleep_ms(1)
    print("Initiating SBD Session...")
    # Wait 50 seconds -> Will be async once ready to implement
    await sleep(50)
    final_status = comms_object.read()
    final_status = final_status.decode() #Decode from bytes obj to string
    return(final_status)
    ## Note this leaves the object on -> Need to call status_proc
    ## Otherwise, powering off will clear MT buffer


# Trigger if there is a message in the MT buffer
def irid_read(comms_object):
    comms_object.write('AT\r')
    while not comms_object.txdone():
        sleep_ms(1)
    sleep_ms(50)
    comms_object.write('AT+SBDRB\r') # Read data from MT buffer
    while not comms_object.txdone():
        sleep_ms(1)
    sleep_ms(100)
    message = comms_object.read()
    sleep_ms(100)
    print(message)
    # Trim around message
    start_index = message.find(b'SBDRB\r')
    start_index += len(b'SBDRB\r') + 2 # Omit the command, and the message length
    end_index = message.find(b'\r\n', start_index) - 2 # Omit the checksum and subsequent OK

    message = message[start_index: end_index]

    # Need to convert from hex to integer, but create outlier response for messages that read wrong ??
    return(message)


def SBDIX_status_proc(status):
    # Editing status into list format, preserving only numbers
    numbers_str = ''.join(filter(lambda x: x.isdigit() or x in ('<', '>', ','), status))
    number_strings = numbers_str.split(',')
    status_new = [int(num.strip('<> ')) for num in number_strings if num.strip('<> ').isdigit()]

    # Return the status
    return(status_new)
