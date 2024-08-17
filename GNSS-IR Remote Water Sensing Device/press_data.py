# Taking pressure data from external sensor via I2C Bus

# Importing libraries
from utime import sleep_ms

# Call the following function every 30 minutes -> Should take in new data and overwrite old measurement
# Takes in I2C object and device address
def take_press_data(i2c, press_addr):

    # Set start condition 
    i2c.writeto_mem(press_addr, 0x11, '11')
    sleep_ms(100)
    
    # Pull pressure data and remove everything but 2 digit hex value

    press_out_xl = i2c.readfrom_mem(press_addr, 0x28, 1)
    press_out_xl = ''.join('{:02x}'.format(byte) for byte in press_out_xl)
    press_out_l = i2c.readfrom_mem(press_addr, 0x29, 1)
    press_out_l = ''.join('{:02x}'.format(byte) for byte in press_out_l)
    press_out_h = i2c.readfrom_mem(press_addr, 0x2a, 1)
    press_out_h = ''.join('{:02x}'.format(byte) for byte in press_out_h)

    #Create 24 bit "word"
    press_hex = press_out_h + press_out_l + press_out_xl

    bin = int(press_hex, 16)

    if (bin & (1 << (32 - 1))): # if sign bit is set e.g., 8bit: 128-255
        bin = bin - (1 << 32)        # compute negative value            

    # Divide by 4096 to find press value in hpa
    press = bin / 4096  # [hPa]

    # offset = 270 # Seems absurdly high but whatever idk
    # press = press + offset
    return press
