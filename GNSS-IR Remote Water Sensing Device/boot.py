# sdcard class from micropythonlib
from sdcard import SDCard
# Machine is a standard micropython library
from machine import Pin, SPI
# Python basic os functions
from os import mount, VfsFat

# Assign CS pin (and start it high) to GPIO 9 (pin 12)
cs = Pin(9, Pin.OUT)
# Assign MOSI-TX pin to GPIO 11 (pin 15)
mosiPin = Pin(11)
# Assign MISO->RX pin to GPIO 8 (pin 11)
misoPin = Pin(8)
# Assign the clock signal to GPIO 10 (pin 14)
sckPin = Pin(10)

# Intialize SPI peripheral (start with 1 MHz)
# use SPI bus 1: MISO 16, MOSI 15, CLK 14, CS17
spi = SPI(1, baudrate=1320000, polarity=0, phase=0, bits=8, firstbit=SPI.MSB, sck=sckPin, mosi=mosiPin, miso=misoPin)

# Initialize SD card
sd = SDCard(spi, cs)

# Mount filesystem
vfs = VfsFat(sd)
#uos.umount("/sd")
mount(vfs, "/sd")
