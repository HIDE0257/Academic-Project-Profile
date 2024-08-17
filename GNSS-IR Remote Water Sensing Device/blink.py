# Machine is a standard micropython library
from machine import Pin
# utime from micropython
from utime import sleep

# make an instance of the pin class from machine to use our LED
LEDpin = Pin("LED", Pin.OUT)

# Toggle the LED every one second
while True:
    LEDpin.toggle()
    sleep(1)