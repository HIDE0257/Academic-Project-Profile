import micropython
import utime as time

def log_and_print(msg):
    # Save message to a log file along with printing to console
    t = time.time()
    print(str(t) + ": " + str(msg))

    with open("/sd/log.txt", "a") as file:
        file.write(str(t) + ": " + str(msg) + "\n")
        file.close()