def get_cf(sat_type,sigID) : # Returns cf (carrier freq [MHz]) and wl (wavelength [m])
    # Inputs:   sat_type = satallite type i.e. GPGSV = GPS, GAGSV = Gallileo
    #           sigID = signal ID to identify a freq band (usually range 0-8)
    # Outputs:  wl - wavelength [m]

    # cf is an float variable --> carrier frequency [MHz]
    # Initialize the value cf 
    cf = 1575.42 # [MHz] L1 band (default)

    if "GPGSV" in sat_type :
    # Distinguish the bands by signal ID 
    # Refer to https://gpsd.gitlab.io/gpsd/NMEA.html#_nmea_4_11_system_id_and_signal_id
        if sigID in ['0','1','2','3']:
            cf = 1575.42  # [MHz] L1 freq
        elif sigID in ['4','5','6']:
            cf = 1227.6   # [MHz] L2
        elif sigID in ['7','8']:
            cf = 1176.45  # [MHz] L5

    elif "GLGSV" in sat_type: # Glonass
        if sigID in ['0','1','2']:
            cf = 1575.42  # [MHz] L1 freq
        elif sigID in ['3','4']:
            cf = 1227.6   # [MHz] L2

    elif "GAGSV" in sat_type: # Galileo
        if sigID in ['6','7']:
            cf = 1575.42  # [MHz] L1 freq 
        elif sigID in ['1']:
            cf = 1176.45  # [MHz] E5a/L5 band
        elif sigID in ['2']:
            cf = 1207.14  # [MHz] E5b/L5 band
        elif sigID in ['3']:
            cf = 1191.795 # [MHz] E5a+b/L5 band
        elif sigID in ['4','5']:
            cf = 1278.75 # [MHz] E6

    elif ("GBGSV" in sat_type or "BDGSV" in sat_type): # BeiDou 
        if sigID in ['0','1','2']:
            cf = 1561.098   # [MHz] B1I/B1Q
        elif sigID in ['3','4']:
            cf = 1575.42    # [MHz] B1C/B1A
        elif sigID in ['5']:
            cf = 1176.45    # [MHz] B2a
        elif sigID in ['6','B','C']:
            cf = 1207.14    # [MHz] B2b
        elif sigID in ['7']:
            cf = 1191.795   # [MHz] B2a+b
        elif sigID in ['8']:
            cf = 1268.52    # [MHz] B3I/B3Q/B3A

    elif "GQGSV" in sat_type: # QZSS
        if sigID in ['0','1','2','3']:
            cf = 1575.42  # [MHz] L1 freq 
        elif sigID in ['5','6']:
            cf = 1227.6   # [MHz] L2
        elif sigID in ['7','8']:
            cf = 1176.45  # [MHz] L5
        elif sigID in ['9','A'] :
            cf = 1278.75  # [MHz] L6



    c = 299792458       # [m/s] speed of light
    wl = c/(cf*10**6)   # [m] wavelength for a cf

    return wl