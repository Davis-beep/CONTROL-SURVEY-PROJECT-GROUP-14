import math

# Function to round half up
def round_half_up(x, ndigits=0):
    multiplier = 10 ** ndigits
    if x >= 0:
        return math.floor(x * multiplier + 0.5) / multiplier
    else:
        return math.ceil(x * multiplier - 0.5) / multiplier

# Function to convert DMS to seconds
def dms_to_seconds(d, m, s):
    return d * 3600.0 + m * 60.0 + s

# Function to convert seconds to DMS
def seconds_to_dms(total_seconds):
    sign = -1 if total_seconds < 0 else 1
    total_seconds = abs(total_seconds)
    d = int(total_seconds // 3600)
    total_seconds %= 3600
    m = int(total_seconds // 60)
    s = total_seconds % 60
    s = round_half_up(s, 0)
    return sign * d, m, s

# Function to input DMS
def input_dms(label):
    d = int(input(f"{label} Degrees: "))
    m = int(input(f"{label} Minutes: "))
    s = float(input(f"{label} Seconds: "))
    return d, m, s

# Function to format DMS
def format_dms(d, m, s):
    return f"{abs(d)}° {m}' {int(s):02d}\""

                                
def main():
    print("=== TRIANGULATION BOOKING SHEET ===")
    
    print("===NB USER IS SUPPOSED TO INPUT THE CORRECTED RIGHT BACK BEARING===")
    print("===WELCOME= ==")
    
    observer = input("Observer Name: ")
    
    occupied_station = input("Occupied Station: ")
    rounds = int(input("Number of Rounds: "))
    num_stations = int(input("Number of Observed Stations: "))

    station_names = []
    for i in range(1, num_stations + 1):
        name = input(f"Name of station {i}: ")
        station_names.append(name)

    first_round_first_station_avg = None

    for r in range(1, rounds + 1):
        print(f"\n--- ROUND {r} ---")
        avg_readings = []
        for i, station in enumerate(station_names, start=1):
            print(f"\nStation {i}: {station}")
            print("Face Left (FL):")
            FL = input_dms("FL")
            print("Face Right (FR):")
            FR = input_dms("FR")

            FL_sec = dms_to_seconds(*FL)
            FR_sec = dms_to_seconds(*FR)
            avg_sec = round_half_up((FL_sec + FR_sec) / 2.0, 0)
            avg_readings.append(avg_sec)

        first_avg = avg_readings[0]
        last_avg = avg_readings[-1]
        misclosure = first_avg - last_avg
        e = misclosure / (num_stations - 1) if num_stations > 1 else 0
        print(f"\nMisclosure = {round_half_up(misclosure, 0)} seconds")
        print(f"Orientation correction (e) = {round_half_up(e, 3)} seconds")

        corrected_readings = []
        for i, avg_val in enumerate(avg_readings, start=1):
            correction = e * (i - 1)
            corrected_val = avg_val + correction
            corrected_val = round_half_up(corrected_val, 0)
            corrected_readings.append(corrected_val)

        station_error_corrs = []
        if r == 1:
            first_round_first_station_avg = corrected_readings[0]
            station_error_corrs = [0 for _ in corrected_readings]
        else:
            station_error_value = first_round_first_station_avg - corrected_readings[0]
            station_error_corrs = [station_error_value for _ in corrected_readings]
            for i in range(len(corrected_readings)):
                corrected_readings[i] = round_half_up(corrected_readings[i] + station_error_value, 0)

        print("\n--- STATION READINGS SUMMARY ---")
        print(f"{'Station':<12} {'Average Reading':<20} {'Correction (e)':<15} {'Corrected Reading':<20} {'Station Error Corr (DMS)':<25}")
        print("-" * 100)
        for idx, (station, avg_val, corr_val, err_corr) in enumerate(zip(station_names, avg_readings, corrected_readings, station_error_corrs), start=1):
            d_a, m_a, s_a = seconds_to_dms(avg_val)
            d_c, m_c, s_c = seconds_to_dms(corr_val)
            d_err, m_err, s_err = seconds_to_dms(err_corr)
            correction_display = e * (idx - 1)
            print(f"{station:<12} {format_dms(d_a, m_a, s_a):<20} {f'{correction_display:+.3f}':<15} {format_dms(d_c, m_c, s_c):<20} {format_dms(d_err, m_err, s_err):<25}")
        print("\n=== COMPUTATION COMPLETE ===")
        print(f"Observer: {observer} | Occupied Station: {occupied_station}") 
      
if __name__ == "__main__":
    main()