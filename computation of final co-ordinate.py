import math
import numpy as np
import matplotlib.pyplot as plt

def dms_to_decimal(d, m, s):
    return d + m/60.0 + s/3600.0

def get_float(p):
    return float(input(p))

def get_dms(p):
    d, m, s = map(float, input(f"{p} (D M S): ").split())
    return dms_to_decimal(d, m, s)

def compute_provisional():
    print("\n=== PROVISIONAL COORDINATES ===")
    name = input("Station name: ").strip() or "P"
    
    print("\nStation A:"); N_A = get_float("  N: "); E_A = get_float("  E: ")
    print("\nStation B:"); N_B = get_float("  N: "); E_B = get_float("  E: ")
    
    if input("Have bearings? (y/n): ").lower() == 'y':
        brg_A = get_dms("Bearing A→P")
        brg_B = get_dms("Bearing B→P")
    else:
        dN, dE = N_B - N_A, E_B - E_A
        brg_AB = math.degrees(math.atan2(dE, dN)) % 360
        brg_BA = (brg_AB + 180) % 360
        brg_A = (brg_AB - get_dms("Angle at A")) % 360
        brg_B = (brg_BA + get_dms("Angle at B")) % 360
    
    dN, dE = N_B - N_A, E_B - E_A
    sin_diff = math.sin(math.radians(brg_B - brg_A))
    if abs(sin_diff) < 1e-12:
        print("Error: Parallel"); return None, None, None
    
    dAP = (dN * math.sin(math.radians(brg_B)) - dE * math.cos(math.radians(brg_B))) / sin_diff
    brg_A_rad = math.radians(brg_A)
    N_P = N_A + dAP * math.cos(brg_A_rad)
    E_P = E_A + dAP * math.sin(brg_A_rad)
    
    print(f"\n✅ Provisional {name}: N={N_P:.3f}, E={E_P:.3f}")
    return N_P, E_P, name

def compute_cut(N_C, E_C, N_P, E_P, brg_deg):
    """Original coordinate-cut method."""
    brg_rad = math.radians(brg_deg)
    delta_N = N_P - N_C
    delta_E = E_P - E_C
    sin_b = math.sin(brg_rad)
    cos_b = math.cos(brg_rad)
    
    if abs(sin_b) < 1e-12 or abs(cos_b) < 1e-12:
        return float('inf'), float('inf')
    
    cot_b = cos_b / sin_b
    tan_b = sin_b / cos_b
    cut_N = cot_b * delta_E - delta_N
    cut_E = tan_b * delta_N - delta_E
    
    return cut_N, cut_E

def least_squares_adjust(stations, N_P, E_P):
    """Find best intersection using least squares."""
    A = []
    b = []
    
    for N_C, E_C, brg in stations:
        brg_rad = math.radians(brg)
        sin_b = math.sin(brg_rad)
        cos_b = math.cos(brg_rad)
        
        if abs(sin_b) < 1e-12 or abs(cos_b) < 1e-12:
            continue
            
        tan_b = sin_b / cos_b
        delta_N = N_P - N_C
        delta_E = E_P - E_C
        cut_E = tan_b * delta_N - delta_E
        
        A.append([-tan_b, 1])
        b.append(cut_E)
    
    if len(A) < 2:
        print("Error: Need at least 2 valid stations")
        return 0, 0
    
    A = np.array(A)
    b = np.array(b)
    
    solution = np.linalg.lstsq(A, b, rcond=None)[0]
    return solution[0], solution[1]

def plot_cuts(stations, cuts, dN_opt, dE_opt, name):
    """Plot with proper layout - no overlap between plot and table."""
    fig = plt.figure(figsize=(14, 10))
    
    # Main plot area
    ax_plot = fig.add_axes([0.1, 0.28, 0.8, 0.68])
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    limit = 1.6
    
    ax_plot.set_xlim(-limit, limit)
    ax_plot.set_ylim(-limit, limit)
    
    ticks = np.arange(-1.6, 1.8, 0.2)
    ax_plot.set_xticks(ticks)
    ax_plot.set_yticks(ticks)
    
    # Plot each station
    for i, ((N_C, E_C, brg), (cut_N, cut_E)) in enumerate(zip(stations, cuts)):
        if abs(cut_N) > 1e6:
            continue
            
        color = colors[i % len(colors)]
        
        # Plot cut point
        ax_plot.plot(cut_E, cut_N, 's', color=color, markersize=10, 
                markeredgecolor='black', markeredgewidth=1.5)
        
        # Calculate line slope from bearing
        brg_rad = math.radians(brg)
        sin_b = math.sin(brg_rad)
        cos_b = math.cos(brg_rad)
        
        if abs(sin_b) > 1e-12:
            slope = -cos_b / sin_b
            
            x_vals = np.linspace(-limit, limit, 100)
            y_vals = cut_N + slope * (x_vals - cut_E)
            ax_plot.plot(x_vals, y_vals, '--', color=color, linewidth=2, 
                   label=f"{chr(99+i)} line", alpha=0.8)
        
        # Labels
        station = chr(99+i)
        if station == 'c':
            ax_plot.annotate(f"ΔN={cut_N:.2f}\nΔE={cut_E:.2f}", 
                       (cut_E, cut_N), xytext=(-40, -40), 
                       textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                edgecolor=color, linewidth=1.5),
                       fontsize=9, ha='center')
        elif station == 'd':
            ax_plot.annotate(f"ΔN={cut_N:.2f}\nΔE={cut_E:.2f}", 
                       (cut_E, cut_N), xytext=(0, 40), 
                       textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                edgecolor=color, linewidth=1.5),
                       fontsize=9, ha='center')
        else:
            ax_plot.annotate(f"ΔN={cut_N:.2f}\nΔE={cut_E:.2f}", 
                       (cut_E, cut_N), xytext=(40, 0), 
                       textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                edgecolor=color, linewidth=1.5),
                       fontsize=9, ha='left')
    
    # Plot adjusted position
    ax_plot.plot(dE_opt, dN_opt, 's', color='purple', markersize=10,
           markeredgecolor='black', markeredgewidth=2, label='True P')
    
    # Formatting
    ax_plot.axhline(y=0, color='k', linewidth=1.5)
    ax_plot.axvline(x=0, color='k', linewidth=1.5)
    ax_plot.grid(True, linestyle='--', alpha=0.5)
    ax_plot.set_xlabel('EASTING (m)', fontsize=12, fontweight='bold')
    ax_plot.set_ylabel('NORTHING (m)', fontsize=12, fontweight='bold')
    ax_plot.set_title(f'CUTS CARTESIAN PLOT\n(ΔN vs ΔE for Each Station)', 
                fontsize=14, fontweight='bold', pad=10)
    
    ax_plot.set_aspect('equal')
    
    # North arrow
    ax_plot.annotate('N', xy=(0.55, 0.25), xycoords='axes fraction',
               fontsize=14, fontweight='bold', ha='center')
    ax_plot.arrow(0.55, 0.15, 0, 0.08, head_width=0.03, head_length=0.03, 
            fc='black', ec='black', transform=ax_plot.transAxes)
    
    ax_plot.legend(loc='upper right', fontsize=9)
    
    # Table area (separate)
    ax_table = fig.add_axes([0.1, 0.05, 0.8, 0.18])
    ax_table.axis('off')
    
    # Create table
    table_data = [[chr(99+i), f"{cut_N:.2f}", f"{cut_E:.2f}"] 
                  for i, (cut_N, cut_E) in enumerate(cuts) if abs(cut_N) < 1e6]
    
    table = ax_table.table(cellText=table_data,
                    colLabels=['STATION', 'CUT (N)', 'CUT (E)'],
                    loc='center',
                    cellLoc='center',
                    bbox=[0, 0, 1, 1])
    
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2)
    
    # Style header
    for i in range(3):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
        table[(0, i)].set_height(0.3)
    
    # Style rows
    for i in range(len(table_data)):
        for j in range(3):
            table[(i+1, j)].set_facecolor('#E7E6E6' if i % 2 == 0 else 'white')
            table[(i+1, j)].set_height(0.25)
    
    plt.show()

def main():
    print("="*60 + "\nSURVEYING CUT COMPUTATION (CORRECTED)\n" + "="*60)
    
    if input("Have provisional coordinates? (y/n): ").lower() == 'y':
        name = input("Station name: ").strip() or "P"
        N_P = get_float("Northing: ")
        E_P = get_float("Easting: ")
    else:
        N_P, E_P, name = compute_provisional()
        if N_P is None: return
    
    n = int(input("\nNumber of known stations: "))
    stations = []
    
    print("\n--- Enter station data ---")
    for i in range(n):
        print(f"\nStation {chr(99+i)}:")
        stations.append((get_float("  Northing: "), 
                        get_float("  Easting: "), 
                        get_dms("  Bearing: ")))
    
    print("\n" + "="*60 + "\nCUTS\n" + "="*60)
    cuts = []
    for i, (N_C, E_C, brg) in enumerate(stations):
        cut_N, cut_E = compute_cut(N_C, E_C, N_P, E_P, brg)
        cuts.append((cut_N, cut_E))
        if abs(cut_N) < 1e6:
            print(f"{chr(99+i)}: ΔN = {cut_N:+.2f}, ΔE = {cut_E:+.2f}")
    
    dN_opt, dE_opt = least_squares_adjust(stations, N_P, E_P)
    
    print("\n" + "="*60 + "\nADJUSTMENT\n" + "="*60)
    print(f"Provisional: N = {N_P:.4f}, E = {E_P:.4f}")
    print(f"Adjusted:    N = {N_P + dN_opt:.4f}, E = {E_P + dE_opt:.4f}")
    print(f"Corrections: ΔN = {dN_opt:+.4f}, ΔE = {dE_opt:+.4f}")
    
    plot_cuts(stations, cuts, dN_opt, dE_opt, name)
    
    print(f"\n✅ Complete: {name} = ({N_P + dN_opt:.4f}, {E_P + dE_opt:.4f})")

if __name__ == "__main__":
    main()