# Unsicherheitsrechnung
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Daten
# allgemeine Parametern
rho = 1000.0          # Dichte (kg/m³)
g = 9.81              # Erdbeschleunigung (m/s²)
t = 20.0              # Zeit (s)

# Kapillare1Daten(duenn)
h1 = np.array([0.536,0.765,1.017,1.318,1.738])  # Hoehenunterschied (m)
V1 = np.array([0.6376e-6, 0.9012e-6, 1.1699e-6, 1.4882e-6, 1.9376e-6])  # Volumen (m³)
l1 = 0.044               # Laenge (m)
r1 = 0.00019             # Radius (m)

# Kapillare2Daten(dick)
h2 = np.array([1.458,1.16,0.995,0.794,0.504])  # Hoehenunterschied (m)
V2 = np.array([3.3327e-6, 2.6321e-6, 2.3424e-6, 1.8606e-6, 1.1868e-6])  # Volumen (m³)
l2 = 0.031               # Laenge (m)
r2 = 0.000215            # Radius (m)

# Unsicherheiten beim Messen
delta_rho = 50.0      # (kg/m³)
delta_h = 0.0005      # (m)
delta_V = 1e-7        # (m³)
delta_t = 0.083       # (s)
delta_l = 0.0001      # (m)
delta_r = 5e-6        # (m)

# Funktion
def calculate_viscosity(h, V, t, l, r, label):
    
    # Berechnung der Stromstaerke i = V/t und der Unsicherheit
    i = V / t
    delta_i = i * np.sqrt((delta_V/V)**2 + (delta_t/t)**2)
    
    # Berechnung der Druckdifferenz Δp = ρgh und der Unsicherheit
    delta_p = rho * g * h
    delta_delta_p = delta_p * np.sqrt((delta_rho/rho)**2 + (delta_h/h)**2)
    
    # lineare Regression W = Δp/i
    slope, intercept, r_value, _, std_err = linregress(i, delta_p)
    W = slope
    delta_W = np.sqrt(std_err**2 + (0.05 * W)**2)  # Regressionsfehler + 5% systematischer Fehler
    
    # Berechnung der Viskositaet η und der Unsicherheit
    eta = (W * np.pi * r**4) / (8 * l)
    delta_eta = eta * np.sqrt(
        (delta_W/W)**2 + 
        (4 * delta_r/r)**2 + 
        (delta_l/l)**2
    )
    
    # Rueckgabe Ergebnisse
    return {
        'label': label,
        'i': i,
        'delta_i': delta_i,
        'delta_p': delta_p,
        'delta_delta_p': delta_delta_p,
        'W': W,
        'delta_W': delta_W,
        'eta': eta,
        'delta_eta': delta_eta,
        'r_value': r_value
    }

# Daten berechnen
result1 = calculate_viscosity(h1, V1, t, l1, r1, "Kapillare1 (duenn)")
result2 = calculate_viscosity(h2, V2, t, l2, r2, "Kapillare2 (dick)")

# Visualisierung
plt.figure(figsize=(14, 5))

# Graph1: Δp-i
plt.subplot(1, 2, 1)
for result in [result1, result2]:
    plt.errorbar(
        result['i'], result['delta_p'], 
        xerr=result['delta_i'], yerr=result['delta_delta_p'], 
        fmt='o', label=f"{result['label']} (R²={result['r_value']**2:.3f})"
    )
    plt.plot(
        result['i'], result['W'] * result['i'] + 0, 
        '--', label=f"W={result['W']:.2e}±{result['delta_W']:.2e} Pa·s/m³"
    )
plt.xlabel('Stromstaerke i (m³/s)')
plt.ylabel('Druckdifferenz Δp (Pa)')
plt.title('Vergleich der Druckdifferenz-Stromstaerke-Beziehungen für zwei Kapillaren')
plt.legend()
plt.grid()

# Graph2: Viskositaet
plt.subplot(1, 2, 2)
for i, result in enumerate([result1, result2]):
    plt.errorbar(
        i, result['eta'], yerr=result['delta_eta'], 
        fmt='o', capsize=5, label=result['label']
    )
plt.xticks([0, 1], [result1['label'], result2['label']])
plt.ylabel('Viskositaet η (Pa·s)')
plt.title('Vergleich der Ergebnisse von Viskositätsmessungen')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()

# Ausgabe der Ergebnisse
print("\nErgebnisse:")
for result in [result1, result2]:
    print(f"\n{result['label']}:")
    print(f"Stroemungswiderstand W = （{result['W']:.3e} ± {result['delta_W']:.3e}）Pa·s/m³")
    print(f"Viskositaet η = （{result['eta']:.3e} ± {result['delta_eta']:.3e}）Pa·s")
    print(f"Unsicherheiten = {result['delta_eta']/result['eta']*100:.1f}%")
    print(f"Bestimmtheitsmaß R² = {result['r_value']**2:.4f}")
