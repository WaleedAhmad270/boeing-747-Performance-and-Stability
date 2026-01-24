# Boeing 747-400 Flight Performance and Stability Analysis

This repository contains numerical models and analysis scripts developed for a complete flight-mission performance and longitudinal stability study of the **Boeing 747-400**.

The project models all major phases of flight using first-principles aerospace engineering methods and validated atmospheric, aerodynamic, and propulsion relations.

---

## Project Overview

The analysis covers a realistic long-haul mission (Lufthansa Flight LH471: Toronto → Frankfurt) and includes:

- Takeoff performance with variable thrust and ground-run integration  
- Segmented climb performance with fuel burn tracking  
- Cruise performance at multiple flight levels  
- Descent performance under idle thrust conditions  
- Landing performance including flare, braking, runway slope, and friction effects  
- Longitudinal static stability analysis (neutral point and maneuver point)

All calculations are based on publicly available aircraft data for the **Boeing 747-400 equipped with GE CF6-80C2 engines**, standard aerodynamic theory, and ISA atmospheric models.

---

## Flight Phases Modeled

### ✈️ Takeoff
- Variable thrust modeling using a constant compression approach  
- Ground run solved via numerical integration  
- Headwind, runway slope, rolling friction, and ground effect included  

### ⛰️ Climb
- Segmented climb profile:
  - Sea level → 10,000 ft (250 KIAS)
  - 10,000 ft → 25,000 ft (340 KIAS)
  - 25,000 ft → cruise altitude (Mach 0.84)
- Fuel burn and mass updated continuously

### 🌍 Cruise
- Steady, level flight at constant Mach number  
- Time-marching integration of fuel burn  
- Step-climb modeled between cruise flight levels

### ⬇️ Descent
- Idle thrust descent modeling  
- Two-segment profile (Mach → KIAS)  
- Constant vertical rate consistent with operational procedures

### 🛬 Landing
- Approach speed from stall margin  
- Numerical integration of landing roll  
- Sensitivity analysis for:
  - Runway friction coefficient
  - Runway slope

---

## Stability Analysis

The repository also includes a **longitudinal static stability analysis**, featuring:

- Mean Aerodynamic Chord (MAC) calculation  
- Wing and tail aerodynamic center estimation  
- Neutral point determination  
- Maneuver point variation with altitude and aircraft mass  

This provides insight into how stability margins evolve throughout the mission.

---


