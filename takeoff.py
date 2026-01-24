import math
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
W = 6.166667                                                    #headwind
g = 9.81
rho = 1.225 * (1-0.02256*0.173) **4.256
print("Density:",rho, "kg/m^3")
k = 1.4
R = 287
c_p = 1005
P_o = 101.325 *1000*(1-0.02256 *0.173) **5.256                 #Pastandard
print("Pressure:",P_o,"Pa")
T_o = 15+ 273.15                                               #Kstandard
S_nozzle = 1.64217                                             #m^2,nozzle area
F_To = 4.4482216153*57160                                      #N,staticthrust
S = .16                                                        #m^2,wing area
V_R = 0.51444444444*171                                        #m/s,rotationspeed
MTOW = 0.453592 *875000                                        #kg
PR = ((R* F_To)/(2*c_p* P_o*S_nozzle)+1)**(k/(k-1))
mu = 0.03                                                      #frictioncoefficient
print("FanPressure Ratio:",PR)

#----------thrust model---------

def calc_thrust(T_o, V,c_p,P_o,k,PR,R,S_nozzle,W):
    Va = V + W
    T_t0 = T_o +(Va**2)/(2* c_p)
    P_t0 = P_o *(T_t0/T_o)**(k /(k-1))
    P_t18 = PR *P_t0

    T_t18 = T_t0 *PR**((k-1) /k)
    T_18 = T_t18 *(P_o/P_t18)**((k-1) /k)

    V_18 = math.sqrt(2*c_p*(T_t18-T_18))
    F_T = (P_o /(R*T_18))*S_nozzle*V_18*(V_18-Va)
    return F_T

def get_C_Lg(MTOW, g,rho,S,V_R):
 return 2* (MTOW*g)/(rho* (V_R**2)*S)

print("C_Lg:",get_C_Lg(MTOW,g,rho,S,V_R))

#----------ground run integrand---------
def f(V):
 T = 4* calc_thrust(T_o,V,c_p,P_o, k, PR,R,S_nozzle,W)
 C_Lg = get_C_Lg(MTOW,g,rho,S,V_R)
 L = 1.2 *0.5*rho*S* C_Lg*(V+W) ** 2
 C_Dg = 0.076 +0.0435*C_Lg**2
 D = 0.5 *rho*S*C_Dg* (V+W)**2
 if (mu* (MTOW*g-L))>=0:
     a = (T-D-mu*(MTOW*g-L)
         -MTOW*g *math.sin(0.0009999996667))/ MTOW
 else:

     a = (T-D-MTOW*g*math.sin(0.0009999996667))/MTOW
 return 2* V/a

#----------numerical integration---------

x,error = quad(f,0,V_R-W)
x1 = 0.5* x

print("Ground speed atrotation:",V_R-W, "m/s")
print("Airspeedat rotation:",V_R,"m/s")
print(f"Ground-rundistancex1={x1:.2f}m(±{error:.2e})")

def time_integrand(V):
    T = 4* calc_thrust(T_o,V,c_p,P_o, k, PR,R,S_nozzle,W)
    C_Lg = get_C_Lg(MTOW,g,rho,S,V_R)
    L = 1.2 *0.5*rho*S* C_Lg*(V+W) ** 2
    C_Dg = 0.076 +0.0435*C_Lg**2
    D = 0.5 *rho*S*C_Dg* (V+W)**2
    if(mu* (MTOW*g-L))>=0:
      a = (T-D-mu*(MTOW*g-L)
     -MTOW*g *math.sin(0.0009999996667))/ MTOW
    else:
     a = (T-D-MTOW*g*math.sin(0.0009999996667))/MTOW
     return 1/ a

takeoff_time_sec,_ =quad(time_integrand,0,V_R-W)
print(f"Take-offroll time:{takeoff_time_sec:.2f}s")

avg_thrust_lbf = (4*calc_thrust(T_o,V_R,c_p,P_o,k,PR,
R,S_nozzle,W) /4.4482216153)
SFC = 0.316                                                      #lb/(lbf·hr)
fuel_burn_lb=SFC*avg_thrust_lbf*(takeoff_time_sec/3600)
print(f"Fuelburned (groundrun):{fuel_burn_lb:.2f} lb")

print("Aircraftweightatrotation:",
MTOW-fuel_burn_lb *0.453592,"kg")

#----------rotation to screen height---------

T_unstick = 4 *calc_thrust(T_o,V_R,c_p, P_o,k,PR,R,S_nozzle, W)
C_Lg = get_C_Lg(MTOW,g,rho,S,V_R)
L = 1.2* 0.5 *rho*S* C_Lg*V_R**2
C_Dg = 0.076 +0.0435*C_Lg**2
D_unstick = 0.5 *rho*S*C_Dg *V_R**2
theta = math.asin((T_unstick-D_unstick)/ (MTOW*g))
x2 = 15.24 /math.tan(theta)
print("Pitchangle:", math.degrees(theta),"deg")
print("Climbdistance to50ft:",x2,"m")
print("Totaltake-off distance:",x1+x2,"m")

climb_time = x2 /V_R
fuel_burn_lb2 = SFC*avg_thrust_lbf*((takeoff_time_sec+climb_time)/3600)
print("Fuelburned (fullTO):",fuel_burn_lb2,"lb")
print("Weightat screenheight:",
MTOW-fuel_burn_lb2*0.453592,"kg")
print("Staticthrust perengine:",
calc_thrust(T_o,V_R, c_p,P_o,k,PR,R,S_nozzle,0),"N")

#----------thrust vs.ground speedp lot---------

ground_speeds = np.linspace(1,V_R-W,200)
thrusts_per_engine = [calc_thrust(T_o,V,c_p,P_o, k,PR,
R,S_nozzle,W)
for V in ground_speeds]

plt.figure()
plt.plot(ground_speeds, thrusts_per_engine)
plt.xlabel("GroundSpeed(m/s)")
plt.ylabel("ThrustperEngine(N)")
plt.title("ThrustEvolutionperEnginevs.Ground Speed")
plt.grid(True)
plt.tight_layout()
plt.show()

#---------- speed vs. distance during ground run---------

V_vals = np.linspace(0.1, V_R- W, 300)
s_vals, t_vals = [0], [0]
for i in range(1, len(V_vals)):
 V1, V2 = V_vals[i- 1], V_vals[i]
 dV, V_mid = V2- V1, 0.5 * (V1 + V2)
 T = 4 * calc_thrust(T_o, V_mid, c_p, P_o, k, PR,
 R, S_nozzle, W)
 C_Lg = get_C_Lg(MTOW, g, rho, S, V_R)
 L = 1.2 * 0.5 * rho * S * C_Lg * (V_mid + W) ** 2
 C_Dg = 0.076 + 0.0435 * C_Lg ** 2
 D = 0.5 * rho * S * C_Dg * (V_mid + W) ** 2
if (mu * (MTOW * g- L)) >= 0:
 a = (T- D- mu * (MTOW * g- L)- MTOW * g * math.sin(0.0009999996667)) / MTOW
else:
 a = (T- D- MTOW * g * math.sin(0.0009999996667)) / MTOW
dt = dV / a
ds = V_mid * dt
t_vals.append(t_vals[-1] + dt)
s_vals.append(s_vals[-1] + ds)
airspeeds = V_vals + W
plt.figure()
plt.plot(s_vals, V_vals, label="Ground Speed")
plt.plot(s_vals, airspeeds, "--", label="Airspeed (Ground + Wind)")
plt.xlabel("Distance (m)")
plt.ylabel("Speed (m/s)")
plt.title("Ground Speed and Airspeed vs. Distance During Take-Off Roll")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()