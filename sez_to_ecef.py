# sez_to_ecef.py
#
# Usage: python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km
# Converts SEZ frame parameters with respect to the origin to ECEF components
# which consists of two matrix rotations and a translation from SEZ vector
# Parameters:
#  o_lat_deg: SEZ origin latitude in degrees
#  o_long_deg: SEZ origin longitude in degrees
#  o_hae_km: SEZ origin height above ellipsoid in km
#  s_km: SEZ s-component (South) in km
#  e_km: SEZ e-component (East) in km
#  z_km: SEZ z-component in km
# Output:
#  ecef_x_km: ECEF x-component in km after two rotations and a translation
#  ecef_y_km: ECEF y-component in km after two rotations and a translation
#  ecef_z_km: ECEF z-component in km after two rotations and a translation
#
# Written by Mandar Ajgaonkar
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
import math  # math module
import sys   # argv
import numpy as np # arrays, matrices, and operations between them

# "constants"
R_E_KM = 6378.1363 # radius of Earth in km
E_E = 0.081819221456 # unitless

# helper functions

# function description

## calculate denominator 
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0 - (ecc**2.0) * (math.sin(lat_rad))**2)

# initialize script arguments
o_lat_deg = float('nan')  
o_lon_deg = float('nan')  
o_hae_km = float('nan')  
s_km = float('nan')
e_km = float('nan')
z_km = float('nan')

# parse script arguments 
if len(sys.argv) == 7:
    try:
        o_lat_deg = float(sys.argv[1])
        o_long_deg = float(sys.argv[2])
        o_hae_km = float(sys.argv[3])
        s_km = float(sys.argv[4])
        e_km = float(sys.argv[5])
        z_km = float(sys.argv[6])
    except ValueError:
        print("Error: lat_deg, long_deg and hae_km must be numeric.")
        exit()
else:
    print('python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km')
    exit()

# convert longitude and latitude to radians
o_lat_rad = o_lat_deg * math.pi / 180.0
o_long_rad = o_long_deg * math.pi / 180.0

# perform calculations to get components of r in ECEF reference frame
denominator = calc_denom(E_E, o_lat_rad)
c_E = R_E_KM / denominator
s_E = R_E_KM * (1.0 - E_E * E_E) / denominator
r_x_km = (c_E + o_hae_km) * math.cos(o_lat_rad) * math.cos(o_long_rad)
r_y_km = (c_E + o_hae_km) * math.cos(o_lat_rad) * math.sin(o_long_rad)
r_z_km = (s_E + o_hae_km) * math.sin(o_lat_rad)
r = np.array([r_x_km, r_y_km, r_z_km])

# transpose r vector to a column vector 
r = r.reshape(-1,1)

# define SEZ vector 
r_SEZ = np.array([s_km, e_km, z_km])

# transpose the SEZ vector to a column vector
r_SEZ = r_SEZ.reshape(-1, 1)

# first matrix rotation using latitude
Ry_lat = np.array([[math.sin(o_lat_rad), 0, math.cos(o_lat_rad)], 
      [0, 1, 0], [-math.cos(o_lat_rad), 0, math.sin(o_lat_rad)]]) 
Rz_long = np.array([[math.cos(o_long_rad), -math.sin(o_long_rad), 0], 
           [math.sin(o_long_rad), math.cos(o_long_rad), 0], [0, 0, 1]])

# calculate result of first rotation
rot_1 = np.dot(Ry_lat, r_SEZ)

# calculate result of second rotation
rot_2 = np.dot(Rz_long, rot_1)

# final results with a translation from ECEF vector calculated from LLH
ecef_x_km = r[0] + rot_2[0]
ecef_y_km = r[1] + rot_2[1]
ecef_z_km = r[2] + rot_2[2]

# reformat as a float rather than an array of float
ecef_x_km = ecef_x_km[0]
ecef_y_km = ecef_y_km[0]
ecef_z_km = ecef_z_km[0]

# print final results 
print(ecef_x_km)
print(ecef_y_km)
print(ecef_z_km)