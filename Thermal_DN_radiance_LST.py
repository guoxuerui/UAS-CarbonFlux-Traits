# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 09:49:23 2025

@author: xu.guo
"""
#only keep 12:00 one moment but not the triple mean
import numpy as np
import rasterio
from rasterio import Affine, MemoryFile
NODATA_VAL = -10000.0  # original NoData
# --- Constants & Parameters ---
epsilon_obj = 0.98
Dist = 100  # meters

# Atmospheric parameters
TAtmC = 7.89
TAtmK = TAtmC + 273.15
TextC = 7.89
TextK = TextC + 273.15
RH = 0.4528
Ld = 323.53

# Transmission calculation
X = 1.9
alpha1 = 0.006569
alpha2 = 0.012620
beta1 = -0.002276
beta2 = -0.006670

H2O = RH * np.exp(1.5587 + 0.06939 * TAtmC - 0.00027816 * TAtmC**2 + 0.00000068455 * TAtmC**3)
Tau = X * np.exp(-np.sqrt(Dist) * (alpha1 + beta1 * np.sqrt(H2O))) + (1 - X) * np.exp(-np.sqrt(Dist) * (alpha2 + beta2 * np.sqrt(H2O)))

#H2O_ext = RH * np.exp(1.5587 + 0.06939 * TextC - 0.00027816 * TextC**2 + 0.00000068455 * TextC**3)
Tau_ext = 1
# Camera calibration
R1 = 17096.453
R2 = 0.045601364
R = R1 / R2
O = -235
B = 1428
F = 1

# Radiance and emissivity
sigma = 5.670374419e-8
def compute_radiance(T_kelvin):
    return R / (np.exp(B / T_kelvin) - F) - O

def calculate_sky_emissivity(Ld, Ta, sigma=5.670374419e-8):
    return Ld / (sigma * Ta**4)

epsilon_sky = calculate_sky_emissivity(Ld, TAtmK)
rad_sky = compute_radiance(TAtmK)
rad_air = rad_sky
rad_ext = compute_radiance(TextK)

def calculate_rad_obj(rad_tot, tau, epsilon_obj, tau_ext, epsilon_sky, rad_sky, rad_air, rad_ext=None):
    if rad_ext is None:
        rad_ext = rad_air
    term1 = (1 / (tau * epsilon_obj )) * rad_tot
    term2 = ((1 - epsilon_obj) / epsilon_obj) * epsilon_sky * rad_sky
    term3 = ((1 - tau) / (tau * epsilon_obj)) * rad_air
    return term1 - term2 - term3

def VueProR_rad2temp(B, R, O, F, rad_obj):       
    return (B / np.log(R / (rad_obj + O) + F)) - 273.15

def planck_radiance(T_K, lambda_um):
    h = 6.626e-34
    c = 2.998e8
    kB = 1.381e-23
    lambda_m = lambda_um * 1e-6
    numerator = 2 * h * c**2 / lambda_m**5
    exponent = (h * c) / (lambda_m * kB * T_K)
    denominator = np.exp(exponent) - 1
    return (numerator / denominator) * 1e-6


# --- Apply to a raster ---
input_path = r"C:\Users\xu.guo\OneDrive - Forschungszentrum Jülich GmbH\Desktop\AgroC_SCOPE project\2023_processing\all-dates_FLIR-Radiance\radiance_temperature\reprojection\radiance_0328_transparent_reflectance_grayscale_modified.tif"
output_temp_path = r"C:\Users\xu.guo\OneDrive - Forschungszentrum Jülich GmbH\Desktop\AgroC_SCOPE project\2023_processing\all-dates_FLIR-Radiance\20230328_radiance_temperature_C.tif"
output_radiance_path = r"C:\Users\xu.guo\OneDrive - Forschungszentrum Jülich GmbH\Desktop\AgroC_SCOPE project\2023_processing\all-dates_FLIR-Radiance\20230328_radiance_10um.tif"

with rasterio.open(input_path) as src:
    dn_array = src.read(1).astype(np.float64)
    profile = src.profile.copy()
    srcmask = (src.read_masks(1) == 0) | (dn_array == NODATA_VAL)
    dn_array = np.where(srcmask, np.nan, dn_array)
    
    # Calculate object radiance
    rad_obj = calculate_rad_obj(
        dn_array, Tau, epsilon_obj, Tau_ext,
        epsilon_sky, rad_sky, rad_air, rad_ext
    )

    # Convert to temperature
    temperature_C = VueProR_rad2temp(B, R, O, F, rad_obj)
    
    # Convert to spectral radiance at 10 µm
    radiance_10um = planck_radiance(temperature_C + 273.15, lambda_um=10.0)

# --- Save results ---
profile.update(driver="GTiff",dtype=rasterio.float32)

with rasterio.open(output_temp_path, "w", **profile) as dst:
    dst.write(temperature_C.astype(np.float32), 1)

with rasterio.open(output_radiance_path, "w", **profile) as dst:
    dst.write(radiance_10um.astype(np.float32), 1)

print("✅ Temperature and radiance rasters written.")
