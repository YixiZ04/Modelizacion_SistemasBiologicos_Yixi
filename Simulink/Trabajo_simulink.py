#################################################################################################################
# Author: Yixi Zhang
# Version: 1.0.0
# Description: A program to simulate a bioreactor: bacterial growth, substrate usage and temperature variation.
# Usage: Run the entire program with a Python IDE or on the terminal.
#################################################################################################################


## Importing the needed libraries

import matplotlib.pyplot as plt
import numpy as np
import time
import os
import re

## 1. Environmental temperature variation

def env_temp_variation (sim_time, A, bias, frequency):
    """
    Input: time, A, bias, frequency
    Output: sine wave
    """
    return A * np.sin (frequency * sim_time) + bias

## 2. Batch culture. Supposing exceeding substrate (always at the max growth rate)

def bacterial_growth (sim_time,X0, mu):
    """
    Input: X0 (DO600, g, g/l...); mu (1/t), time (t) where t == any time unit;
    Output: A list with bacteria's population at time t.
    """
    return X0*np.exp(mu*sim_time)

def substrate_usage (biomass, Yxs):
    """
    Input: biomass (g), Yxs (substrate usage ratio, g biomass produced/ substrate used).
    Output: Usage of substrate.
    """
    return biomass / Yxs

## 3.Temperature variation due to bacterial growth.
def celcius_elevation_per_hour (heat_production_cal, cp, batch_vol):
    """
    Input: heat production (cal/g biomass/s), Cp (specific heat) in J/kg H2O/ºC,
    batch_vol (L)
    Output: celcius/g biomass / h
    """
    cal_to_jules = 4.184 # J/cal
    heat_production_jules = heat_production_cal * cal_to_jules # J/g biomass/s
    heat_production_per_hour = heat_production_jules * 3600 # J/g biomass / h
    # Suppose 150L Batch == 150 Kg H2O.
    specific_heat_batch = cp * batch_vol # J/ºC
    celcius_per_biomass_per_hour = heat_production_per_hour / specific_heat_batch #ºC/g_biomass/h
    return round (celcius_per_biomass_per_hour,6)

def temp_var_ratio_due2_bacterial_growth (biomass, ta_elevation_per_hour):
    """
    Input: biomass (g). celcius elevation per hour (ºC/g biomass/h)
    Output: temperature variation per hour (ºC/h)
    """
    return biomass * ta_elevation_per_hour


def temp_var_due2_bacterial_growth (sim_time, temp_var_ratio):
    """
    Input: the ratio for temperature variation due to bacterial growth (ºC/h). Time (h)
    Output: temperature variation due to bacterial growth.
    """
    step = sim_time [1] - sim_time [0]
    return temp_var_ratio * step
## 4. Effect of environmental temperature to medium.
def compute_cylinder_area (base_diameter, height):
    """
    Input: height, base_diameter of a cylinder. In any unit, but both unit has to be the same.
    Output: area of the cylinder in the given unit square.
    """
    base_radius = base_diameter / 2
    bases_area = 2 * (np.pi*base_radius)
    wall_area = 2*np.pi*base_radius *height
    total_area = bases_area + wall_area
    return total_area

def heat_loss_bioreactor (temp_0, area_bioreactor, thermal_conductivity, wall_thickness, temp_env, sim_time, mass, cp):
    """
    Input: initial temperature of bioreactor, the area of the bioreactor, the thermal conductivity, wall thickness, environmental temperature variation, simulation time, mass, cp
    Output: temperature variation of the bioreactor with thermostat off
    """
    temp_bioreactor = [temp_0] # List to store temperature variation.
    # Defining dt
    step = sim_time [1] - sim_time [0]
    dt = step/1000
    #Calculate the resistence
    resistence = wall_thickness / (area_bioreactor * thermal_conductivity)
    # The temperature variation for every step.
    for i in range (len (temp_env) - 1):
        temp_diff = temp_bioreactor [i] - temp_env [i] #ºC or K
        dqdt = temp_diff / resistence
        temp_change_ratio = dqdt / (mass * cp) # ºC/s
        temp_change_ratio_per_hour = temp_change_ratio * 3600
        temp_var_per_step = temp_change_ratio_per_hour * dt
        new_temp = temp_bioreactor[i] - temp_var_per_step
        temp_bioreactor.append(new_temp)
    return temp_bioreactor

def temp_no_thermo_with_growth (temp_no_growth, temp_elevated):
    """
    Input: temperature variation with no thermo and no bacterial growth; and temperature variation due2 bacterial growth
    Output: temperature variation with no thermo and with bacterial growth.
    """
    temp_var = []
    for i in range (len (temp_no_growth)):
        real_temp = temp_no_growth [i] + temp_elevated [i]
        temp_var.append(real_temp)
    return temp_var

## 5. Thermostat, warming only.
def warming_thermostat (temp_diff, wall_thick, area_thermostat,conduct_bioreactor, mass, cp, sim_time):
    """
    Input: temperature difference between bioreactor and thermostat; the thickness of the bioreactor,
    contact area bioreactor-thermostat; conductivity bioreactor; mass and specific heat for the medium; and time of simulation
    Output: the temperature augmentation due to thermostat, and Watts consumed.
    """
    step = sim_time [1] - sim_time [0]
    dt = step / 1000 # (h)
    req_thermostat = wall_thick / (area_thermostat * conduct_bioreactor) #Thermal resistence
    dQdt = temp_diff / req_thermostat #Watts (W)
    dTdt = dQdt / (mass * cp) # ºC/s
    dTdt_hour = dTdt * 3600 #ºC/h
    temp_var = dTdt_hour * dt # ºC
    return temp_var,dQdt

def cooling_without_thermostat (temp_diff, wall_thick, area_bioreactor, conduct_bioreactor, mass, cp, sim_time):
    """
       Input: temperature difference between bioreactor and the environment; the thickness of the bioreactor,
       contact area bioreactor-thermostat; conductivity bioreactor; mass and specific heat for the medium; and time of simulation
       Output: cooling due to thermostat
    """
    dt = (sim_time [1] - sim_time [0]) / 1000 # (h)
    Req_bioreactor = wall_thick / (area_bioreactor * conduct_bioreactor) # Thermal resistence (K/W)
    dQdt = temp_diff / Req_bioreactor # (W)
    dTdt = dQdt / (mass * cp) # ºC/s
    dTdt_hour = dTdt * 3600 # ºC/h
    temp_var = dTdt_hour * dt #ºC
    return temp_var

def thermostat (temp0, temp_thermostat, temp_env, temp_bac, temp_opt, wall_thick, area_thermostat, area_bioreactor,conduct_bioreactor, mass, cp, sim_time):
    """
    Inputs: the initial temperature of the bioreactor; temperature set for the thermostat, environmental temperature variation,
    temperature variation due to bacterial growth, the thickness of the wall, contact area bioreactor-thermostat, contact area
    bioreactor-environment, conductivity of bioreactor, mass and specific heat of the medium; and time of simulation.
    Output: Temperature variation of the bioreactor. And Watts consumed per time.
    """
    # Defining variable
    step = sim_time [1] - sim_time [0]
    scaling_factor = 0.001 / step
    ddt = 0.015 /scaling_factor
    temp_bioreactor = [temp0]
    count_32 = 0 # Increases every time the temperature >= 32º. This variable makes sure working if temp0 > 32ºC
    dQdt_total = [0] #To store the Watts consumed per time.
    for i in range (len (sim_time) - 1):
        # keep_working boolean: True if tª < temp_opt - 2; else == False
        if temp_bioreactor [i] >= temp_opt + 2:
            count_32 += 1
            keep_working = False
        elif temp_bioreactor[i] <= temp_opt - 2:
            keep_working = True
        real_temp = temp_bioreactor[i] + temp_bac [i]
        # If the thermostat is on and temp0 < temp_opt + 2; it starts to work.
        if count_32 == 0 :
            temp_diff_thermostat = temp_thermostat - real_temp
            temp_diff_bioreactor = temp_env [i] - real_temp
            temp_var, dQdt = warming_thermostat(temp_diff_thermostat, wall_thick, area_thermostat, conduct_bioreactor,mass, cp, sim_time)
            temp_loss = cooling_without_thermostat(temp_diff_bioreactor, wall_thick, area_bioreactor, conduct_bioreactor, mass, cp, sim_time)
            real_temp_var = temp_var + temp_loss
            real_temp = real_temp + real_temp_var
            temp_bioreactor.append(real_temp)
            dQdt_total.append(dQdt*ddt)
        # If thermostat has worked before and the tª > temp_opt - 2. It stops working.
        elif count_32 != 0 and not keep_working:
            temp_diff_env = temp_env [i] - real_temp
            temp_var = cooling_without_thermostat(temp_diff_env, wall_thick, area_bioreactor, conduct_bioreactor, mass, cp, sim_time)
            temp_var = temp_var + real_temp
            temp_bioreactor.append(temp_var)
            dQdt_total.append(0)
        # If thermostat has worked before and the tª < temp_opt - 2. It restarts working.
        elif count_32 != 0 and keep_working:
            temp_diff_thermostat = temp_thermostat - real_temp
            temp_diff_bioreactor = temp_env[i] - real_temp
            temp_var, dQdt = warming_thermostat(temp_diff_thermostat, wall_thick, area_thermostat, conduct_bioreactor,mass, cp, sim_time)
            temp_loss = cooling_without_thermostat(temp_diff_bioreactor, wall_thick, area_bioreactor,conduct_bioreactor, mass, cp, sim_time)
            real_temp_var = temp_var + temp_loss
            real_temp = real_temp + real_temp_var
            temp_bioreactor.append(real_temp)
            dQdt_total.append(dQdt * ddt)
        # The else condition is actually redundant, since all possible cases are already covered.
        else:
            continue
    return temp_bioreactor, dQdt_total

## 6. Costs and energy consumed

def thermostat_cost (jules, euro_per_Jule):
    """
    Input: List of Jules used in real time, and euro/Jule (kWh/3.6e6)
    Output: total cost for warming, and real time cost.
    """
    cost_per_time = []
    total_cost = 0
    for jule in jules:
        cost = jule * euro_per_Jule
        total_cost += cost
        cost_per_time.append(total_cost)
    return round (total_cost,2), cost_per_time
def kWh_consumed (cost, euro_per_kWh):
    """
    Input: total cost for warming, €/kWh
    Output: total kWh consumed.
    """
    return round (cost / euro_per_kWh, 2) #kWh

def parse_csv_file (filename):
    """
    Input: a csv file with 3 items for each row in this exact format: name, euro/kg,grams_used
    Output: list of list with all information from the csv file.
    """
    file = open (filename).readlines ()
    good_file = []
    for line in file:
        line = line.rstrip ('\n')
        items = line.split(',')
        if len (items) != 3:
            continue
        else:
            good_file.append(items)
    return good_file

def medium_cost (filename):
    medium_composition = parse_csv_file(filename)
    medium_cost = 0
    medium_cost_per_component = {}
    for items in medium_composition:
        euro_per_kg = np.float16(items [1])
        gram_used = np.float16(items [2])
        cost = euro_per_kg * gram_used / 1000
        medium_cost += cost
        medium_cost_per_component[items[0]] = cost
    return round (medium_cost,2), medium_cost_per_component

def autoclave_cost (power, autoclaving_time, euro_per_kWh):
    """
    Input: power of the autoclave (kW), time autoclaving (min)
    Output: cost in euros and energy consumed (kWh)
    """
    time_in_hours = autoclaving_time / 60 # (h)
    power_used = power * time_in_hours #kWh
    cost = power_used * euro_per_kWh
    return round (cost,2), round (power_used, 2)

def air_injection_cost (pump_power, injecting_time, euro_per_kWh):
    """
    Input: power of the pump (W), time for injection (h), €/kWh
    Output: cost in euros and energy consumed (kWh)
    """
    power_used = pump_power/1000 * injecting_time
    cost = power_used * euro_per_kWh
    return round (cost,2), power_used
def calculate_impeller_power (Np, rho, rps, diameter):
    """
    Input: Number of power, density (rho): g/L, revolutions per second, diameter of the impeller (m)
    Output: impeller power (W)
    """
    return round (Np*rho*(rps**3)*(diameter**5), 2)

def stirring_cost (impeller_power, stirring_time, euro_per_kWh):
    """
    Input: impeller power (W), stirring time (h), €/kWh
    Output: cost in euros and energy consumed (kWh)
    """
    power_used = impeller_power/1000 * stirring_time
    cost = power_used * euro_per_kWh
    return round (cost, 2), round (power_used, 2)

def total_cost (thermostat_cost, medium_cost, autoclave_cost, air_injection_cost, stirring_cost):
    """
    Input: thermostat_cost, medium_cost, autoclave_cost, stirring cost
    Output: total cost in real time.
    """
    other_costs = np.float16(medium_cost) + np.float16 (autoclave_cost) + np.float16 (air_injection_cost) + np.float16 (stirring_cost)
    total_cost_in_time = []
    for cost in thermostat_cost:
        cost_in_total = cost + other_costs
        total_cost_in_time.append (cost_in_total)
    return total_cost_in_time

## 7. Plots
def save_plots (fig_name, dir_name):
    actual_dir = '.'
    path2dir = actual_dir + '/'+ dir_name
    os.makedirs (path2dir, exist_ok=True)
    image_id = path2dir + '/' + fig_name + '.png'
    print (f'The plots will be saved in {path2dir}')
    plt.savefig (image_id, format = 'png')


def plot_bacterial_growth_substrate_usage (sim_time, biomass, substrate_consumed, S0):
    """
    Input: time (h), biomass, substrate_consumed, S0 (initial substrate)
    Output: 2 plots: one for bacterial growth and one for substrate usage and resting substrate.
    """
    rest = S0 - substrate_consumed
    fig, axis = plt.subplots (1,2, figsize = (12, 5))
    axis[0].plot (sim_time, biomass, label = "Bacterial Growth (g)")
    axis[0].set_title ("Exponential growth")
    axis[0].set_xlabel ("Time (h)")
    axis[0].set_ylabel ("Biomass(g)")
    axis[0].set_xlim ([0, 168])
    axis[0].grid (True)
    axis[0].legend ()
    axis[1].plot (sim_time, substrate_consumed, label = "Substrate Usage (g)")
    axis[1].plot (sim_time, rest, label = "Restant Substrate (g)")
    axis[1].set_title ("Substrate variation")
    axis[1].set_xlabel ("Time (h)")
    axis[1].set_ylabel ("Substrate (g)")
    axis[1].set_xlim ([0, 168])
    axis[1].grid (True)
    axis[1].legend ()

def plot_temp_due2_bacterial_growth (sim_time, ta_elevation):
    """
    Input: biomass (g). celcius elevation per hour (ºC/h), time(h)
    Output: graph for temperature variation ratio.
    """
    plt.plot (sim_time, ta_elevation, 'r--', label = "Tª variation due to bacterial growth (ºC)")
    plt.title ("Temperature variation due to bacterial growth")
    plt.xlabel ("Time (h)")
    plt.ylabel ("ºC")
    plt.xlim ([0, 168])
    plt.grid (True)
    plt.legend ()

def plot_temperature_variation (sim_time, env_temp, temp_no_thermo, real_temp):
    """
    Input: simulation time, environmental temperature variation, bioreactor temperature variation with thermostat off and on.
    Output: a single plot for temperature variation
    """
    plt.figure(figsize=(10, 8))
    plt.plot (sim_time, env_temp, 'b-', label = "Environmental temperature variation")
    plt.plot (sim_time, temp_no_thermo, 'r-', label = "Tempearature with thermo off")
    plt.plot (sim_time, real_temp, 'y-', label ="Temperature with thermostat on")
    plt.title ("Temperature plots")
    plt.xlabel ("Time (h)")
    plt.ylabel ("Temperature (ºC)")
    plt.xlim([0, 168])
    plt.grid (True)
    plt.legend ()

def plot_total_costs (sim_time, thermostat_cost, total_cost_in_time):
    """
    Input: simulation time, the cost of the thermostat in real time and total cost in time.
    Output: 2 plots: 1 for thermostat only (we can better appreciate the form), and other for both thermostat cost and total cost in time.
    """
    fig, axis = plt.subplots (1,2, figsize = (12, 5))
    # Thermostat only curve
    axis [0].plot (sim_time, thermostat_cost, 'g-', label = "Thermostat cost")
    axis [0].set_title ("Thermostat cost")
    axis [0].set_xlabel ("Time (h)")
    axis [0].set_ylabel ("Cost (€)")
    axis [0].set_xlim ([0, 168])
    axis [0].grid (True)
    axis [0].legend ()
    # Total cost and thermostat cost curves.
    axis [1].plot (sim_time, total_cost_in_time, 'r-', label = "Total costs")
    axis [1].plot (sim_time, thermostat_cost, 'b-', label = "Thermostat cost")
    axis [1].set_title ("Costs plots")
    axis [1].set_xlabel ("Time (h)")
    axis [1].set_ylabel ("Costs (€)")
    axis [1].set_xlim([0, 168])
    axis [1].grid (True)
    axis [1].legend ()

## 8. Defining default values
default_simulation_time = np.linspace (0, 168, 16800) #hours (h)

# 8.1 Environmental temperature variation parameters
default_env_temp = 21 #ºC
default_env_var_per_day = 3 #ºC
default_frequency = np.pi / 12

# 8.2 Bacterial growth parameters
default_X0 = 3 #grams (g)
default_mu = 0.0315 # h^-1
default_S0 = 5000 #g sucrose
default_Yxs = 0.176 # g biomass/g sucrose (for Azotobacter vinelandii)
default_optimal_temp = 30 # ºC
default_heat_production_bacs = 0.0403 #cal/g biomass/s

# 8.3 Bioreactor dimensions (300L)
default_base_diameter = 0.67 #meter (m)
default_bioreactor_height = 0.85 # (m)
default_bioreactor_area = compute_cylinder_area(default_base_diameter, default_bioreactor_height) #m^2
default_wall_thickness = 0.003 # (m) == 3 mm
default_thermostat_area = default_bioreactor_area #m^2. Thermo jacket that covers half of the bioreactor.

# 8.4 Heating parameters
default_temp0 = 25 #ºC
default_cp = 4184 # J/kgºC. It is Cp(H2O).
default_conductivity = 16.2 # W/mK. Stainless steal
default_mass_medium = 150 #kg. 150L batch, assuming density~1kg/L
default_temp_thermostat = 50 #ºC

# 8.5 Impeller parameters
default_np = 5 # Adimensional
default_rho = 1000 # kg/m^3
default_rps = round (200/60, 4) # rps
default_impeller_diameter = round (default_base_diameter / 3, 4) # m
default_impeller_power = calculate_impeller_power(default_np, default_rho, default_rps, default_impeller_diameter) # W
default_impeller_time = 168 #h (1 week)

# 8.6 Cost parameters
default_cost_kWh = 0.16 # euro/kWh
default_cost_Jules = default_cost_kWh/3.6e6  # euro/Jule
default_medium_csv_file = "medium.csv"
default_autoclave_time = 21 #minutes
default_autoclave_power = 27 #kW
default_pump_power = 185 # W
default_injection_time = 168 # hours


## 9. Default setting plots

def show_default_values_env_parameters ():
    """
    Usage: show default values of environment temperature variation
    """
    print ("Showing default values for environmental temperatures...")
    print (f'simulation time: 168 h')
    print (f'environmental temp: {default_env_temp} ºC')
    print (f'daily temperature variation: {default_env_var_per_day} ºC')

def show_default_values_bacteria ():
    """
    Usage: show default values of bacterial growth parameters
    """
    print ("Showing default values for Azetobacter vinelandii...")
    print (f'X0: {default_X0} g')
    print (f'mu: {default_mu} h^-1')
    print (f'S0: {default_S0} g sucrose')
    print (f'Yxs: {default_Yxs} g biomass / g sucrose')
    print (f'optimal temperature for bacterial growth: {default_optimal_temp} ºC')
    print (f'heat_production_bacs: {default_heat_production_bacs} cal / g biomass / s')

def show_default_values_bioreactor ():
    """
    Usage: show default values of bioreactor dimensions and parameters
    """
    print ("Showing default values for bioreactor dimensions...")
    print (f'base_diameter: {default_base_diameter} m')
    print (f'bioreactor height: {default_bioreactor_height} m')
    print (f'bioreactor area: {default_bioreactor_area} m2')
    print (f'wall thickness: {default_wall_thickness} m')
    print (f'thermostat_area: {default_thermostat_area} m2')

def show_default_values_warming ():
    """
    Usage: show default values for warming parameters
    """
    print ("Showing default values for heating...")
    print (f'initial bioreactor temperature: {default_temp0} ºC')
    print (f'specific heat capacity of the medium: {default_cp} J/kg H2O/ºC')
    print (f'bioreactor conductivity: {default_conductivity} W/m/K')
    print (f'mass of the medium: {default_mass_medium} kg (150L, assuming density~1 kg/L)')
    print (f'temperature of thermostat: {default_temp_thermostat} ºC')

def show_default_impeller_parameters ():
    """
    Usage: show default values for impeller parameters
    """
    print ("Showing default values for the impeller...")
    print (f'number of power (Np): {default_np} (adimensional)')
    print (f'density of the medium: {default_rho} kg/m^3')
    print (f'rps: {default_rps} rps (200 rpm)')
    print (f'impeller diameter: {default_impeller_diameter} m')
    print (f'impeller power: {default_impeller_power} W')
    print (f'impeller working time: {default_impeller_time} h')

def show_default_values_cost ():
    """
    Usage: show default values for cost parameters
    """
    print ("Showing cost parameters...")
    print (f'cost per kWh: {default_cost_kWh} €/kWh')
    print (f'medium file used: {default_medium_csv_file}')
    print (f'autoclave power: {default_autoclave_power} kW')
    print (f'autoclave time: {default_autoclave_time} min')
    print (f'injection pump power: {default_pump_power} W')
    print (f'air injection time: 168 h')

def default_plots ():
    """
    Usage: to show the plots with the default values.
    """
    # Environmental temperature variation
    temp_env = env_temp_variation (default_simulation_time,default_env_var_per_day, default_env_temp, default_frequency)

    # Bacterial growth and substrate usage
    biomass = bacterial_growth (default_simulation_time, default_X0, default_mu)
    consumed_substrate = substrate_usage (biomass, default_Yxs)

    # Temperature variation due to bacterial growth
    ta_elevation_per_biomass_hour = celcius_elevation_per_hour (default_heat_production_bacs, default_cp, default_mass_medium)
    ta_elevation_per_hour = temp_var_ratio_due2_bacterial_growth(biomass, ta_elevation_per_biomass_hour)
    ta_elevation = temp_var_due2_bacterial_growth(default_simulation_time, ta_elevation_per_hour)

    # Temperature variations values
    temp_with_no_thermostat = heat_loss_bioreactor (default_temp0, default_thermostat_area, default_conductivity, default_wall_thickness, temp_env, default_simulation_time, default_cp, default_mass_medium)
    temp_var_with_no_thermo_with_growth = temp_no_thermo_with_growth(temp_with_no_thermostat, ta_elevation)
    real_temp_bioreactor, dQdt = thermostat (default_temp0, default_temp_thermostat, temp_env, ta_elevation, default_optimal_temp, default_wall_thickness, default_thermostat_area, default_bioreactor_area, default_conductivity, default_mass_medium, default_cp,default_simulation_time)

    # Costs and energy consumed by each system
    total_thermostat_cost, cost_in_time = thermostat_cost (dQdt, default_cost_Jules)
    energy_consumed_thermostat = kWh_consumed(total_thermostat_cost, default_cost_kWh)
    total_medium_cost, cost_per_component = medium_cost (default_medium_csv_file)
    total_autoclave_cost, energy_consumed_autoclave = autoclave_cost (default_autoclave_power, default_autoclave_time, default_cost_kWh)
    total_stirring_cost, energy_consumed_stirring = stirring_cost(default_impeller_power, default_impeller_time, default_cost_kWh)
    total_injection_cost, energy_consumed_pump = air_injection_cost(default_pump_power,default_injection_time, default_cost_kWh)

    # Total costs of the simulation
    sum_costs = round (total_thermostat_cost + total_medium_cost + total_autoclave_cost + total_injection_cost + total_stirring_cost, 2)
    total_energy_consumed = round (energy_consumed_thermostat + energy_consumed_autoclave + energy_consumed_pump + energy_consumed_stirring, 2)
    real_time_costs = total_cost (cost_in_time, total_medium_cost, total_autoclave_cost, total_injection_cost, total_stirring_cost)

    # Defining the directory name with the cost of the simulation.
    dir_name = "default_values_plots_" + str (sum_costs)
    # Bacterial growth and substrate usage plots
    print ("Press Enter to see bacterial growth and substrate usage plots.")
    enter = str(input()) #Any input will work as long as the user press any bottom.
    if len (enter) == 0 or len (enter) != 0:
        print ("Showing the plots...")
        plot_bacterial_growth_substrate_usage(default_simulation_time, biomass, consumed_substrate, default_S0)
        save_plots("bacterial_growth", dir_name)
        plt.show ()
        print("=======================================================================================")
        time.sleep(1)

    # Temperature variation due to bacterial growth.
    print ("Press Enter to see the temperature variation due to bacterial growth")
    enter = str(input())
    if len(enter) == 0 or len(enter) != 0:
        print("Showing the plots...")
        plot_temp_due2_bacterial_growth(default_simulation_time, ta_elevation)
        save_plots("temperature_variation_due2_growth", dir_name)
        plt.show()
        print("=======================================================================================")
        time.sleep(1)

    # Temperature variation plot
    print ("Press Enter to see temperature variation plots.")
    enter = str(input())
    if len (enter) == 0 or len (enter) != 0:
        print("Showing the plots...")
        plot_temperature_variation(default_simulation_time, temp_env, temp_var_with_no_thermo_with_growth, real_temp_bioreactor)
        save_plots("temperature_variation", dir_name)
        plt.show ()
        print("=======================================================================================")
        time.sleep(1)

    # Cost plots and information
    print ("Press Enter to see costs of this simulation.")
    enter = str(input())
    if len (enter) == 0 or len (enter) != 0:
        print (f'The total cost of the thermostat is: {total_thermostat_cost} €.')
        print (f'The energy consumed of the thermostat is: {energy_consumed_thermostat} kWh.')
        print (f'The total cost of the medium is: {total_medium_cost} €.')
        for key, value in cost_per_component.items():
            price = round (value, 2)
            print (f'The cost for {key} is: {price} €.')
        print (f'The total cost of the autoclave is {total_autoclave_cost} €.')
        print (f'The energy consumed of the autoclave is: {energy_consumed_autoclave} kWh.')
        print (f'The total cost of stirring is: {total_stirring_cost} €')
        print (f'The total energy consumed of stirring is: {energy_consumed_stirring} kWh.')
        print (f'The total cost of aeration is {total_injection_cost} €.')
        print (f'The energy consumed of aeration is {energy_consumed_pump} kWh.')
        print ('-------------------------------------------------------------------')
        print (f'The total energy consumed is {total_energy_consumed} kWh.')
        print (f'The total cost of this simulation is: {sum_costs} €.')
        print("Showing the plots...")
        plot_total_costs(default_simulation_time, cost_in_time, real_time_costs)
        save_plots("total_costs", dir_name)
        plt.show ()

## 10. Simulink's style simulation
def get_num_input ():
    """
    To get a number input
    """
    pattern = r"^[0-9]+(\.[0-9]+)*$" # One or more digits followed by an optional decimal part. "\." means a literal point
    while True:
        num = input ()
        m = re.search(pattern, num)
        if m: # If a number is detected just return it
            return num
        else: # Else ask for another input
            print ('Please try again entering a number.')



def simulink_plots ():
    """
    Usage: to simulate the simulink GUI.
    """
    # Get initial temperature value.
    print ("Please introduce the initial temperature (ºC) for the bioreactor:")
    simulink_temp0 = np.float16 (get_num_input ())

    #Get step value.
    print ("Please introduce the step value (0.001, e.g.)")
    step = np.float16 (get_num_input ())

    #Defining new variables according to the step value given
    points = int (168 / step)
    sim_sim_time = np.linspace (0, 168, num = points)

    # Environmental temperature variation
    temp_env = env_temp_variation(sim_sim_time, default_env_var_per_day, default_env_temp, default_frequency)

    # Bacterial growth and substrate usage.
    biomass = bacterial_growth(sim_sim_time, default_X0, default_mu)
    consumed_substrate = substrate_usage(biomass, default_Yxs)

    # Temperature variation due to bacterial growth
    ta_elevation_per_biomass_hour = celcius_elevation_per_hour(default_heat_production_bacs, default_cp,default_mass_medium)
    ta_elevation_per_hour = temp_var_ratio_due2_bacterial_growth(biomass, ta_elevation_per_biomass_hour)
    ta_elevation = temp_var_due2_bacterial_growth(sim_sim_time, ta_elevation_per_hour)

    # Temperature variation
    temp_with_no_thermostat = heat_loss_bioreactor(simulink_temp0, default_thermostat_area, default_conductivity, default_wall_thickness, temp_env, sim_sim_time, default_cp, default_mass_medium)
    temp_var_with_no_thermo_with_growth = temp_no_thermo_with_growth(temp_with_no_thermostat, ta_elevation)
    real_temp_bioreactor, dQdt = thermostat(simulink_temp0, default_temp_thermostat, temp_env, ta_elevation, default_optimal_temp, default_wall_thickness, default_thermostat_area, default_bioreactor_area, default_conductivity, default_mass_medium, default_cp, sim_sim_time)

    # Costs and energy consumed by systems
    total_thermostat_cost, cost_in_time = thermostat_cost(dQdt, default_cost_Jules)
    energy_consumed_thermostat = kWh_consumed(total_thermostat_cost, default_cost_kWh)
    total_medium_cost, cost_per_component = medium_cost(default_medium_csv_file)
    total_autoclave_cost, energy_consumed_autoclave = autoclave_cost(default_autoclave_power, default_autoclave_time, default_cost_kWh)
    total_stirring_cost, energy_consumed_stirring = stirring_cost(default_impeller_power, default_impeller_time, default_cost_kWh)
    total_injection_cost, energy_consumed_pump = air_injection_cost(default_pump_power, default_injection_time, default_cost_kWh)

    # Total cost and energy consumed
    sum_costs = round (total_thermostat_cost + total_medium_cost + total_autoclave_cost + total_injection_cost + total_stirring_cost, 2)
    total_energy_consumed = round (energy_consumed_thermostat + energy_consumed_autoclave + energy_consumed_pump + energy_consumed_stirring, 2)
    real_time_costs = total_cost(cost_in_time, total_medium_cost, total_autoclave_cost, total_injection_cost, total_stirring_cost)

    # Defining the directory name to save the plots: identified by the user given values and total costs
    dir_name = "simulink_like_plots_" + str(simulink_temp0) + "_" + str(step) + "_" + str (sum_costs)

    # Bacterial growth and substrate usage plots
    print("Press Enter to see bacterial growth and substrate usage plots.")
    enter = str(input())
    if len(enter) == 0 or len(enter) != 0:
        print("Showing the plots...")
        plot_bacterial_growth_substrate_usage(sim_sim_time, biomass, consumed_substrate, default_S0)
        save_plots("bacterial_growth", dir_name)
        plt.show()
        print("=======================================================================================")
        time.sleep(1)

    # Temperature variation due to bacterial growth.
    print("Press Enter to see the temperature variation due to bacterial growth")
    enter = str(input())
    if len(enter) == 0 or len(enter) != 0:
        print("Showing the plots...")
        plot_temp_due2_bacterial_growth(sim_sim_time, ta_elevation)
        save_plots("temperature_variation_due2_growth", dir_name)
        plt.show()
        print("=======================================================================================")
        time.sleep(1)

    # Temperature variation plot
    print("Press Enter to see temperature variation plots.")
    enter = str(input())
    if len(enter) == 0 or len(enter) != 0:
        print("Showing the plots...")
        plot_temperature_variation(sim_sim_time, temp_env, temp_var_with_no_thermo_with_growth, real_temp_bioreactor)
        save_plots("temperature_variation", dir_name)
        plt.show()
        print("=======================================================================================")
        time.sleep(1)

    # Costs plot and information
    print("Press Enter to see costs of this simulation.")
    enter = str(input())
    if len(enter) == 0 or len(enter) != 0:
        print(f'The total cost of the thermostat is: {total_thermostat_cost} €.')
        print(f'The energy consumed of the thermostat is: {energy_consumed_thermostat} kWh.')
        print(f'The total cost of the medium is: {total_medium_cost} €.')
        for key, value in cost_per_component.items():
            price = round(value, 2)
            print(f'The cost for {key} is: {price} €.')
        print(f'The total cost of the autoclave is {total_autoclave_cost} €.')
        print(f'The energy consumed of the autoclave is: {energy_consumed_autoclave} kWh.')
        print(f'The total cost of stirring is: {total_stirring_cost} €')
        print(f'The total energy consumed of stirring is: {energy_consumed_stirring} kWh.')
        print(f'The total cost of aeration is {total_injection_cost} €.')
        print(f'The energy consumed of aeration is {energy_consumed_pump} kWh.')
        print('-------------------------------------------------------------------')
        print(f'The total energy consumed is {total_energy_consumed} kWh.')
        print(f'The total cost of this simulation is: {sum_costs} €.')
        print("Showing the plots...")
        plot_total_costs(sim_sim_time, cost_in_time, real_time_costs)
        save_plots("costs_plot", dir_name)
        plt.show()

## 11. Menu
def display_menu ():
    """
    Usage: display the menu
    """
    print ("This is a program to simulate a bioreactor with thermo jacket.")
    print ("Please indicate the type of simulation you want to do:")
    print ("1) Using all default values")
    print ("2) Simulink-like simulation: you can set the initial temperature for bioreactor and fixed_step value.")

    # Get option
    not_over = True
    while not_over:
        print ("Please enter your choice:")
        option = str(input())
        if option == "1":
            print ("We are going to simulate the bioreactor using all default values...")
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_env_parameters()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_bacteria()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_bioreactor()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_warming()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_cost()
            print("=======================================================================================")
            time.sleep(1)
            show_default_impeller_parameters()
            print("=======================================================================================")
            time.sleep(1)
            default_plots()
            print("=======================================================================================")
            time.sleep(1)
            print ("The simulation is over.")
            not_over = False
        elif option == "2":
            print ("We are going to simulate the bioreactor simulink-like method...")
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_env_parameters()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_bacteria()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_bioreactor()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_warming()
            print("=======================================================================================")
            time.sleep(1)
            show_default_impeller_parameters()
            print("=======================================================================================")
            time.sleep(1)
            show_default_values_cost()
            print("=======================================================================================")
            time.sleep(1)
            simulink_plots()
            print("=======================================================================================")
            time.sleep(1)
            print("The simulation is over.")
            not_over = False
        else:
            print ("Please introduce a valid option.")

# Show the menu after running the program

display_menu()





















