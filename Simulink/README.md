# 1. Introduction 
This program is built to simulate a bioreactor with a thermo 
jacket for temperature control. For more detail please check the pdf file.

# 2. Environment requirement
To run the program, two libraries need to be installed: 
```
numpy
matplotlib.pyplot      (Python modules)
```

# 3. Usage
The user only need to run the entire program file in an IDE, such as Vs code, PyCharm or Spyder.
And follow the instructions shown in the console. After running the program, this information
should appear on the screen:

This is a program to simulate a bioreactor with thermo jacket.
Please indicate the type of simulation you want to do:
1) Using all default values
2) Simulink-like simulation: you can set the initial temperature for bioreactor and fixed_step value.
Please enter your choice:  

After finishing one simulation, if user want to do another, just run this line in the console:
``` 
display_menu ()      (python console)
```

### 3.1. Default mode
When the user chooses this mode, default values will be used (Azetobacter vinelandii, batch 
culture in 150L medium in 300L bioreactor with thermo jacket). These default values will appear on the screen.

### 3.2. Simulink mode
This mode is designed to approximate the simulink GUI (Graphical User Interface). 
Here the user can tune 2 parameters:
```
1. The starting temperature of the bioreactor
2. Fixed step value.   
``` 
As this simulation is for 1 week (168h), there is no option to change the simulation time. 
However, if the user really want to change it, they should go to the part of the program changing
the default values.

# 4. Program description

We can decompose this program in 12 parts:
0. Modules importation: numpy, matplotlib.pyplot, time and os.
1. Environmental temperature variation functions
2. Bacterial growth and substrate usage functions.
3. Temperature variation functions.
4. Environmental temperature effects functions.
5. Thermostat-associated functions.
6. Costs and energy consumed functions.
7. Plotting results and figure saving functions.
8. Defining default parameters
9. Default-mode functions
10. Simulink-like mode functions.
11. Menu  

NOTE: The folder where to save the plots will be created after running the program and it is identified
by the mode and the total cost of the simulation. For example:
``` 
default_values_plots_89.76  (directory)
simulink_like_plots_15.0_0.001_90.72 (directory)
```

NOTE2: If the user want to use another medium, they should provide a .csv file in thie format:
``` 
compond_name,€/kg,grams_used        (csv file)
```
In any other format, this program won't calculate the cost of the medium.  
And please make sure this csv file is in the same directory that the .py file.  
For an example see:  "medium.csv"

# 5. Credits
This program is developed by Yixi Zhang for a task in "Modelización de Sistemas Biológicos"
(Biotecnología, ETSIAAB, UPM). For any questions, please contact me: yixi.z@alumnos.upm.es.  

The culture medium and conditions collected from this paper:  
García-Molina, G., Echavarri-Erasun, C., del Barrio, M., 
González-García, S., Rubio, L. M., Pita, M., 
& De Lacey, A. L. (2025). 
Electroenzymatic N₂ reduction to ammonia by nitrogenase from Azotobacter vinelandii immobilized on low density graphite electrodes under direct electron transfer regime. 
Electrochimica Acta, 529, 146342. 
https://doi.org/10.1016/j.electacta.2025.146342