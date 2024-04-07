# MEEN 432 Spring 2024
Repository for MEEN 432 Automotive Engineering
Luke Poirot - UIN: 730000185

# Project 3 Week 1 Update
Instructions for Running Week 1:
1. Download the files from the week 1 folder
2. Run init.m
3. Run initDriveCycle.m
4. Run initDriveCycleHighway.m
5. Run p3_runsim.m

Work Done Week 1:
1. Imported demo code for init.m, initDriveCycle.m, and p3_runsim.m
2. Extended the Urban drive cycle test to contain the entire data set and simulated the data with our model
3. Simulated our model on the Highway drive cycle data
4. In p3_w1_car, iteratively adjusted PID controller proportional, integral, and derivative gains until the vehicle error is less than 3 mph throughout the cycle
	- ***Proportional gain: 0.1 (old) -> 0.5 (new)
	- ***Integral gain:     0.0 (old) -> 0.0 (kept same)
	- ***Derivative gain:   0.0 (old) -> 0.0 (kept same)

# Project 3 Week 2 Update
Instructions for Running Week 2:
1. Delete the files downloaded from the week 1 folder (we do not want the upcoming downloads to save as new names)
2. Download the files from the week 2 folder
3. Run init.m
4. Run initDriveCycle.m
5. Run initDriveCycleHighway.m
6. Run P3init.m
7. Run p3_runsim.m

Work Done Week 2:
1. Used code from week 1 and imported demo code for P3init.m
2. In P3init.m, removed the commenting out of the plotting code for electric motor/transmission data (contour plot)
3. Outputted the angular velocity and torque from the Simulink model
4. Added in a calculation to predict the energy consumed by the vehicle when driving on the two EPA cycles

## Week 1/2 Feedback (5/5)
 good job! Looking further into Project 3, I would advise the team to do the following: 1) In the driver subsystem, start adding/developing logic for regen braking. 2) In the braking subsystem, there are more things that needto be added regarding the logic of how the brake torque/brake force is calculated. I would suggest looking at the long dynamics lecture and  taking assumptions to simplify the calculations for the brake torque but I advise to utilize some of the braking logic (locked and unlocked) 3) Start developing an Electric Motor Drive subsystem that has the components needed for the electric motor (battery, inverter, electric motor) and develop the logic that is needed to calculate the motor torque. The team should also add a Drive subsystem that takes in the outputs from the Electric Motor Drive Subsystem and multiplies the outputs by some Final Drive Gear Ratio (FDG) which will replace the current powertrain block. 

 # Project 3 Week 3 Update
 Instructions for Running Week 3:
 1. Open Matlab and Simulink
 2. Delete the files downloaded from the week 1 and 2 folder (we do not want the upcoming downloads to save as new names)
 3. Download the files from the week 3 folder and open each of the matlab and simulink files
 5. Run init.m
 6. Run initDriveCycle.m
 7. Run Finalinit.m
 8. Run p3_runsim.m

Work Done Week 3:
1. Used code from week 2 and imported the Finalinit.m demo code
2. Implemented an electric motor drive block that includes the electric motor subsystem and the battery subsystem
3. Implemented the electric motor subsystem
4. Implemented the battery subsystem
5. Updated the rest of the simulation for both urban and highway to work with the new electric motor drive block
6. Changed the powertrain block to a drive block and added additional functionalities
7. Organized the simulink models for clarity and removed unnecessary files from the submission
