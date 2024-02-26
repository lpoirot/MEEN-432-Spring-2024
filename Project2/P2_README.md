# MEEN 432 Spring 2024
Repository for MEEN 432 Automotive Engineering
Luke Poirot - UIN: 730000185

# Project 2 Week 1 Update
Begin by downloading the file in the project 2 folder called "P2_Week1_Update" as well as the rotate function.
Run the code.
You should see a car that makes its way around a NASCAR-like track, with its path traced behind it.

## Week 1 Feedback (3.8/5)
The rotate function was used in your script but it was not added to the submission so it did not run the simulation properly. I assume you have the function and forgot to submit it. For Week 2, start developing a lateral dynamic model of a vehicle that contains subsystems that are listed in the Week 2 document.

# Project 2 Week 2 Update
Instructions for running:
1. Begin by downloading all files in the project 2 folder.
2. Run the initialization file "init" first.
3. After, run the file labeled "P2_Week2_Update".
4. After entering the number of desired laps, the code will produce the car's speed while its driving by plotting it.

Work done:
1. Imported the Simulink demo file including the driver model, the lateral dynamics model, and transformation/rotation subsystems. 
2. Defined the tloops function in the raceStat.m file.
3. Adjusted parameter to ensure the car remains on the track (delta equation).
4. Adjusted parameters to ensure fastest lap time (delta equation and sped up the car on straight aways and slowed it down on curves).
5. Integrated the Simulink file with the P2 week 1 code to plot the results to visually display the model.
6. Added in prompts and an input feature for the user.
