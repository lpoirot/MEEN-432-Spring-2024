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

## Week 2 Feedback (5/5)
Good Job Team 19! This team hit all the points that I wanted to make about the final submission as you included the animation and raceStat function so kuddos! I really enjoyed the user input of how many laps that was cool! The one thing I would suggest about the user input is to use the response to determine the what the simtime should be instead of using it as the loop index. It seems that your model is able to accomplish two laps in 132 seconds, so if a user inputs 6 laps, you should adjust the simtime to accomplish those 6 laps (396 seconds) that way the raceStat function is properly displaying the correct number of loops.

# Project 2 Week 3 Update
Instructions for running:
1. Begin by downloading all files in the project 2 folder.
2. Run the initialization file "init" first.
3. After, run the file labeled "P2_Week3_Update".
4. After entering the number of desired laps, the code will produce the car's speed while its driving by plotting it.

Work done:
1. Changed the outputted time variable to be the simtime for each lap from the simulink (previously, it was a recorded real life time for each lap).

Summary of Contributions:
Because we are wrapping up a project this week, we decided to not submit each team member's contributions as seperate files as our matlab/simulink arcitecture is not currently built like this; however, we will do this going forward in future project, and we have included a summary of each member's contributions below
- Luke Rehfuss: Helped with plotting the track and track/car visuals, created user prompts and input retrieval, determined the delta_f variable for steering around the track, and helped with various coding/simulation errors
- Luke Poirot: Created a way for the car to rotate around the track, and developed ways to keep the car on the loop after the first turn, led the submissions for the team to the github, and helped with various coding/simulation errors
- Matthew Valtierra: Developed the physical track for the car to drive around, collected the data from the lateral dynamic model for outputting to the track plot and analyzing in the raceStat.m file, and helped with various coding/simulation errors
