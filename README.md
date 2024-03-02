# Gravity/Coulumbic Force Simulator
This repository contains the code for a rudimentary gravity simulator created using OOP in C++. 

The NBodys2,cpp file is the C++ file for my code. It takes a CSV file as an input with a header line and the subsequent lines containing the required information:(mass, x, y, z, vx, vy, vz, charge). A file named 'settings.csv' has the data from the sun, all eight planets, and Pluto in the correct format. The data was taken from JSL Horizon on January 1, 2022, 00:00:00. The information can be found under resources. 

The outputs are N+1 CSV files('nplanet. txt' or 'nelectron.txt.'). N files contain the information of each planet/electron after a given time, and the remaining file is a 'data .tmp' fil3, with the data dumped. It can be used to plot our data using the GNU plot. The code will also say 'DONE' and give the time taken once finished. 

The file labelled "documentation.zip" contains the doxygen documentation. Click on the file named "index.html" to access the HTML.
There are two code files. NbodyCmd takes commands directly from the command line as arguments. NbodySim(input) runs the program and then takes in inputs. 
