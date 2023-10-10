#Small mammals and associated infections in China: a systematic review and spatial modelling analysis

##Compiled standalone software and/or source code

###1. System requirements
####1.1 All software dependencies and operating systems (including version numbers)
Software R v4.0.3. R is available as Free Software under the terms of the Free Software Foundation's GNU General Public License in source code form.

####1.1 Versions the software has been tested on
Software R v4.0.3.

####1.3 Any required non-standard hardware
None.

###2.Installation guide
####2.1 Instructions
Open the R programming language environment by running the R interpreter in a terminal or RStudio. Select the R file or code snippet to be executed, either using an editor or entering it through the command line. Execute the code by pressing Enter or clicking the run button for the line of code currently selected. R will output the resulting data based on the executed code. Debug your code using control flow statements and debugging tools to check for errors and correct them. Save the results of R programming to a file for future use or sharing. Close the R programming language environment by exiting the R interpreter and closing any related windows or programs.

####2.2 Typical install time on a "normal" desktop computer
Installation of R programming is estimated to take up to 20 minutes in the fastest case.

###3 Demo
###3.1 Instructions to run on data
(01)rat_fenbu_study.csv

The distribution of small mammals at the county level.
**ADCODE99** represents the specific county

**Mus musculus**represents the absence of Mus musculus

**Rattus norvegicus**represents the absence of Rattus norvegicus

The number "1" indicated the specific small mammal detected in the county, while "0" without the small mammal in the county.

(02)Data of disease.xlsx

**ADCODE99**represents the specific county

**Presence of hantaviridae in rat**represents the presence of hantavirus in small mammals
	
**Incidence of hantaviridae**represents the incidence of hantavirus in the county
	
**Presence of Leptospira in rat**represents the presence of leptospira in small mammals

**Incidence of Leptospira**represents the presence of leptospira in small mammals

###3.2 Expected output

(03)pred_gamma.csv

**ADCODE99**represents the specific county

**incidence**represents the incidence county of hantavirus

**pred1**represents the predicted incidenceof hantavirus in the first time
**pred2**represents the predicted incidenceof hantavirus in the second time
**pred3**represents the predicted incidenceof hantavirus in the third time
......**pred100**represents the predicted incidenceof hantavirus in the 100 time

###3.3 Expected run time for demo on a "normal" desktop computer
About5 hours.

###4 Instructions for use
###4.1 How to run the software on your data
Load Supplementary Code 1. R code on R, then read in the data and follow the code flow.Run your code in R or RStudio. This can be done by pressing "Ctrl+Enter" or clicking the "Run" button on the toolbar.

###4.2 Reproduction instructions
The code can run smoothly and be copied completely.
Open R or RStudio.
Create a new script file or open an existing script file.
Enter your code in the file.
Save your script file.