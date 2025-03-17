rm(list = ls())
#witch_folder = "../witch" #Where you're WITCH code is located
witch_folder = '/Users/cindyazuero/Documents/WITCH model/witch-master'
#main directory of your results files
#main_folder <- witch_folder # by default, the witch source folder
# With DAC branch
main_folder <- '/Users/cindyazuero/CMCC Dropbox/Cindy Azuero/WITCH/Results/CDR_curves-Oct2024/'
# With master
#main_folder <- "/Users/cindyazuero/CMCC Dropbox/Cindy Azuero/WITCH/Results/CDR_curves-Oct2024/ctax-oct15_2024"
subdir = c("./DAC_branch/DAC_branch_ctax_Dec09_2024", "./ctax-oct15_2024") #can be multiple directories

restrict_files = c("results_") #to all scenarios matching partly one of its arguments
exclude_files = c("")
removepattern = c("")

yearmin = 1980
yearmax = 2100

#If you want to have significant separations or parts of file names, specify file_separate <- c(type="first|last|separate", sep="_", names="c("file_new"))
#file_separate <- c("last", "_", c("specification"))
#Name scenarios (also subsets to the ones given (otherwise it takes gdx filename) as a mapping
#scenlist <- c("results_ssp2_asia_curpol"="Current policies")

#c(lsf.str()) #show all available functions

#Initialize default options, load all witch and other functionsget
source('R/witch_functions.R')

#gdxcompaR (Standard gdxcompaR based on typical variables, otherwise edit in gdxcompaR/server.R)
runApp(appDir = "gdxcompaR/witch")
