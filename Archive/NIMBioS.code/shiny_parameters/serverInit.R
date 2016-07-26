###
### These are initialization steps that occur within the shinyServer function
### split them out to streamline files 
###
    
VolumeRoots = c(home="~",getVolumes()(),temp=tempdir(),wd="./")   #function from shinyFiles

rValues <- reactiveValues(ssClass=ssClassInit(),
                          scenarioNumber=1,
                          lstScenario=1,
                          history=NULL,
                          msg=NULL)

####this reactiveValue contains values that are needed for file operations
####cannot use reactive due to constraints imposed by shinyFiles
supportValues <- reactiveValues(ssLoadEnv=new.env(),  #environment to load an rdata file into
                                objLabel=NULL,        #name of ssClass object for saving
                                simroot = NULL,
                                OS = .Platform$OS.type, #used when setting up simulation runs
                                simexec = c("fsc251","fsc252","fsc25221","fsc251.exe","fsc252.exe","fsc25221.exe")
                                )
