from ConfigParser import SafeConfigParser 
import sys

#Connect to  the config file
config= SafeConfigParser()
config.read("Macro_Qushape.Config")


#Get parameters
Task=config.get("Program","Task")
Path=config.get("Paths","InputFile")
OutputPath=config.get("Paths","Outputreactivityfile")
FileExtension=config.get("Paths","Iofile_Extenstion")
SelectedNuc=(config.get("Conditions","Nucleotides")).split(',')
Nucreadout=(config.get("Conditions","Nucreadout")).split(',')
Method=(config.get("Normalization","Method"))
Lowervalue=float(config.get("Conditions","Lowervalue"))
Threshold=float(config.get("Conditions","Threshold"))
Desactiv_threshold=float(config.get("Conditions","Desactiv_threshold"))
