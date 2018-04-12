from ConfigParser import SafeConfigParser 
import sys

#Connect to  the config file
config= SafeConfigParser()
config.read("Macro_Qushape.Config")

#Get parameters
Path=config.get("Paths","InputFile")
OutputPath=config.get("Paths","Outputreactivityfile")
FileExtension=config.get("Paths","Iofile_Extenstion")
Selected_Nucleotide=(config.get("Conditions","Nucleotides")).split(',')
Lowervalue=float(config.get("Conditions","Lowervalue"))
Threshold=float(config.get("Conditions","Threshold"))
Desactiv_threshold=float(config.get("Conditions","Desactiv_threshold"))
