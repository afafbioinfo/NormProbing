@ saaidi_LIX_2018, Updated October 2019

# NormProbing

NormProbing is a code for parsing QuShape outputs.
It computes and outputs normalized reactivity data while giving the user the possibility of choosing the nature of the nucleotides to be taken into account 
for the normalization.


## Executing NormProbing

NormProbing can be invoked through the command line: 

      python2.7 NormProbing.py

The method will run with a configuration specified within `NormProbing.Config` as explained bellow.

## Input files

### QuShape output folder

The required name/format for files is RNA_reagent_number.txt 

**Example:**

an example of input data: [mSP RNA probed with three reagents three times](https://github.com/afafbioinfo/Macro_CE/tree/master/Qu-shape-output) is provided 
       
Parameters can be set/modified through the file Macro_Qushape.Config as follow:


## Outputs

###  CE_reactivity_profiles folder:
contains the normalized reactivities for different experiments.

#### Reactivity folder:
contains final normalized reactivities by RNA and by reagent.


## Configuration
Most configuration options are set by modifying the content of the `NormProbing.Config` file.

### Main options


`Program`: 
Task: 2
==> The program by default performs normalization
if Task set to 2, the reactivity average calculation is performed in addition.

[Paths]
InputFile: Qu-shape-output
Iofile_Extenstion: txt
Outputreactivityfile: CE_reactivity_profiles

[Sequence]
Start:1
End:259
==> Different experiments come with different ranges, specify one COMMON window for all the considered experiments.

[Normalization]
Method: 
==> For short RNA (length<=300 bps) [resp.long RNA (length>300 bps] the Norm1 (Norm2) is set by default. 
The user can also either choose Norm1 or Norm2.
Norm1: The 2% of the most reactive peaks are removed. The normalization term correponds to the average
of the next 8% of the most reactive peaks.
Norm2: Quantile normalization
Reactivity values above 1.5 times the interquartile range are removed.
Then, the 10% of the highest remaining reactivities are averaged to serve as a nomalization term.

[Conditions]
Nucleotides: G,C
==> Specify nucleotides to be considered for the normalization, if empty all nuclotides are selected by default.
 nucleotides should be separated by ','
 
Nucreadout: A,C
==> Specify the nucleotides to be read out, if empty the reactivity will be calculated for all nucleotides by default.

Lowervalue: -0.3 
==>The max value for null normalized values.

Threshold=0.2
==>The average distance between two values from individuel experiment.

Desactiv_threshold=0.6
==> The reativity value from which the previous filetring by threshold is no more applied. 

############# Rules for reactivity calculation:
# if at least two values in the set are different from -10 , compute the average value. otherwise, return -10.
# if only one value equals to -10 , the mean is calculated from the other remaining values 
Get couple of 2 by 2 values then calculate their distance to chek if there is one distance above the 'Threshold'.
This filter is applied only for values below 'Desactiv_threshold'.

###################################################################################################################


Output: 
