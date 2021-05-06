# NormProbing (previously named Macro-CE)

NormProbing is a code for parsing QuShape outputs.
It computes and outputs normalized reactivity data while giving the user the possibility of choosing the nature of the nucleotides to be taken into account 
for the normalization.

@saaidi_LIX_2018

## Dependencies
      - python 2.7
      - numpy

## Executing NormProbing

NormProbing can be invoked through the command line: 

      python2.7 NormProbing.py

The method will run with a configuration specified within `NormProbing.Config` as explained bellow.

## Input files

### QuShape output folder

The required name/format for files is RNA_reagent_number.txt 

Note that the filename is important, as it will be used as a base name for the normalization and the output files.

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


`Task`: The default program performs normalization.
If in addition `Task` is set to  2, the calculation of the reactivity average is carried out in addition.

### Paths options

`InputFile`: Main input directory 

`Iofile_Extenstion`: format file in the input directory

`Outputreactivityfile`:  Main output directory 

### Sequence options

For the same RNA sequence, different experiments may cover different ranges of the sequence. The following options guarantee a COMMON range for all the considered experiments.

`Start`: Position of the first nucleotide to consider from the RNA sequence

`End`:Position of the last nucleotide to consider from the RNA sequence

### Normalization options

`Method`: For short RNAs (length<=300 bps) [resp.long RNA (length>300 bps] the Norm1 (resp.Norm2) is set by default. 

The user can also either choose Norm1 or Norm2.

      - Norm1: The 2% of the most reactive peaks are removed. The normalization term correponds to the average of the next 8% of the most reactive peaks.

      - Norm2: Parametric normalization; Reactivity values above both 1.5 times the interquartile  range and the 75th percentile are removed. Then, the 10% of the highest remaining reactivities are averaged to serve as a nomalization term.

###  Conditions options

`Nucleotides`: Specify the category of nucleotides to be considered for the normalization. If no value is specified, all nucleotides are taken into account by default. Nucleotides should be separated by ','.
 
`Nucreadout`: Specify the nucleotides to be read out. If no value is specified the reactivity will be calculated for all nucleotides by default.


`Threshold`: The distance between two reactivity values from two experiments for a given nucleotide

`Desactiv_threshold`: The reactivity value from which the filetring by threshold is no more applied. 

## Rules used to compute reactivity 

-  If at least two values in the set are different from -10 , compute the average value. otherwise, return -10.

-  Get couple of 2 by 2 values then calculate their distance to chek if there is one distance above the 'Threshold'.
This filter is applied only for values below 'Desactiv_threshold'.


