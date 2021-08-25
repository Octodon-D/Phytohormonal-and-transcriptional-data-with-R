# Phytohormonal and transcriptional data with R
*Within this notebook I will give examples for the data processing and analysis I did with R as part of my doctoral thesis.*

In the first part of my doctoral thesis I investigated the activation of different signalling pathways in response to oviposition-mediated priming of anti-herbivore defences (in a nutshell: how the plant defence against feeding larvae is prepared by prior insect egg deposition, of which the larvae will hatch) in different plant species interacting with different lepidopteran herbivores. 
Therefore I measured the accumulation of phytohormones and associated transcripts of defence related genes after different treatments including insect egg deposition and/or larval feeding. Here I show the script of the processing and analysis I did with these transcriptional and phytohormonal data. In both scripts the data is processed (organized, arranged, value calculation), then statistical analyses (linear mixed models) are conducted as well as visualizations (boxplots, barplots). 

## Transcriptional analysis: [qPCR_processing_analysis.R](qPCR_processing_analysis.R)
Transcripts of defence related genes were quantified with a real time PCR (qPCR). The output of the measurement machine was pre-processed with LinRegPCR to have the ct-values for each sample (96 samples per run). For each biological sample three technical replicates were conducted, why also an outlier correction was included in the script.  
### Script overview:
* Setup
* Data processing (part 1)
* Outlier correction 
* Data processing (part 2)
* PCR quality evaluation 
* Calculations
* Statistics
* Visualization: Boxplots, Barplots

## Phytohormone analysis: [pytohormone_preprocessing_analysis.R](pytohormone_preprocessing_analysis.R)
The phytohormones of interest: salicylic acid (SA), abscisic acid, jasmonic acid (JA) and jasmonic acid-isoleucine (JA-Ile) were measured were using an Ultra-High-Performance-Liquid-Chromatography (UPLC) coupled to a Time Of Flight Mass Spectrometer (Q-ToF-ESI). The output of the UPLC for each sample is a chromatogram. Phytohormone levels were quantified according to the peak area of the plant-derived phytohormones (defined peaks within the chromatogram) relative to the peak area of the internal standards. The concentrations per sample were normalized according to the fresh weight of the leaf tissue samples. With the peak areas as input data, the script conducts these calculations. 
### Script overview:
* Setup
* Data processing
* Statistics
* Visualization: Boxplots, Barplots, Heatmap
