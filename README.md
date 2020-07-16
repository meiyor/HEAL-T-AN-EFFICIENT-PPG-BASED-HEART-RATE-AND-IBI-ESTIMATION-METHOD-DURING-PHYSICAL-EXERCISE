# HEAL-T-AN-EFFICIENT-PPG-BASED-HEART-RATE-AND-IBI-ESTIMATION-METHOD-DURING-PHYSICAL-EXERCISE
This repository contains a functional code composed by Matlab and bash codes for Interbeat-interval (IBI), Heart-rate (HR) calculation, and HRV metrics from an artifactual Blood Volume Pulse (BVP) signal,  and a corresponding accelerometer. Please refer to our paper from EUSIPCO 2016 Torres, J. M. M., Ghosh, A., Stepanov, E. A., &amp; Riccardi, G. (2016, August). Heal-t: An efficient ppg-based heart-rate and ibi estimation method during physical exercise. In 2016 24th European Signal Processing Conference (EUSIPCO) (pp. 1438-1442). IEEE.

# HEAL-T: An Efficient PPG based Heart-Rate and IBI Estimation During Physical Exercise

Photoplethysmography (PPG) is a simple, unobtrusive and low-cost technique for measuring blood volume pulse (BVP) used in heart-rate (HR) estimation. However, PPG based
heart-rate monitoring devices are often affected by motion artifacts in on-the-go scenarios, and can yield a noisy BVP signal reporting erroneous HR values. Recent studies have proposed spectral decomposition techniques (e.g. M-FOCUSS,  Joint-Sparse-Spectrum) to reduce motion artifacts and increase
HR estimation accuracy, but at the cost of high computationalload. The singular-value-decomposition and recursive calculations present in these approaches are not feasible for the implementation in real-time continuous-monitoring scenarios. In
this paper, we propose an efficient HR estimation method based on a combination of fast-ICA, RLS and BHW filter stages that avoids sparse signal reconstruction, while maintaining a high HR estimation accuracy. The proposed method outperforms
the state-of-the-art systems on the publicly available TROIKA data set.

## Getting Started

HEAL-T application can be executed  in bash following the command:
```bash
 sh runHEALT.sh DATA_PATH NAME_OUTPUT_folder METHOD_SELECTOR FS Fp1 Fz1 Fp2 Fz2 inc1 inc2 overlap WIN S_sel fileo.out MATLAB_PATH
```
or - **if DATA_PATH will be the same output path***
```bash
 sh runHEALT.sh DATA_PATH METHOD_SELECTOR FS Fp1 Fz1 Fp2 Fz2 inc1 inc2 overlap WIN S_sel fileo.out MATLAB_PATH
```

** YOU CAN ALSO RUN THIS CODE IN MATLAB directly **

1) Please open matlab command window prompt and execute

```matlab
   addpath(genpath(HEAL-T_path))
```

2) Now you can run

```matlab
heal_t_call(DATA_PATH,NAME_OUTPUT_folder,METHOD_SELECTOR,FS,{[Fp1 Fz1],[Fp2 Fz2]},[inc1,inc2],overlap,WIN,S_sel);
```

3) You can see the same notifications that appears in the log file from bash but in this case from the command window. (Please see below)


### Prerequisites:

Before running this code be sure adding the following dependencies folders in the thirdparty directory

1) eeglab: in this code we only use runica functions, download the last eeglabversion from https://sccn.ucsd.edu/eeglab/downloadtoolbox.php
2) MFOCUSS: This code library calculate a sparse spectrum reconstruction for stationary signals, download this from Zhiling Zhang code repository http://dsp.ucsd.edu/~zhilin/Software.html

## Built With

* MATLAB > R2016a
* bash

### Running the tests

###### The input parameters should be set as follow:

###### 1.DATA_PATH: the folder path that will be processed heal_t_call.m *

###### 2. NAME_OUTPUT_folder : ibifolder selection based on the method 0-> IBIHEALPEAK 1-> ibi.txt *

###### 3. METHOD_SELECTOR : 0 -> EMBC method 1 -> new HEAL-T method *        

###### 4. SAMPLING FREQUENCY (FS): sampling frequency value (synchronize BVP with Accel) [Hz].

###### 5. Fp1 (Hz) low cut-off first filter BHW [Hz].

###### 6. Fz1 (Hz) high cut-off first filter BHW [Hz]. 

###### 7. Fp2 (Hz) low cut-off second filter BHW [Hz].

###### 8. Fz2 (Hz) high cut-off second filter BHW [Hz].      

###### 9. inc1 filter one (Hz)

###### 10. inc2 filter two (Hz)

###### 11. overlap per each BVP segment (sec or percentage)

###### 12. windows size (WIN) BVP segment (sec or number of windows)

###### 13. spline selection (S_sel) BVP segment (sec or number of windows)

###### 14. file.out : name of the nohup output file

##### 15. MATLAB_PATH: this is the matlab binary path defined by user (i.e., in Matlab /usr/local/../bin/matlab and in mac /Applications/Matlab.../bin/matlab)

__IF YOU RUN HEAL-T FROM BASH:__

Take into accoun the folder tree description:

Main -> main code, i.e, pwd or a folder that includes .m

functionsHR -> auxiliary functions for HR processing and some functions necessary from MFOCUSS baseline evaluation and the first HEAL-T using time overlap

Data -> defined by user

eeglab -> appears in thirdparty folder please unzip it before run the code (i.e. EEGlab ICA)

### You can run the HEAL-T code based on this bash command:
```bash 
   sh runHEALT.sh DATA_FOLDER 0 1 32 0.7 2.5 0.7 3.5 0 0 1 50 10 fileo.out MATLAB_PATH
```

### Your input BVP and ACC  files with the corresponding time series should be included in your Data folder as:

1. bvp.txt
2. acc.txt

bvp and acc files with the input time series must contain at least 2 columns (one -> unix_time, second->non-normalized values). In case of the accelerometer the file should
include the time column and the three axis columns [t ; x ; y ; z] in this order. You can create this as a custom csv file delimited by colon. 

** OUTPUT: outputfile will be generated in the given log folder, with the name and the ID of each process as prefix


## The Outputs files structure

The output files will be formatted as follows_

1. The average HR output file will be a .csv file with with three columns: Col1: window time [index], Col2: HR [bpm], and Col3: HR smooth [bpm]

2. The output file containing the IBI peak-to-peak will be another .csv file with three columns Col1: Time [s] (take into account your sample frequency), Col2: HR values [bpm] , Col3: HR smoothed values [bpm] (after Smoothing spline).

e.g. PID_IBIHEALPEAK.txt

Please review the output in the log folder and you can check the following items when the script runs successfully:

* The code is capable to read inputs when the log folder is already updated by the corresponding method ("File Exist!" string appears)
* Code will end when the string "Process PID has ended successfully" appears. 
* If you want to know the main workflow of this approach in detail, please refer to workflow.eps file.

#### Note1: change the DATA_PATH as you desire : please take into account that inside this folder you should have a file called bvp.txt in which you should define your bvp with the same parameters and instructions explained above

#### Note2: Please run this code from the main directory and any input, data, or output data folder should be defined by user respecting the folder tree.

Some outputs time-series calculated from ICASSP2015 cup (http://archive.signalprocessingsociety.org/community/sp-cup/ieee-sp-cup-2015/) data are:

![alt text](https://github.com/meiyor/HEAL-T-AN-EFFICIENT-PPG-BASED-HEART-RATE-AND-IBI-ESTIMATION-METHOD-DURING-PHYSICAL-EXERCISE/blob/master/sub9train.jpg)
![alt text](https://github.com/meiyor/HEAL-T-AN-EFFICIENT-PPG-BASED-HEART-RATE-AND-IBI-ESTIMATION-METHOD-DURING-PHYSICAL-EXERCISE/blob/master/res7testprobe.jpg)

The Bland-Altman plots obtained from the ICASSP2015 data were:

![alt text](https://github.com/meiyor/HEAL-T-AN-EFFICIENT-PPG-BASED-HEART-RATE-AND-IBI-ESTIMATION-METHOD-DURING-PHYSICAL-EXERCISE/blob/master/blandalt2.png)
![alt text](https://github.com/meiyor/HEAL-T-AN-EFFICIENT-PPG-BASED-HEART-RATE-AND-IBI-ESTIMATION-METHOD-DURING-PHYSICAL-EXERCISE/blob/master/IBIest.png)

## Contributing

Please read the paper Torres, J. M. M., Ghosh, A., Stepanov, E. A., & Riccardi, G. (2016, August). Heal-t: An efficient ppg-based heart-rate and ibi estimation method during physical exercise. In 2016 24th European Signal Processing Conference (EUSIPCO) (pp. 1438-1442). IEEE (https://ieeexplore.ieee.org/document/7760486) for more details

### Authors

* **Juan Manuel Mayor Torrres**, **Arindam Ghosh**, **Evgeny Stepanov**, and **Giuseppe Riccardi**

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

see the [LICENCE.md](https://github.com/meiyor/HEAL-T-AN-EFFICIENT-PPG-BASED-HEART-RATE-AND-IBI-ESTIMATION-METHOD-DURING-PHYSICAL-EXERCISE/blob/master/LICENCE.md) file for details

### Acknowledgments

University of Trento - Signals and Interactive Sistems Lab members and alumni (http://sisl.disi.unitn.it/)
