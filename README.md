# cern24_TestBeam

## Repository for JLab Hall D March 2024 Large GEMTRD Test Beam

Workflow on JLab Gluon compute nodes:   
Log in to Gluon nodes with JLab computing account (must have 2FA)
```
ssh -XY [$USERNAME]@scilogin.jlab.org  
Password: [Pin+OTP]  
ssh hallgw  
Password: [Pin+OTP]  
ssh gluon[100-150]  
Password: [JLab CUE]  
```
After cloning repo, make directory links   
```
./make_gluon_links.sh  
```
And use executable files for run analysis   
```
./trdclass_halld24.sh [$RUNNUMBER] [$MAXNUMBEROFEVENTS] [$FIRSTEVENT]  
./trd_mlp_halld24.sh [$RUNNUMBER]  
```
With this workflow, raw .evio data files have already been processed into .root files that now live in the `ROOT/` directory. These .root data files are analyzed with the `trdclass_halld24.C` analysis macro. Output from there is saved in the `RootOutput/halld24` directory. The output from this in the form of a .root TTree file is passed on to `trd_mlp_halld24.C` where a rejection factor calculation is done for different particle efficiencies. Output from this NN macro is saved in the `mlpOutput/halld24` directory. Run Numbers corresponding to the HallD 2024 test beam are in the 4000-4999 range.  

The [$FIRSTEVENT] option should not be less than one. For example, if you want to analyze the first 10K events in run 4123, you would use:  
```
./trdclass_halld24.sh 4123 10000 1  
```
If you want to analyze all events in this run, the [$MAXNUMBEROFEVENTS] options should be zeroed:  
```
./trdclass_halld24.sh 4123 0 1  
```
To run with an event-by-event clustering display for the detector, simply uncomment the indicated definitions at the top of the `trdclass_halld24.C` macro:  
```
...  
//-- For single evt clustering display, uncomment BOTH:  
#define SHOW_EVTbyEVT  
#define SHOW_EVT_DISPLAY  
...  
```


See JLab's [2FA documentation](https://jlab.servicenowservices.com/sp?id=kb_article_view&sysparm_article=KB0012313&sys_kb_id=a8caee091b990910a552ed3ce54bcbe3&spa=1.)  
See JLab analysis code repository: [trd_root](https://github.com/JeffersonLab/trd_root/tree/main)  
See JLab JANA repository: [JANA4ML4FPGA](https://github.com/JeffersonLab/JANA4ML4FPGA/tree/main)  

