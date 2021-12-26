# CardiSee
CardiSee: a supervised classification model for cardiovascular disease detection.

Overview of workflow, as well as descriptions of each file:

1) Obtaining data from MIT and University of Rochester Hospital (Physionet database). 
2) Process data (noise detection and artifact removal), done with the nt1.m file. The nt1.m program is called by the hrv1b.m program, which calculates heart rate variability metrics for each ECG measurement, by also calling the HRV.m file. The HRV.m file is a script written by Dr. Marcus Vollmer at the University of Greifswald in Germany, that calculates a series of standard HRV metrics for a single ECG measurement.
3) Feature selection using statistical analysis
4) Random forest classification model training and testing, done in the hrv1e.m file. 

<img width="786" alt="Screen Shot 2021-12-25 at 4 33 49 PM" src="https://user-images.githubusercontent.com/9027401/147396087-ef68bc7d-97fc-4506-8ced-83f464cf52df.png">

