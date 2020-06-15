Updates and fixes of V2:

1) Created SS-Version of DATA to read it in for all model versions.
*NOT USED CURRENTLY INSTEAD USED "IF" term in ScriptFile.m
**REASON: for creating Data file used different Script called Data_6hourly_Creator
SO May be no need to use this approach.

2) Created if statement in DETECT SS and FAKE SS,
 which alows to set SWC/Gness/AtmPressure/SoilT mean and constant.

3) Fixed input of Catm to be 356 at all times without change.

4) Main_Program currently not updated to acomodate these changes, however,
plotting will work well.