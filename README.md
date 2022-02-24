# MATLAB/Octave Implementation of Recommendation ITU-R P.528-5

This code repository contains a MATLAB/Octave software implementation of  [Recommendation ITU-R P.528-5](https://www.itu.int/rec/R-REC-P.528/en) with a propagation prediction method for aeronautical mobile and radionavigation services in the frequency range 100 - 30000 MHz.  

This is a translation of the original reference C++ implementation of this Recommendation available at [NTIA/p528](https://github.com/NTIA/p528) provided by the US National Telecommunications and Information Administration [NTIA](https://www.ntia.gov). This version of the code corresponds to the reference MATLAB/Octave version approved by ITU-R Working Party 3K and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx).

The following table describes the structure of the folder `./matlab/` containing the MATLAB/Octave implementation of Recommendation ITU-R P.528.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`tl_p528.m`                | MATLAB function implementing Recommendation ITU-R P.528-5          |
|`validate_p528.m`          | MATLAB script used to validate the implementation of Recommendation ITU-R P.528-5 against the reference results provided in the subfolder `./Data Tables`            |



## Function Call

~~~
result = tl_p528(d__km, h_1__meter, h_2__meter, f__mhz,  T_pol, p);
~~~


## Required input arguments of function `tl_p528`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `d__km`               | scalar double | km   | 0 < `d`   | Great circle path distance between terminals  |
| `h_1__meter`      | scalar double | m    | 1.5 ≤ `h_1__meter` ≤ 20000 | Height of the low terminal |
| `h_2__meter`      | scalar double | m    | 1.5 ≤ `h_2__meter` ≤ 20000 | Height of the high terminal |
| `f__mhz`          | scalar double | MHz    | 100 ≤ `f__mhz` ≤ 30000   | Frequency|
| `T_pol`           | scalar int    |       |             |  Polarization <br> 0 = horizontal <br> 1 = vertical |
| `p`          | scalar double | %    | 1 ≤ `p` ≤ 99   | Time percentage|



 
## Outputs ##

Outputs are contained within a defined `result` structure:

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `A__db`    | double | dB    | Basic transmission loss |
| `d__km`	| double  |	km	|Great circle path distance. Could be slightly different than specified in input variable if within LOS region |
| `A_fs__db`    | double | dB    | Free-space basic transmission loss |
| `A_a__db`    | double | dB    | Median atmospheric absorption loss |
| `theta_h1__rad`    | double | rad    | Elevation angle of the ray at the low terminal|
| `propagation_mode`    | int |    | Mode of propagation <br>1 = Line of Sight<br> 2 = Diffraction<br> 3 = Troposcatter|
| `rtn`    | int |    | Return flags / Error codes|


## Software Versions
The code was tested and runs on:
* MATLAB versions 2017a and 2020a
* Octave version 6.1.0

## References

* [Recommendation ITU-R P.528](https://www.itu.int/rec/R-REC-P.528/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)

* [NTIA/p528](https://github.com/NTIA/p528) 


