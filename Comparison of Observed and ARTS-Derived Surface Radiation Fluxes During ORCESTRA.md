# Comparison of Observed and ARTS-Derived Surface Radiation Fluxes During ORCESTRA

## Scientific Goal
To compare surface radiation fluxes measured at the BCO to surface radiation fluxes derived from ARTS during the western Atlantic leg of ORCESTRA in order to access the performance of the instruments.

## Method
- [x] Collect radiosonde launches from BCO during western Atlantic leg of ORCESTRA
- [x] Group radionsondes into 'clearsky' and 'non-clearsky'
- [x] The mean midday clear-sky profile during the Campaign was chosen. 
- [x] Calculate zenith sun position with respect to BCO for each sonde and their launch time.
- [x] Ingest clearsky sondes into ARTS and retrieve estimates of surface radiation fluxes for each sonde
- [x] Collect retrievals of surface radiation fluxes measured at the BCO
- [x] Compare ARTS-derived surface fluxes to measured surface fluxes on clearsky days

## Results
- Rough clear sky radiosonde plot
![](https://pad.gwdg.de/uploads/ae8e0f2c-9f09-481c-9bec-1fe40d3c445e.png)

- Rough cloudy radiosonde plot
![](https://pad.gwdg.de/uploads/a8775db2-1a6d-49a6-b2ce-98e2d3ef589e.png)

I may have to compile noon time cases or set the angle of the Sun.

Create a function that calculates the longitude of the sun based on the time of the day at the BCO.

## Midday Clear-sky Profiles:
1. 2024:09:07 17:02
2. 2024:09:10 16:48
3. 2024:09:11 16:50
4. 2024:09:12 13:47
5. 2024:09:12 16:46
6. 2024:09:13 16:49
7. 2024:09:14 14:02
8. 2024:09:16 13:49
9. 2024:09:20 14:08
10. 2024:09:20 16:49
11. 2024:09:22 14:04

Mean Midday Sounding

![](https://pad.gwdg.de/uploads/31d110fe-4ac7-4bde-b10e-87ca1e6e7cfd.png)


## Code that converts launch time of sounding to Sun longitude
```python=
import numpy as np
from pprint import pprint


class TimeOfSounding:
    """Creates a sounding time object and allows for conversion to sun angle and other things."""

    def __init__(self, time: str, rotation_rate: int | float = 0.25):
        """
        Initial the time of the sounding.
        :param time: The launch time of the sounding.
        :param rotation_rate: The rate of rotation of the planet in degrees per minute. Default is 0.25 degrees/min.
        """

        self.date = np.datetime64(time)
        self.rotation_rate = rotation_rate  # minutes per hour

        self.sun_longitude = 0
        self.sun_latitude = 0

        self.noon = self.date.astype('datetime64[D]') + np.timedelta64(12, 'h')

    def time_to_solar_angle(self):
        """Gets the time with respect to noon for a given date and calculates the sun angle above zenith."""

        time_difference = self.date - self.noon

        solar_angle = (time_difference / np.timedelta64(1, 'm')) * self.rotation_rate

        return solar_angle


radiosonde_1 = TimeOfSounding(time="2024-09-07T09:00")
print(radiosonde_1.time_to_solar_angle())  # -45.0°
```

## Radiation Fluxes
Below is a plot of surface fluxes using the model atmosphere. It follows that we can take just the flux at the surface and call it a day. But one thing we need to retrieve are profiles of atmospheric constituents in order for the ARTS calculation to be accurate.

![](https://pad.gwdg.de/uploads/32653a16-e82b-41d2-b3e4-f76badbafaf9.png)

### ARTS derived flux profiles from radiosondes

![](https://pad.gwdg.de/uploads/c267e3e6-b341-4557-a94c-c629686d06d4.png)

ARTS is happy only with positive sun angles, i.e. times beyond noon. 

# Comparison to Observed (First Attempt)

Below, we compare the observed shortwave downwelling radiation to the ARTS-derived surface radiation fluxes.

![](https://pad.gwdg.de/uploads/27760d22-2339-4f23-851b-c13c6daae082.png)

The main takeaway from this result is the challenge that comes from aligning the timezone and zenith position of the sun.

Take mean of time during sonde lifetime from obs

### Updated Plot
![](https://pad.gwdg.de/uploads/2263b94b-4e14-4e15-943d-8004612782e8.png)

Clouds likely impacted the soundings from the 13th and the 16th so we can remove them.





## PowerPoint Presentation Link
https://mycavehilluwi-my.sharepoint.com/:p:/g/personal/chavez_pope_mycavehill_uwi_edu/IQDhdA0arYSuRpx5PdABWvYXAW2zXod5IKl9HIFsuzgajPw?e=GBjuMn

### Daniel's Presentation Plan (feel free to edit Chavez)
***Introduction*** 
- Project Name
- Overview of BCO
- Overview of ORCESTRA (and BCO's role in ORCESTRA)
***Motivation***
- BCO representing 'half a climate' of data
- Discuss how quality control is now an important task
- The ORCESTRA campaign provides a solid foundation for using ARTS to test the accuracy of observations
***Data and Method***
- Explanation of BCO instruments that we compare to ARTS derived variables
- Classification of Launched Radiosondes as either clearsky or non-clearsky (and how that determination was made)
- Mention how this classification can be made through the satellite images and regions of high humidity in the radiosonde profiles
- Acquisition of CO2 from https://cds.climate.copernicus.eu/datasets/satellite-carbon-dioxide?tab=overview value of 422 ppm from https://www.climate.gov/news-features/understanding-climate/climate-change-atmospheric-carbon-dioxide#:~:text=In%20May%202024%2C%20carbon%20dioxide,ppm%2C%20also%20a%20new%20record.
- Ingestion of clearsky profiles into ARTS to calculate atmospheric flux profiles (and of course, the surface)
- Mention that radiosonde profiles have to arranged such that they increase monotonically, (Invoked with types: pyarts3.arts.GeodeticField3, kwargs = { name: str, data: ndarray, grid_names: list, grids: list }), because the ARTS GeodeticField requires a monotonically increasing profile.
- 
***Results and Discussion***
- Plot of all clearsky radiosondes (could go into Method if so desired) and their mean
- Plot of atmospheric flux profiles and their mean
- Plot of surface radiation fluxes compared to BCO observed radiation fluxes
- Any statistical comparisons that could be useful
- Segway into why things exist the way they do in the results
- Discussion of limitatations (not using all of the trace gases, aerosols, designation of clearsky conditions)
    - The way ARTS was used leaves much room for more realistic calculations of OLR and SW downward flux
***Conclusions and Assessment***
***References***
***Supplementary Figures*** (if any)









