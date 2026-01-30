# %%

# Imports
import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pa
import xarray as xr

# %%

ATLANTIC_STANDARD_TIME = np.timedelta64(4, "h")

THERMAL_SURFACE_REFLECTIVITY = 0.02  # From MODIS
VISIBLE_SURFACE_REFLECTIVITY = 0.22  # Value for healthy grass


def calc_rcemip_ozone(plev, g1=3.6478, g2=0.83209, g3=11.3515):
    r"""Compute the ozone volumetric mixing ratio from pressure.

    .. math::
        O_3 = g_1 \cdot p^{g_2} e^\frac{-p}{g_3}

    Parameters:
        plev (ndarray): Atmospheric pressure [Pa].
        g1, g2, g3 (float): Fitting parameters for gamma distribution
            according to Wing et al. (2017).

    Returns:
          ndarray: Ozone profile [VMR].

    Reference:
        Wing et al., 2017, Radiative-Convective Equilibrium Model
        Intercomparison Project

    """
    p = plev / 100
    return g1 * p ** g2 * np.exp(-p / g3) * 1e-6


class AtmosphericProfile:
    """An atmospheric profile for use with ARTS."""

    def __init__(self, input_data_path: str, sonde_index: int = 0):
        self.sonde_index = sonde_index
        self.input_data_path = input_data_path
        self.atmospheric_profiles = xr.open_dataset(self.input_data_path)

        self.atmospheric_profile = self.atmospheric_profiles.isel(sonde_id=self.sonde_index)
        self.atmospheric_profile_launch_time = None
        self.atmospheric_profile_final_time = None
        self.atmospheric_profile_duration = None
        self.noon = None
        self.arts_atmospheric_profile = None

        self.rotation_rate = 0.25  # Degrees per minute
        self.solar_longitude = 0

        self.solar = None
        self.thermal = None
        self.altitude = None

        self.downwelling_sw_rad = None
        self.surface_olr = None

    def get_launch_time(self):
        """Extracts the launch time of the radiosonde and stores it."""

        self.atmospheric_profile_launch_time = self.atmospheric_profile.launch_time.values - ATLANTIC_STANDARD_TIME

        self.noon = self.atmospheric_profile_launch_time.astype("datetime64[D]") + np.timedelta64(12, "h")

        return self.atmospheric_profile_launch_time

    def get_launch_duration(self):
        """Gets the duration of the radiosonde. Use to calculate time mean for pyranometer reading."""

        self.atmospheric_profile_final_time = self.atmospheric_profile.flight_time.max(skipna=True).values

        self.atmospheric_profile_duration = self.atmospheric_profile_final_time - self.atmospheric_profile_launch_time

        return self.atmospheric_profile_duration

    def calculate_solar_longitude(self):
        """Computes the solar angle using the launch time of the sounding."""

        time_difference = self.atmospheric_profile_launch_time - self.noon

        self.solar_longitude = np.abs((time_difference / np.timedelta64(1, "m")) * self.rotation_rate)

        return self.solar_longitude

    def create_arts_atmospheric_profile(self):
        """Create atmospheric profile that is ARTS friendly"""

        self.arts_atmospheric_profile = xr.Dataset(
            data_vars={
                "p": (("alt"), self.atmospheric_profile.p.values),
                "t": (("alt"), self.atmospheric_profile.ta.values),
                "H2O": (("alt"), self.atmospheric_profile.mr.values),
                "CO2": ("alt", np.ones_like(self.atmospheric_profile.p.values) * 422 / 1e6),
                "O3": ("alt", calc_rcemip_ozone(self.atmospheric_profile.p.values / 100)),
                "O2": ("alt", np.ones_like(self.atmospheric_profile.p.values) * 0.21),
                "N2": ("alt", np.ones_like(self.atmospheric_profile.p.values) * 0.78)
            },
            coords={
                "alt": self.atmospheric_profile.alt.values,
                "lat": self.atmospheric_profile.launch_lat.values,
                "lon": 0
            }
        )

        self.arts_atmospheric_profile = self.arts_atmospheric_profile.sortby("alt")

        return self.arts_atmospheric_profile

    def calculate_atmospheric_fluxes(self):
        """Uses ARTS to calculate the atmospheric flux profiles."""

        pa.data.download()

        fop = pa.recipe.AtmosphericFlux(
            species=["H2O-161", "O2-66", "N2-44", "CO2-626", "O3-XFIT", "H2O-ForeignContCKDMT400", "H2O-SelfContCKDMT400"],
            remove_lines_percentile={"H2O": 70},
            atmospheric_altitude=self.arts_atmospheric_profile.alt.max(skipna=True).values,
            solar_longitude=self.solar_longitude,
            surface_temperature=self.arts_atmospheric_profile.t.values[0],
            thermal_surface_reflectivity=THERMAL_SURFACE_REFLECTIVITY,
            visible_surface_reflectivity=VISIBLE_SURFACE_REFLECTIVITY
        )

        self.solar, self.thermal, self.altitude = fop(atmospheric_profile=self.arts_atmospheric_profile)

        self.downwelling_sw_rad, self.surface_olr = self.solar.down[-1], self.thermal.up[-1]

        return self.solar, self.thermal, self.altitude

    def plot_atmospheric_solar_down_flux_profile(self):
        """Plots the atmospheric solar down flux profile"""

        fig, ax = plt.subplots(figsize=(8, 12))

        flux_profile = ax.plot(self.solar.down, self.altitude / 1e3)

        ax.set(
            title="Profile of Downward Shortwave Flux",
            xlabel="Flux $Wm^{-2}$",
            ylabel="Altitude [km]"
        )

        # Remove top/right spines and move bottom/left outward
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('outward', 5))
        ax.spines['bottom'].set_position(('outward', 5))

        return flux_profile


class SurfaceFluxes:
    """An object containing the outputs from all the atmospheric profiles"""

    def __init__(self, atmospheric_profiles: list[AtmosphericProfile]):
        self.atmospheric_profiles = atmospheric_profiles

        self.profile_launch_times = []
        self.profile_final_times = []
        self.profile_durations = []

        self.sun_longitudes = []

        self.profiles_altitude = []
        self.profiles_down_solar_flux = []
        self.profiles_olr = []

        self.surface_down_solar_fluxes = []
        self.surface_olr_fluxes = []

    def _extract_fluxes(self, atmospheric_profile: AtmosphericProfile):
        """Gets surface and profile fluxes for downwelling SW and OLR."""

        atmospheric_profile.get_launch_time()
        self.profile_launch_times.append(atmospheric_profile.atmospheric_profile_launch_time)

        atmospheric_profile.get_launch_duration()
        self.profile_durations.append(atmospheric_profile.atmospheric_profile_duration)
        self.profile_final_times.append(atmospheric_profile.atmospheric_profile_final_time)

        atmospheric_profile.calculate_solar_longitude()
        self.sun_longitudes.append(atmospheric_profile.solar_longitude)

        atmospheric_profile.create_arts_atmospheric_profile()

        atmospheric_profile.calculate_atmospheric_fluxes()
        self.surface_down_solar_fluxes.append(atmospheric_profile.downwelling_sw_rad)
        self.surface_olr_fluxes.append(atmospheric_profile.surface_olr)

        self.profiles_altitude.append(atmospheric_profile.altitude)
        self.profiles_down_solar_flux.append(atmospheric_profile.solar.down)
        self.profiles_olr.append(atmospheric_profile.thermal.up)

    def populate(self):
        """Fills the lists with data."""

        for atmospheric_profile in self.atmospheric_profiles:
            self._extract_fluxes(atmospheric_profile=atmospheric_profile)


filepath = "C:\\Users\\drowe\\BCO_midday_clear_sky_soundings_ORCESTRA.nc"

profiles_list = []

profiles_num = list(range(0, 11, 1))

for profile_num in profiles_num:
    profile = AtmosphericProfile(filepath, sonde_index=profile_num)
    profiles_list.append(profile)

data = SurfaceFluxes(profiles_list)
data.populate()

print(
    f"radiosonde_launch_times = {data.profile_launch_times}\n\nradiosonde_final_times = {data.profile_final_times}\n\n"
    f"down_sw = {data.surface_down_solar_fluxes}\n\nolr = {data.surface_olr_fluxes}"
)

# Plot flux profiles of all radiosondes
fig1, ax1 = plt.subplots(figsize=(7, 10), dpi=400)

for i in range(0, len(data.profiles_down_solar_flux), 1):
    ax1.plot(data.profiles_down_solar_flux[i], data.profiles_altitude[i] / 1e3)

    ax1.set(
        title="Profiles of Global Downwelling Solar Radiation",
        xlabel="Downward Solar Flux W/m$^2$",
        ylabel="Altitude (km)".title()
    )

    ax1.set_xlim(850, 1350)

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_position(('outward', 0))
    ax1.spines['bottom'].set_position(('outward', 0))

fig2, ax2 = plt.subplots(figsize=(7, 10), dpi=400)

for i in range(0, len(data.profiles_olr), 1):
    ax2.plot(data.profiles_olr[i], data.profiles_altitude[i] / 1e3)

    ax2.set(
        title="Profiles of Outgoing Longwave Radiation",
        xlabel="Outgoing Thermal Flux W/m$^2$",
        ylabel="Altitude (km)".title()
    )

    ax2.set_xlim(225, 500)

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_position(('outward', 0))
    ax2.spines['bottom'].set_position(('outward', 0))

plt.show()

