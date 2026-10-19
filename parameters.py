import datetime
import numpy as np
from typing import Final

START_DATE: Final = datetime.datetime(2026, 10, 11) # arbitrary start date at 12:00 am

MAX_SIM_HOURS: Final = 1
DT_CONTROLLER: Final = 0.25     # sampling time
DUTY_CYCLE = 0.6                # expressed as decimal ratio of sampling time
DT_SIMULATOR: Final = 0.05      # physics engine update frequency

INITIAL_ORIENTATION = [1, 0, 0, 0]              # scalar, q1, q2, q3
INITIAL_ANGULAR_VELOCITY = np.array([1.5, 2.5, 1])  # x, y, z
INERTIA_MATRIX = np.array([[0.0022, 0, 0],           # in satellite body frame - must be diagonalised
                           [0, 0.0022, 0],
                           [0, 0, 0.0022]])
INITIAL_KEPLER_ELEMENTS = {
    "SEMI": 6700,   # semimajor_axis; unit: km
    "INCL": 97,     # inclination; unit: degrees
    "RAAN": 310,    # right ascension of the ascending node; unit: degrees
    "ARGP": 0,      # argument_of_periapsis; unit: degrees
    "ANOM": 0,      # true_anomaly; unit: degrees
    "ECCT": 0       # eccentricity; ratio - assumed 0
    }
MAGNETOMETER_HARD_IRON_BIAS = np.array([5, 2, 3])  # microtesla
MAX__MAGNETORQUER_MOMENT: Final = np.array([0.2, 0.2, 0.2])  # A m^2
