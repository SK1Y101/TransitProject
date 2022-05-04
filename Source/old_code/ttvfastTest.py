from ttvfast import models
import ttvfast

gravity = 0.000295994511                        # AU^3/day^2/M_sun
stellar_mass = 0.95573417954                    # M_sun

planet1 = models.Planet(
    mass=0.00002878248,                         # M_sun
    period=1.0917340278625494e+01,              # days
    eccentricity=5.6159310042858110e-02,
    inclination=9.0921164935951211e+01,         # degrees
    longnode=-1.1729336712101943e-18,           # degrees
    argument=1.8094838714599581e+02,            # degrees
    mean_anomaly=-8.7093652691581923e+01,       # degrees
)

planet2 = models.Planet(
    mass=0.00061895914,
    period=2.2266898036209028e+01,
    eccentricity=5.6691301931178648e-02,
    inclination=8.7598285693573246e+01,
    longnode=4.6220554014026838e-01,
    argument=1.6437004273382669e+00,
    mean_anomaly=-1.9584857031843157e+01,
)

planets = [planet1, planet2]
Time = -1045                                    # days
dt = 0.54                                       # days
Total = 1700                                    # days

results = ttvfast.ttvfast(planets, stellar_mass, Time, dt, Total)

import numpy as np

print(np.array(results["positions"]))
