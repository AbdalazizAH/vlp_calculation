import vlp as v


test = v.Vlp(
    oil_rate=100,
    water_rate=50,
    GOR=300,
    gas_grav=0.65,
    oil_grav=35,
    wtr_grav=1.07,
    dimameters=2.441,
    angle=90,
    thp=150,
    tht=100,
    twf=150,
    depth=5000,
    sample_size=51,
    pressure=4000.0,
    thickness=10,
    k=50.0,
    visc=0.5,
    GasGrav=0.65,
    themp=150.0,
    rw=0.328,
    re=1053.0,
    s=5,
    oilfvf=1.2,
)
vlps = test.vlp()

import matplotlib.pyplot as plt

plt.plot(vlps[0], vlps[1])
plt.show()
