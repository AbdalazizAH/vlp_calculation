import pandas as pd
from matplotlib import pyplot as plt
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
vlps =test.vlp()


vlp_data = pd.DataFrame({"Flow Rate": vlps[0], "Bottomhole Pressure": vlps[1]})
print (vlp_data)
# Plotting
plt.plot(vlp_data["Flow Rate"], vlp_data["Bottomhole Pressure"])
plt.xlabel("Flow Rate")
plt.ylabel("Bottomhole Pressure")
plt.title("VLP Curve")
plt.show()
