import numpy as np
import BeggsandBrill as BB

class Vlp:
    """
    Calculates the vertical lift performance (VLP) of a well.
    """

    def __init__(
        self,
        oil_rate: float=None,
        GOR:      float=None,
        gas_grav: float=None,
        oil_grav: float=None,
        wtr_grav: float=None,
        dimameters: float=None,
        angle: float=None,
        thp: float=None,
        tht:float=None,
        twf:float=None,
        depth:float=None,
        sample_size:float=None,
        pressure:float=None,
        thickness:float=None,
        k:float=None,
        visc:float=None,
        GasGrav:float=None,
        themp:float=None,
        rw:float=None,
        re:float=None,
        s:float=None,
        oilfvf:float=None,
        water_rate:float=None,
    ):
        """
        Initializes the Vlp class with the given parameters.
        """

        self.oil_rate = oil_rate
        self.water_rate = water_rate
        self.GOR = GOR
        self.gas_grav = gas_grav
        self.oil_grav = oil_grav
        self.wtr_grav = wtr_grav
        self.dimameters = dimameters
        self.angle = angle
        self.thp = thp
        self.tht = tht
        self.twf = twf
        self.depth = depth
        self.sample_size = sample_size

        # Additional parameters (potentially unused?)
        self.pressure = pressure
        self.thickness = thickness
        self.k = k
        self.visc = visc
        self.GasGrav = GasGrav
        self.themp = themp
        self.rw = rw
        self.re = re
        self.s = s
        self.oilfvf = oilfvf

        # Define range of flow rates to evaluate (ASSYM #)
        self.rates = np.linspace(0.1, 3992, 50)

    def t_grad(self):
        """
        Calculates the temperature gradient.
        """

        if self.depth == 0:
            return 0
        else:
            return abs(self.tht - self.twf) / self.depth

    def depths(self) -> list:
        """
        Creates a list of depths at which to calculate pressure.
        """

        return np.linspace(0, self.depth, 51)

    def temps_gradient(self) -> list:
        """
        Calculates the temperature at each depth based on the gradient.
        """

        return self.tht + self.t_grad() * self.depths()

    def pressure_traverse(self, q):
        """
        Calculates the pressure profile for a given flow rate.
        """

        p = [self.thp]  # Initialize with tubing head pressure
        dpdz = []

        for i in range(1, len(self.depths())):
            dz = self.depths()[i] - self.depths()[i - 1]
            dpdz_step = BB.Pgrad(
                p[-1],
                self.temps_gradient()[i],
                q,
                self.water_rate,
                self.GOR,
                self.GasGrav,
                self.oil_grav,
                self.wtr_grav,
                self.dimameters,
                self.angle,
            )
            dpdz.append(dpdz_step)
            pressure = p[-1] + dz * dpdz_step
            p.append(pressure)

        return p, dpdz

    def vlp(self):
        """
        Calculates the VLP curve (bottomhole pressure vs. flow rate).
        """

        bhps = []
        ratess =[]
        for q in self.rates:
            p, _ = self.pressure_traverse(q)
            bhp = p[-1]
            print(bhp)
            bhps.append(bhp)
            ratess.append(q)
            vlp = [ratess,bhps ]
        # return bhps
        return vlp
