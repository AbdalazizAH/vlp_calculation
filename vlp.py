import BeggsandBrill as BB


class Vlp:
    """
    Calculates the vertical lift performance (VLP) of a well.
    """

    def __init__(
        self,
        oil_rate=None,
        GOR=None,
        gas_grav=None,
        oil_grav=None,
        wtr_grav=None,
        dimameters=None,
        angle=None,
        thp=None,
        tht=None,
        twf=None,
        depth=None,
        sample_size=None,
        pressure=None,
        thickness=None,
        k=None,
        visc=None,
        GasGrav=None,
        themp=None,
        rw=None,
        re=None,
        s=None,
        oilfvf=None,
        water_rate=None,
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
        self.rates = self.linspace(0.1, 3992, 50)

    def linspace(self, start, stop, num):
        """
        Emulates numpy's linspace function.
        """
        step = (stop - start) / (num - 1)
        return [start + i * step for i in range(num)]

    def t_grad(self):
        """
        Calculates the temperature gradient.
        """

        if self.depth == 0:
            return 0
        else:
            return abs(self.tht - self.twf) / self.depth

    def depths(self):
        """
        Creates a list of depths at which to calculate pressure.
        """

        return self.linspace(0, self.depth, 51)

    def temps_gradient(self):
        """
        Calculates the temperature at each depth based on the gradient.
        """
        depths = self.depths()
        return [self.tht + self.t_grad() * depth for depth in depths]

    def pressure_traverse(self, QASume):
        """
        Calculates the pressure profile for a given flow rate.
        """

        p = [self.thp]  # Initialize with tubing head pressure
        dpdz = []

        depths = self.depths()
        temps_gradient = self.temps_gradient()

        for i in range(1, len(depths)):
            dz = depths[i] - depths[i - 1]
            dpdz_step = BB.Pgrad(
                p[-1],
                temps_gradient[i],
                QASume,
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
        ratess = []
        for QASume in self.rates:
            p, _ = self.pressure_traverse(QASume)
            bhp = p[-1]
            # print(bhp)
            bhps.append(bhp)
            ratess.append(QASume)
            vlp = [ratess, bhps]
        # return bhps
        return vlp
