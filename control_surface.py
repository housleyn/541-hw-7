class ControlSurface:
    def __init__(self):
        self.u = None
        self.orientation = None 
        self.area = None
        self.position = None
        self.b = 0.0
        self.aE = 0.0
        self.aW = 0.0
        self.aP = 0.0
        self.u_old = None
        self.d = 0.0
        self.u_prime = 0.0
        self.u_corr = 0.0
        self.Fe = 0.0
        self.Fw = 0.0
        

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        self._u = value
        self.orientation = 'vertical'  # u is aligned with x-direction → vertical face

  
    def calculate_x_coefficients(self, De, Fe, Dw, Fw,pe, pw,Aw, Ae):
        self.Aw = Aw
        self.Ae = Ae
        self.Fe = Fe
        self.Fw = Fw
        self.aE = De + max(-Fe, 0) 
        self.aW = Dw + max(Fw,0)
        self.aP = self.aE + self.aW + (Fe-Fw)  
        self.b = (pe-pw)*.5*(Aw+Ae) 
        self.d = self.area / self.aP
    def calculate_x_coefficients_first(self, De, Fe, Dw, Fw,pe, pw, Aa):
        self.Aa = Aa
        self.Fe = Fe
        self.Fw = Fw
        A1 = self.area
        self.aE = 0 
        self.aW = 0
        self.aP = self.aE + self.aW + Fw*(.5*(A1/Aa)**2) +Fe 
        self.b = (pw-pe)*A1 + Fw *(A1/Aa)*self.u_old 
        self.d = A1 / self.aP
        self.ua = A1/Aa * self.u_old
    def calculate_x_coefficients_last(self, De, Fe, Dw, Fw, pe, pw):
        self.Fe = Fe
        self.Fw = Fw
        A = self.area
        self.aE = De + max(-Fe, 0) 
        self.aW = Dw + max(Fw,0)
        self.aP = self.aE + self.aW + (Fe-Fw)  
        self.b = (pw-pe)*(A) 
        self.d = self.area / self.aP

    def calculate_area(self, domain):
        y_lower = domain.lower_boundary_func(self.position, 0)
        y_upper = domain.upper_boundary_func(self.position, 0)

        self.area = abs(y_upper - y_lower)
        