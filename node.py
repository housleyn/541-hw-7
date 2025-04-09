class Node():
    def __init__(self):
        self.aW = 0.0
        self.aE = 0.0
        self.aP = 0.0
        self.b = 0.0
        self.p = 0.0
        self.p_prime = 0.0
        self.p_new = 0.0
        self.position = None 
        self.area = 0.0
        self.p_corr = 0.0
        self.p_old = 0.0
        self.area = 0.0

    def calculate_coefficients(self,rho,Fw, Fe, Aw, Ae, dw, de):
        self.Fe = Fe
        self.Fw = Fw
        self.aE = rho * Ae *de
        self.aW = rho * Aw *dw


        self.aP = self.aE + self.aW 

        
        self.b = Fw - Fe

    def calculate_area(self, domain):
        y_lower = domain.lower_boundary_func(self.position, 0)
        y_upper = domain.upper_boundary_func(self.position, 0)

        self.area = abs(y_upper - y_lower)
            
            