import numpy as np
import matplotlib.pyplot as plt

class SIMPLE():
    def __init__(self, mesh, boundary, material):
        self.mesh = mesh
        self.boundary = boundary
        self.material = material
        self.domain = mesh.domain
        self.rho = material.rho
        self.mu = material.mu
        self.nx = mesh.nx
        self.alphap = .8
        self.alphau = .8
        self.iteration = 0
        self.tol = 1e-1
        self.runtime = 0.0

    def generate_initial_guesses(self):
        mdot = 1

        for i in range(self.nx-1):
            self.mesh.u_faces[i].u_old = mdot / (self.rho * self.mesh.u_faces[i].area)
        for i in range(self.nx):
            self.mesh.nodes[i].p = self.boundary.left_boundary - self.mesh.nodes[i].position * (self.boundary.left_boundary - self.boundary.right_boundary) / (self.mesh.x_max - self.mesh.x_min)


    def calculate_u_new(self):
        De = 0 
        Dw = 0
        self.mesh.u_faces[0].calculate_x_coefficients_first(De, self.rho*(self.mesh.u_faces[0].u_old + self.mesh.u_faces[1].u_old)*.5 * self.mesh.nodes[1].area, 
                                                      Dw,self.rho*self.mesh.u_faces[0].u_old*self.mesh.u_faces[0].area, 
                                                      self.mesh.nodes[1].p, self.mesh.nodes[0].p, self.mesh.nodes[0].area)
        for i in range(1,self.nx-2):
            pe = self.mesh.nodes[i].p
            pw = self.mesh.nodes[i+1].p
            Aw = self.mesh.nodes[i].area
            Ae = self.mesh.nodes[i+1].area
            Fw = .5 * self.rho * (self.mesh.u_faces[i].u_old + self.mesh.u_faces[i-1].u_old) * self.mesh.nodes[i].area
            Fe = .5 * self.rho * (self.mesh.u_faces[i].u_old + self.mesh.u_faces[i+1].u_old) * self.mesh.nodes[i+1].area
            self.mesh.u_faces[i].calculate_x_coefficients(De, Fe, Dw, Fw, pe, pw, Aw, Ae)
        self.mesh.u_faces[-1].calculate_x_coefficients_last(De, self.rho*self.mesh.u_faces[-1].u_old*self.mesh.u_faces[-1].area, 
                                       Dw, self.rho*(self.mesh.u_faces[-1].u_old + self.mesh.u_faces[-2].u_old)*.5 * self.mesh.nodes[-2].area, 
                                       self.mesh.nodes[-1].p, self.mesh.nodes[-2].p)
        A, b = self.mesh.construct_A_matrix_from_control_surfaces()
        u_solutions = np.linalg.solve(A, b)
        for i, u_value in enumerate(u_solutions):
            self.mesh.u_faces[i].u = u_value
        return A, b

    def calculate_p_prime(self):
        self.mesh.nodes[0].b = 0
        # self.mesh.nodes[-1].b = 0
        for i in range(1,self.nx-1):
            Fw = self.rho * self.mesh.u_faces[i-1].u * self.mesh.u_faces[i-1].area
            Fe = self.rho * self.mesh.u_faces[i].u * self.mesh.u_faces[i].area
            Aw = self.mesh.u_faces[i-1].area
            Ae = self.mesh.u_faces[i].area
            dw = self.mesh.u_faces[i-1].d
            de = self.mesh.u_faces[i].d
            self.mesh.nodes[i].calculate_coefficients(self.rho, Fw, Fe, Aw, Ae, dw, de)

        A, b = self.mesh.construct_A_matrix_from_nodes()
        p_prime_solutions = np.linalg.solve(A, b)
        for i, p_prime_value in enumerate(p_prime_solutions):
            self.mesh.nodes[i].p_prime = p_prime_value
        return A, b


    def solve_for_new_u_prime(self):
        for i in range(self.nx-1):
            self.mesh.u_faces[i].u_prime = self.mesh.u_faces[i].d *(self.mesh.nodes[i].p_prime - self.mesh.nodes[i+1].p_prime)
    
    def solve_for_corrected_velocities_and_pressures(self):
        for i in range(self.nx-1):
            self.mesh.u_faces[i].u_corr = self.mesh.u_faces[i].u + self.mesh.u_faces[i].u_prime
        self.mesh.nodes[0].p_corr = self.mesh.nodes[0].p - .5 * self.rho* (self.mesh.u_faces[0].u_corr * self.mesh.u_faces[0].area / self.mesh.nodes[0].area)**2 
        for i in range(1,self.nx):
            self.mesh.nodes[i].p_corr = self.mesh.nodes[i].p + self.mesh.nodes[i].p_prime
    
    def calculate_new_pressures_and_velocities(self):
        for i in range(self.nx-1):
            # self.mesh.u_faces[i].u_old = self.mesh.u_faces[i].u
            self.mesh.u_faces[i].u = (1-self.alphau)*self.mesh.u_faces[i].u_old + self.alphau*self.mesh.u_faces[i].u_corr
        for i in range(self.nx):
            self.mesh.nodes[i].p_old = self.mesh.nodes[i].p
            self.mesh.nodes[i].p = (1-self.alphap)*self.mesh.nodes[i].p_old + self.alphap*self.mesh.nodes[i].p_corr


    def set_old_to_new(self):
        for i in range(self.nx-1):
            self.mesh.u_faces[i].u_old = self.mesh.u_faces[i].u


    def run(self, runtime):
        if self.iteration == 0:
            self.generate_initial_guesses()
        converged = False
        for iteration in range(runtime):
            self.calculate_u_new()
            self.calculate_p_prime()
            self.solve_for_new_u_prime()
            self.solve_for_corrected_velocities_and_pressures()
            self.calculate_new_pressures_and_velocities()
            self.set_old_to_new()
            self.iteration += 1
            print(f"Iteration {self.iteration}:")
            if any(abs(self.mesh.u_faces[i].u - self.mesh.u_faces[i].u_old) < self.tol for i in range(self.nx-1)) and \
               any(abs(self.mesh.nodes[i].p - self.mesh.nodes[i].p_old) < self.tol for i in range(self.nx)) :
                converged = True
        # self.plot_pressure_vs_distance()

    def plot_pressure_vs_distance(self):
        # Extract the position and pressure values from the nodes
        
        positions = np.array([node.position for node in self.mesh.nodes])
        pressures = np.array([node.p for node in self.mesh.nodes])

        # Plot pressure vs axial distance
        plt.plot(positions, pressures, marker='o', linestyle='-', color='b')
        plt.xlabel('Axial Distance (x)')
        plt.ylabel('Pressure (p)')
        plt.title('Pressure vs Axial Distance')
        plt.grid(True)
        plt.show()