from node import Node
from control_surface import ControlSurface
import numpy as np 


class Mesh():
    def __init__(self, domain, nx):
        self.domain = domain
        self.nx = nx
        self.shape = nx
        self.nodes = np.empty(self.shape,dtype=object)
        self.u_faces = np.empty(self.shape-1, dtype=object) 

        x_bounds, y_bounds = self.domain.get_bounds()
        x_min, x_max = x_bounds

        self.x_min = x_min
        self.x_max = x_max

    def construct_nodes(self):
        for i in range(self.nx):
            x = self.domain.x_min + i * (self.domain.x_max - self.domain.x_min) / (self.nx-1)
            self.nodes[i] = Node()
            self.nodes[i].position = x
            self.nodes[i].calculate_area(self.domain)

    def construct_control_surfaces(self):
        for i in range(self.u_faces.shape[0]):
            self.u_faces[i] = ControlSurface()
            self.u_faces[i].position = (self.nodes[i].position + self.nodes[i+1].position) / 2
            self.u_faces[i].calculate_area(self.domain)
            

    def construct_mesh(self):
        self.construct_nodes()
        self.construct_control_surfaces()


    def construct_A_matrix_from_control_surfaces(self):
        N = self.u_faces.shape[0]
        A = np.zeros((N, N))
        b = np.zeros(N)

        for i, face in enumerate(self.u_faces):
           
            aP = face.aP
            aE = face.aE
            aW = face.aW 
            b[i] = face.b
            
            if i > 0:
                A[i, i-1] = -aW
            A[i, i] = aP
            if i < N - 1:
                A[i, i+1] = -aE
        
        return A, b

    def construct_A_matrix_from_nodes(self):
        N = self.nodes.shape[0]
        A = np.zeros((N, N))
        b = np.zeros(N)

        for i, node in enumerate(self.nodes):
            
            aP = node.aP 
            aE = node.aE 
            aW = node.aW 
            b[i] = node.b 
            
            if i > 0:
                A[i, i-1] = -aW
            A[i, i] = aP
            if i < N - 1:
                A[i, i+1] = -aE
        
        return A, b
