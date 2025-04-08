class Boundary():
    def __init__(self, mesh):
        self.mesh = mesh
        self.left_boundary = None
        self.right_boundary = None


    def apply_pressure_boundary(self, left_boundary, right_boundary):
        first_node = self.mesh.nodes[0]
        last_node = self.mesh.nodes[-1]
        self.left_boundary = left_boundary
        self.right_boundary = right_boundary
        first_node.p = left_boundary
        last_node.p = right_boundary

        first_node.aP = 1
        last_node.aP = 1
        first_node.aE = 0
        last_node.aW = 0
        first_node.aW = 0
        last_node.aE = 0
        first_node.b = left_boundary
        last_node.b = right_boundary


    