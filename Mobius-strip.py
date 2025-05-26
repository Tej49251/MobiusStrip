import numpy as np
from scipy.integrate import simps
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class MobiusStrip:
    def __init__(self, R=1.0, w=0.2, n=100):
        self.R = R          # Radius from center
        self.w = w          # Width of strip
        self.n = n          # Resolution (grid size)
        self.u = np.linspace(0, 2 * np.pi, self.n)
        self.v = np.linspace(-w / 2, w / 2, self.n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        self.X, self.Y, self.Z = self._generate_mesh()

    def _generate_mesh(self):
        """Generate the 3D coordinates from parametric equations"""
        U, V = self.U, self.V
        R = self.R

        X = (R + V * np.cos(U / 2)) * np.cos(U)
        Y = (R + V * np.cos(U / 2)) * np.sin(U)
        Z = V * np.sin(U / 2)
        return X, Y, Z

    def compute_surface_area(self):
        """Approximate surface area using numerical integration."""
        du = 2 * np.pi / (self.n - 1)
        dv = self.w / (self.n - 1)

        # Partial derivatives w.r.t. u and v
        Xu = np.gradient(self.X, du, axis=1)
        Yu = np.gradient(self.Y, du, axis=1)
        Zu = np.gradient(self.Z, du, axis=1)

        Xv = np.gradient(self.X, dv, axis=0)
        Yv = np.gradient(self.Y, dv, axis=0)
        Zv = np.gradient(self.Z, dv, axis=0)

        # Cross product of partials
        cross = np.cross(np.stack((Xu, Yu, Zu), axis=-1),
                         np.stack((Xv, Yv, Zv), axis=-1))

        # Norm of cross product
        dA = np.linalg.norm(cross, axis=-1)

        # Surface area using Simpson's rule
        area = simps(simps(dA, self.v), self.u)
        return area

    def compute_edge_length(self):
        """Approximate the length of the edge along v=+w/2"""
        u = self.u
        v_edge = self.w / 2
        x = (self.R + v_edge * np.cos(u / 2)) * np.cos(u)
        y = (self.R + v_edge * np.cos(u / 2)) * np.sin(u)
        z = v_edge * np.sin(u / 2)

        points = np.stack((x, y, z), axis=-1)
        length = np.sum([euclidean(points[i], points[i+1])
                         for i in range(len(points) - 1)])
        return length

    def plot(self):
        """Plot the Möbius strip using matplotlib"""
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.X, self.Y, self.Z, cmap='plasma', edgecolor='k', linewidth=0.2, alpha=0.9)
        ax.set_title("Möbius Strip")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    mobius = MobiusStrip(R=1.0, w=0.2, n=200)
    mobius.plot()

    area = mobius.compute_surface_area()
    edge_length = mobius.compute_edge_length()

    print(f"Surface Area ≈ {area:.4f}")
    print(f"Edge Length ≈ {edge_length:.4f}")
