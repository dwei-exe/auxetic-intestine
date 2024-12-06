import math

class Parametric2D:
    def __init__(self, x, y, z, lowS, highS, lowT, highT):
        self.x = x
        self.y = y
        self.z = z
        self.lowS = lowS
        self.highS = highS
        self.lowT = lowT
        self.highT = highT
    def curvature(self, resolutionS, resolutionT):
        for i in range(resolutionS + 1):
            s = lowS + i * (highS - lowS) / resolutionS
            for j in range(resolutionT + 1):
                t = lowT + j * (highT - lowT) / resolutionT
    def derivatives(self, resolutionS, resolutionT):
        for i in range(resolutionS + 1):
            s = lowS + i * (highS - lowS) / resolutionS
            for j in range(resolutionT + 1):
                t = lowT + j * (highT - lowT) / resolutionT

colon = Parametric2D(
    lambda s, t: (s + math.pi) / 3,
    lambda s, t: math.cos(t) * (1 + 1 / 5 * sin(s)),
    lambda s, t: math.sin(t) * (1 + 1 / 5 * sin(s)),
    -20, 20,
    0, 2 * math.pi
)
