import math
import PyKEP as kep
import numpy as np
from collections import defaultdict

CUBE_RES = 10e3

#"""Returns a list of list of debris which are at time t in the same volume of resolution res."""
def cube(planets, res=CUBE_RES):
    cubemap = defaultdict(list)
    for p in planets:
        try:
            key = tuple(map(lambda foo: int(math.floor(foo/CUBE_RES)), p._r))
            cubemap[key].append(p)
        except:
            pass
    res = [foo for foo in cubemap.values() if len(foo) > 1]
    return res


def collision_prob(p1, p2):
    sigma = (p1.radius + p2.radius) ** 2 * math.pi # crossectional collision area
    dU = CUBE_RES ** 3 # volume of the cube
    Vimp = np.linalg.norm(p1._v - p2._v) # rel. velocity
    return  Vimp/dU * sigma


def setradii(planets, satcat):
    for p in planets:
        try:
            a = float(satcat[p.name.strip()].radarA)
            p.radius = math.sqrt(a/math.pi)
        except Exception as e:
            # print e
            p.radius = 0.5 # TODO: sample from radii distribution / use mean


def update(planets, ep):
    for p in planets[:]:
        try:
            p._r, p._v = map(np.array, p.eph(ep))
        except:
            # satelite decayed
            planets.remove(p)

