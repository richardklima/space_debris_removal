import PyKEP as kep
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import datetime

# from cube import setradii, update


LB = 0.05
# LB = 0.11


def setmass(planets):    
    for p in planets:
        p.mass = 4./3. * math.pi * p.radius ** 3 * 92.937 * (p.radius * 2) ** -0.74 # EQ 1


def _Am(d):
    # Eq (6)
    logd = math.log10(d)
    # alpha
    if logd <= -1.95:
        alpha = 0.
    elif -1.95 < logd < 0.55:
        alpha = 0.3 + 0.4 * (logd + 1.2)
    else: # >= 0.55
        alpha = 1.
    # mu1
    if logd <= -1.1:
        mu1 = -0.6
    elif -1.1 < logd < 0:
        mu1 = -0.6 - 0.318 * (logd + 1.1)
    else: # >= 0
        mu1 = -0.95
    # sigma1
    if logd <= -1.3:
        sigma1 = 0.1
    elif -1.3 < logd < -0.3:
        sigma1 = 0.1 + 0.2 * (logd + 1.3)
    else:
        sigma1 = 0.3
    # mu2
    if logd <= -0.7:
        mu2 = -1.2
    elif -0.7 < logd < -0.1:
        mu2 = -1.2 - 1.333 * (logd + 0.7)
    else:
        mu2 = -2.0

    # sigma2
    if logd <= -0.5:
        sigma2 = 0.5
    elif -0.5 < logd < -0.3:
        sigma2 = 0.5 - (logd + 0.5)
    else:
        sigma2 = 0.3

    N1 = np.random.normal(mu1, sigma1, size=(1,1))
    N2 = np.random.normal(mu2, sigma1, size=(1,1))
    return float(10. ** (alpha * N1 + (1-alpha) * N2))


def _dv(Am):
    # Eq (12)
    mu = 0.9 * math.log10(Am) + 2.9
    sigma = 0.4
    N = np.random.normal(mu, sigma, size=(1,1))
    return float(10 ** N)



class _rv_dist(stats.rv_continuous):
    def _cdf(self, x):
        return np.where(x < LB, 0., 1 - x ** -1.71 / (LB ** -1.71))



def _checksum(line):
    res = 0
    for c in line:
        if 48 < ord(c) <= 58:  
            res += ord(c) - 48
        if c == '-':
            res += 1
    return res % 10


def _create_tles(ep, source, fragments):
    res = []

    ep_date = (datetime.datetime(2000, 1, 1) + datetime.timedelta(ep)).timetuple()
    ep_day = ep_date.tm_yday + ep_date.tm_hour/24. + ep_date.tm_min/24./60. + ep_date.tm_sec/24./60./60.
    ep_str = (str(ep_date.tm_year)[2:] + '{:12.8f}'.format(ep_day))[:14]
    # TODO change satellite ID and posibly alter international identifier?
    line1 = source.line1[:18] + ep_str + source.line1[32:-1]
    line1 += str(_checksum(line1))

    for fragment in fragments:
        r, v = source.eph(ep)
        v = np.array(v) + fragment[-3:] # add dV
        el = kep.ic2par(r, v, source.mu_central_body)
      
        try:    
            n = math.sqrt(source.mu_central_body/(4.*math.pi**2*el[0]**3)) * kep.DAY2SEC # mean motion in days 
        except:
            continue  
            # continue to next fragment 


        M = el[5] - el[1] * math.sin(el[5])
        if M < 0:
            M += 2 * math.pi

        line2 = source.line2[:8]
        line2 += '{:8.4f} '.format(el[2] * kep.RAD2DEG) # inclination (i)
        line2 += '{:8.4f} '.format(el[3] * kep.RAD2DEG) # RA (W)
        line2 += '{:.7f} '.format(el[1])[2:]            # eccentrictiy (e)
        line2 += '{:8.4f} '.format(el[4] * kep.RAD2DEG) # argument of perigee (w)
        line2 += '{:8.4f} '.format(M * kep.RAD2DEG) # mean anomaly (M)
        line2 += '{:11.8f}'.format(n)                  # mean motion

#        line2 += source.line2[63:68] #'{:5d}'.format(1) # revolutions
        line2 += '{:5d}'.format(0) # revolutions
        line2 += str(_checksum(line2))

        # sometimes there is an error 'Eccentricity out of range'
        try:
            res.append(kep.planet.tle(line1, line2))
        except:
            pass

    return res


def breakup(ep, p1, p2):
    # we assume d1 > d2
    if p1.radius < p2.radius:
        p1, p2 = p2, p1
    
    dv = np.linalg.norm(p1._v - p2._v)
    catastrophRatio = (p2.mass*dv**2)/(2*p1.mass*1000)

    # if the ratio is smaller than 40 J/g then it is non-catastrophic collision
    if catastrophRatio<40:
        M = p2.mass*dv
        Lc = np.linspace(LB, 1.)
        num = int(0.1 * M ** 0.75 * LB ** -1.71)
    else:
        # catastrophic collision
        M = p1.mass + p2.mass
        # Eq (2)
        Lc = np.linspace(LB, 1.)
        num = int(0.1 * M ** 0.75 * LB ** -1.71)

    # # Eq (2)
    # Lc = np.linspace(LB, 1.)

    # num = int(0.1 * M ** 0.75 * LB ** -1.71)


        
    ddist = _rv_dist()

    d = ddist.rvs(size=int(num+0.2*num)) # sample 20% more
    d = d[d < 2. * p1.radius][:num] # filter by max size

    # plt.hist(d, bins=100, color='w')
    # plt.axvline(max(d1, d2))
    # plt.xlabel('characteristic length [m]')
    # plt.show()

    
    A = 0.556945 * d ** 2.0047077
    Am = np.array(map(_Am, d))
    dv = np.array(map(_dv, Am))
    m = A / Am

    # create samples on unit sphere
    u = np.random.random(size=(len(dv),)) * 2. -1.
    theta = np.random.random(size=(len(dv),)) * 2. * np.pi
    v = np.sqrt(1 - u**2)
    p = np.array(zip(v * np.cos(theta), v * np.sin(theta), u))
    # dv velocity vectors
    try:
        dv_vec = np.array(zip(p[:,0] * dv, p[:,1] * dv, p[:,2] * dv))
    except:
        # if len(p.shape) != 2:
        return
        #     p = p[None,:]
        #     dv_vec = np.array(zip(p[:,0] * dv, p[:,1] * dv, p[:,2] * dv))
        # else:
        #     pass

#    import matplotlib.pyplot as plt
#    from mpl_toolkits.mplot3d import Axes3D
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.scatter(dv_vec[:,0], dv_vec[:,1], zs=dv_vec[:,2], s=d*100.)
#    plt.show()
#    
    #plt.plot(Lc, d.cdf(Lc), color='red')
#    plt.plot(Lc, d.pdf(Lc), color='blue')
#    plt.plot(Lc, pdf, color='green')

    
#    # Figure 5
#    fig = plt.figure()
#    ax = plt.gca()
#    ax.scatter(Am, d, marker='x', color='k')
#    ax.set_yscale('log')
#    ax.set_xscale('log')
#    plt.ylabel('A/m [m^2/kg]')
#    plt.xlabel('characteristic length [m]')
#    plt.show()
#
#
#    # Figure 6
#    plt.hist(Am, bins = 10 ** np.linspace(np.log10(0.01), np.log10(10), 100), color='w')
#    plt.gca().set_xscale('log')
#    plt.xlabel('A/m [m^2/kg]')
#    plt.show()
#
#
#    # Figure 7
#    plt.hist(dv, bins=100, color='w')
#    plt.xlabel('delta V [m/s]')
#    plt.show()
#
#    # is mass conserved? - short answer: no
#    print M, sum(m)
#
#    plt.hist(m, bins=100, color='w')
#    plt.xlabel('mass [kg]')
#    plt.show()
# 
    fragments = np.array(zip(d, A, Am, m, dv, dv_vec[:,0], dv_vec[:,1], dv_vec[:,2]))

    # let's at least not add mass to the system
    # print len(fragments)
    idx = np.argmin(np.cumsum(m) < M)
    if idx:
        fragments = fragments[:idx]
    # distribute according to size first, then mass
    _foo = fragments[:,0] > p2.radius * 2 # fragments too large to originate from satelite 2
    fragments1 = fragments[_foo]
    _m1 = sum(fragments1[:,3]) 
    _total_m = sum(fragments[:,3])

    # if the second planet is too small and there are no debris from it
    if len(fragments[_foo==False]) == 0:
        fragments2 = []
        # print 'EMPTY'
    else:
        _rest = np.random.permutation(fragments[_foo == False])
        _cumsum = np.cumsum(_rest[:,3]) # calculate the cumulative sum of the mass of the rest
        _idx = np.argmax(_cumsum > (p1.mass/(M) * _total_m) - _m1)
        fragments1 = np.vstack((fragments1, _rest[:_idx]))
        fragments2 = _rest[_idx:]

    # # TODO: match could be improved
    # _m1 = sum(fragments1[:,3])
    # _m2 = sum(fragments2[:,3])
    # print _m1/_m2, m1/m2, len(fragments1), len(fragments2)
    
    debris1 = _create_tles(ep, p1, fragments1)
    debris2 = _create_tles(ep, p2, fragments2)
    return debris1, debris2


def plot_orbits(ep, debris1, debris2):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    ax.set_xlim(-1e7, 1e7)
    ax.set_ylim(-1e7, 1e7)
    ax.set_zlim(-1e7, 1e7)
    for deb in debris1:
        try:
            kep.orbit_plots.plot_planet(deb, ax=ax, t0=kep.epoch(ep), s=0, color='r', alpha=0.2)
        except:
            pass

    for deb in debris2:
        try:
            kep.orbit_plots.plot_planet(deb, ax=ax, t0=kep.epoch(ep), s=0, color='b', alpha=0.2)
        except:
            pass

    plt.show()
