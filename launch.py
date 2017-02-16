import PyKEP as kep
import numpy as np
import random
import math
import datetime
from scipy.stats import norm
import matplotlib.mlab as mlab
from collections import namedtuple

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

importantAssets = []
importantAssetsTLE = [] 
inclination_distributions = []
eccentricity_distributions = []
W_distributions = []
w_distributions = []
monthOfLaunchesDistributions = []
yearLaunchesDistributions = []
meanMotion_distributions = []

# sample: 
# a (semi major axis = altitude + earth radius)
# i (inclination)
# W (progression of the orbital plane)
# w (where is the closest point of the elipse to the center?)
# always sample uniform between 0, 360 deg
# M (where is the debris on the orbit at the reference epoch?)
def loadData_satcat():
    satcat = kep.util.read_satcat('satcat2016.txt')
    return satcat

def loadData_tle():
    debris = kep.util.read_tle('tle2016.tle', with_name=False)
    return debris

def findImportantAssets(planets, currentYear, expectedLifespan):
    importantAssets = dict()
    for i in planets.items():
        beginning = currentYear - expectedLifespan
        # some satcat entries are missing information - e.g. orbital status code NEA - no elements available
        try:
            apogee = int(i[1].apogee)
            perigee = int(i[1].perigee)
            yearOfLaunchNo = int(i[1].launchdate.split('-')[0]) 
        except: 
            apogee = 0
            perigee = 0
            yearOfLaunchNo = 0
        if yearOfLaunchNo >= beginning and i[1].decay.strip() == "" and 'DEB' not in i[1].name and 'R/B' not in i[1].name and perigee > 0 and perigee < 2000:
            importantAssets[i[0]] = i[1]

    return importantAssets

def findTLEs(satcatList, TLEList):
    TLEs = []
    for planetTLE in TLEList:
        for planet in satcatList.items():
            if planet[0].strip() == planetTLE.name.strip():
                TLEs.append(planetTLE)
                break
    return TLEs 

def findDistributions(importantAssets, importantAssetsTLE, agentList):
    global inclination_distributions
    global eccentricity_distributions
    global W_distributions
    global w_distributions
    global monthOfLaunchesDistributions
    global yearLaunchesDistributions
    global meanMotion_distributions
    for agent in agentList:
        inclinPom = []
        eccentrPom = []
        WPom = []
        wPom = []
        meanMotionPom = []
        monthOfLaunches = []
        yearLaunches = []
        for plan in importantAssets.items():
            if plan[1][5].strip() in agent:
                for tleObj in importantAssetsTLE:
                    if plan[0].strip() == tleObj.name.strip():
                        inclinPom.append(float(tleObj.line2[8:16]))
                        # inclinPom.append(float(tleObj.line2[8:16]) * kep.DEG2RAD)
                        eccentrPom.append(float("0."+tleObj.line2[26:33]))
                        WPom.append(float(tleObj.line2[17:25]))
                        wPom.append(float(tleObj.line2[34:42]))
                        meanMotionPom.append(float(tleObj.line2[52:63]))
                        # WPom.append(float(tleObj.line2[17:25])* kep.DEG2RAD)
                        # wPom.append(float(tleObj.line2[34:42])* kep.DEG2RAD)
                        yearLaunches.append(int(plan[1][6][:4]))
                        monthOfLaunches.append(int(plan[1][6][5:7]))
                        break
#        print len(inclinPom)
        inclination_distributions.append(inclinPom)
        eccentricity_distributions.append(eccentrPom)
        W_distributions.append(WPom)
        w_distributions.append(wPom)
        meanMotion_distributions.append(meanMotionPom)
        monthOfLaunchesDistributions.append(monthOfLaunches)
        yearLaunchesDistributions.append(yearLaunches)

    # finding distributions for rest of the agents
    inclinPom = []
    eccentrPom = []
    WPom = []
    wPom = []
    meanMotionPom = []
    monthOfLaunches = []
    yearLaunches = []
    for plan in importantAssets.items():
        if plan[1][5].strip() not in agentList:
            for tleObj in importantAssetsTLE:
                if plan[0].strip() == tleObj.name.strip():
                    inclinPom.append(float(tleObj.line2[8:16]))
                    # inclinPom.append(float(tleObj.line2[8:16]) * kep.DEG2RAD)
                    eccentrPom.append(float("0."+tleObj.line2[26:33]))
                    WPom.append(float(tleObj.line2[17:25]))
                    wPom.append(float(tleObj.line2[34:42]))
                    meanMotionPom.append(float(tleObj.line2[52:63]))
                    # WPom.append(float(tleObj.line2[17:25])* kep.DEG2RAD)
                    # wPom.append(float(tleObj.line2[34:42])* kep.DEG2RAD)
                    yearLaunches.append(int(plan[1][6][:4]))
                    monthOfLaunches.append(int(plan[1][6][5:7]))
                    break
#    print len(inclinPom)
    inclination_distributions.append(inclinPom)
    eccentricity_distributions.append(eccentrPom)
    W_distributions.append(WPom)
    w_distributions.append(wPom)
    meanMotion_distributions.append(meanMotionPom)
    monthOfLaunchesDistributions.append(monthOfLaunches)
    yearLaunchesDistributions.append(yearLaunches)

def findHistograms(agentList):
    # lifespan for finding oe distributions
    expectedLifespan = 10
    planetsTLE = loadData_tle()
    satcat = loadData_satcat()
    importantAssets = findImportantAssets(satcat, 2016, expectedLifespan)
    importantAssetsTLE = findTLEs(importantAssets, planetsTLE)
    findDistributions(importantAssets, importantAssetsTLE, agentList)

    histValues_inclination = []
    histValues_eccentricity = []
    histValues_W = []
    histValues_w = []
    binsEdges_inclination = []
    binsEdges_eccentricity = []
    binsEdges_W = []
    binsEdges_w = []
    histValues_meanMotion = []
    binsEdges_meanMotion = []
    startOfLaunches = int(2016 - expectedLifespan)
    years = range(startOfLaunches, 2016, 1)    
    yearLaunches = []
    # iterate over players plus 1 for other players
    for i in range(0,len(agentList)+1):
        yearLaunchesPom = dict()
        for y1 in years:
            yearLaunchesPom[y1] = 0
        h_i, b_i = np.histogram(inclination_distributions[i], bins=50,normed=True)
        #h_i = h_i/sum(h_i)
        h_i = h_i.astype(np.float32)
        histValues_inclination.append(h_i)
        binsEdges_inclination.append(b_i)

        h_e, b_e = np.histogram(eccentricity_distributions[i], bins=50,normed=True)
        # h_e = h_e/sum(h_e)
        h_e = h_e.astype(np.float32)
        histValues_eccentricity.append(h_e)
        binsEdges_eccentricity.append(b_e)

        h_W, b_W = np.histogram(W_distributions[i], bins=50,normed=True)
        h_W = h_W.astype(np.float32)
        histValues_W.append(h_W)
        binsEdges_W.append(b_W)

        h_w, b_w = np.histogram(w_distributions[i], bins=50,normed=True)
        h_w = h_w.astype(np.float32)
        histValues_w.append(h_w)
        binsEdges_w.append(b_w)

        h_meanMotion, b_meanMotion = np.histogram(meanMotion_distributions[i], bins=50,normed=True)
        h_meanMotion = h_meanMotion.astype(np.float32)
        histValues_meanMotion.append(h_meanMotion)
        binsEdges_meanMotion.append(b_meanMotion)
#        yearLaunchesPom = [0] * len(years)
        for y in yearLaunchesDistributions[i]:
            if y in years:
                yearLaunchesPom[y] += 1
        yearLaunches.append(yearLaunchesPom)
#        print yearLaunchesPom
    oe_histograms = [histValues_inclination,binsEdges_inclination, histValues_eccentricity, binsEdges_eccentricity, histValues_W, binsEdges_W, histValues_w, binsEdges_w, histValues_meanMotion, binsEdges_meanMotion, yearLaunches]
    return oe_histograms    #print len(histValues_inclination)

def sample_oe(player, oe_histograms):
    histValues_inclination = oe_histograms[0]
    binsEdges_inclination = oe_histograms[1]
    histValues_eccentricity = oe_histograms[2]
    binsEdges_eccentricity = oe_histograms[3]
    histValues_W = oe_histograms[4]
    binsEdges_W = oe_histograms[5]
    histValues_w = oe_histograms[6]
    binsEdges_w = oe_histograms[7]
    histValues_meanMotion = oe_histograms[8]
    binsEdges_meanMotion = oe_histograms[9]
    # change player number to index of player
    player = player - 1
    probs_i = histValues_inclination[player]/sum(histValues_inclination[player])
    sampled_i_pom = int(np.random.choice(range(1,len(binsEdges_inclination[player])),1,p=probs_i))
    sampled_i = random.uniform(binsEdges_inclination[player][sampled_i_pom-1], binsEdges_inclination[player][sampled_i_pom])
    
    probs_e = histValues_eccentricity[player]/sum(histValues_eccentricity[player])
    sampled_e_pom = int(np.random.choice(range(1,len(binsEdges_eccentricity[player])),1,p=probs_e))
    sampled_e = random.uniform(binsEdges_eccentricity[player][sampled_e_pom-1], binsEdges_eccentricity[player][sampled_e_pom])
    
    probs_W = histValues_W[player]/sum(histValues_W[player])
    sampled_W_pom = int(np.random.choice(range(1,len(binsEdges_W[player])),1,p=probs_W))
    sampled_W = random.uniform(binsEdges_W[player][sampled_W_pom-1], binsEdges_W[player][sampled_W_pom])
    
    probs_w = histValues_w[player]/sum(histValues_w[player])
    sampled_w_pom = int(np.random.choice(range(1,len(binsEdges_w[player])),1,p=probs_w))
    sampled_w = random.uniform(binsEdges_w[player][sampled_w_pom-1], binsEdges_w[player][sampled_w_pom])
    
    probs_meanMotion = histValues_meanMotion[player]/sum(histValues_meanMotion[player])
    sampled_meanMotion_pom = int(np.random.choice(range(1,len(binsEdges_meanMotion[player])),1,p=probs_meanMotion))
    sampled_meanMotion = random.uniform(binsEdges_meanMotion[player][sampled_meanMotion_pom-1], binsEdges_meanMotion[player][sampled_meanMotion_pom])

    sample_M = random.uniform(0,360)
    
    return sampled_i, sampled_e, sampled_W, sampled_w, sample_M, sampled_meanMotion

def _checksum(line):
    res = 0
    for c in line:
        if 48 < ord(c) <= 58:
            res += ord(c) - 48
        if c == '-':
            res += 1
    return res % 10

def create_tle(oscul_elements, agentNo, year, month, launchIterator):
#    ep = int(math.floor(16*365.25))
    month = int(month)
    if len(str(month)) == 1:
        month = '0' + str(month)

    if len(str(launchIterator)) == 1:
        launchIterator = '00'+ str(launchIterator)
    if len(str(launchIterator)) == 2:
        launchIterator = '0'+ str(launchIterator)

    mu = float(398600)
    n = float(oscul_elements[5])
    a = (mu/(n * 2. * math.pi/(float(24 * 3600))) ** 2.) ** (1./3.)
    e = float(oscul_elements[1])
    apogee = (1+e)*(a - 6378)
    perigee = (1-e)*(a - 6378)

    # satellite name with code of agent
    satellite_number = str(year)[2:]+str(agentNo)+'99'
#    print satellite_number
    launchNo = str(month) + str(launchIterator)
    line1 = "1 "+satellite_number+"  "+str(year)[2:]+launchNo+"  11111.11111111 +.11111111 +00000-0 -11111-1 0  1110"
#    ep_date = (datetime.datetime(2000, 1, 1) + datetime.timedelta(ep)).timetuple()
    ep_date = str(year)
    ep_day = 99.9999999
    ep_str = str(ep_date[2:] + '{:12.8f}'.format(ep_day)[:14])
    # TODO change satellite ID and posibly alter international identifier?
    line1 = line1[:18] + ep_str + line1[32:-1]
    line1 += str(_checksum(line1))
    # print oscul_elements
    line2 = "2 " + satellite_number + " "
    line2 += '{:8.4f} '.format(oscul_elements[0]) # inclination (i)
    line2 += '{:8.4f} '.format(oscul_elements[2])  # RA (W)
    line2 += '{:.7f} '.format(oscul_elements[1])[2:]            # eccentrictiy (e)
    line2 += '{:8.4f} '.format(oscul_elements[3]) # argument of perigee (w)
    line2 += '{:8.4f} '.format(oscul_elements[4]) # mean anomaly (M)
    line2 += '{:11.8f}'.format(oscul_elements[5]) # mean motion (n)
    line2 += '{:5d}'.format(0) # revolutions
    line2 += str(_checksum(line2))
    tlePlanet = kep.planet.tle(line1, line2)
    tlePlanet.name = str(year)+'-'+launchNo
    return tlePlanet

def plotOrbits(planets_tle):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ep = kep.epoch(17)
    for pln in planets_tle:
        e = float("0."+ pln.line2[26:32])
        i = float(pln.line2[8:15])
        W = float(pln.line2[17:24])
        w = float(pln.line2[34:41])
        M = 20
        # print e,i,W,w
        oe = [2000000 + 6378000, e , i * kep.DEG2RAD, W * kep.DEG2RAD, w * kep.DEG2RAD, M * kep.DEG2RAD]
        pl = kep.planet.keplerian(ep, oe, earth.mu_self, 0, 0, 0, '')
        kep.orbit_plots.plot_planet(pl, ax=ax, alpha=0.2, s=0, color='red')
        # kep.orbit_plots.plot_planet(pl, ax=ax, alpha=0.2, color='red')
    plt.show()

def create_satcat(tleObject, agentNo, year, month):
    if len(str(month)) == 1:
        month = '0' + str(month)
    # Standard gravitational parameter for the earth    
    mu = float(398600)
    n = float(tleObject.line2[52:63])
    a = (mu/(n * 2. * math.pi/(float(24 * 3600))) ** 2.) ** (1./3.)
    e = float('0.'+tleObject.line2[26:33])

    agents = ['CIS','US','PRC','EU', 'OTH']
    nameOfNew = tleObject.name.strip()
    satcatentry = namedtuple('satcatentry', 'noradn multnameflag payloadflag operationstatus name ownership launchdate launchsite decay period incl apogee perigee radarA orbitstatus')
    intdsgn = tleObject.name.strip() 
    noradNo = '9999'
    nameFlag = ' '
    payloadFlag = ' '
    statusCode = ' '
    sateliteName = 'NEW_LAUNCH'
    owner = agents[agentNo-1]
    launchDate = str(year)+'-'+str(month)+'-'+'99'
    launchSite = 'none'
    decayDate = ' '
    orbitalPeriod = '9999'
    inclination = float(tleObject.line2[8:15]) 
    apogee = (1+e)*(a - 6378)
    perigee = (1-e)*(a - 6378)
    crossSectionA = '99'
    orbStatusCode = ' '
#    satcat[intdsgn] = satcatentry(noradNo,nameFlag,payloadFlag,statusCode,sateliteName,owner,launchDate,launchSite,decayDate,orbitalPeriod,inclination,apogee,perigee,crossSectionA,orbStatusCode)

    return [intdsgn,satcatentry(noradNo,nameFlag,payloadFlag,statusCode,sateliteName,owner,launchDate,launchSite,decayDate,orbitalPeriod,inclination,apogee,perigee,crossSectionA,orbStatusCode)]

def get_new_tle_and_satcat(agentNo, oe_histograms, year, month, launchIterator):
    tle = create_tle(sample_oe(agentNo, oe_histograms),agentNo, year, month, launchIterator)
    return tle, create_satcat(tle, agentNo, year, month)
    

def getExtrapolationLinearFunction(yearLaunchesDistribution, agentIndex):
    expectedLifespan = 10
    startOfLaunches = int(2016 - expectedLifespan)
    yearsRecent = range(startOfLaunches, 2016, 1)    
#    yearsRecent = range(1996,2016,1)
#    print len(yearsRecent), len(yearLaunchesDistribution[agentIndex].values())
#    print yearsRecent
#    print yearLaunchesDistribution[agentIndex].values()
    a,b = np.polyfit(yearsRecent,yearLaunchesDistribution[agentIndex].values(), 1)
    return a, b

#oe_histograms = []
#oe_histograms = findHistograms()

EU_list = ["EU","ASRA","BEL","CZCH","DEN","ESA","ESRO","EST","EUME","EUTE","FGER","FR","FRIT","GER","GREC","HUN","IT","LTU","LUXE","NETH","NOR","POL","POR","SPN","SWED","SWTZ","UK"]
agentList = ['CIS', 'US', 'PRC', EU_list]
#
oe_histograms = findHistograms(agentList)
##print inclination_distributions[3]
#weights0 = np.ones_like(inclination_distributions[0])/len(inclination_distributions[0])
#weights1 = np.ones_like(inclination_distributions[1])/len(inclination_distributions[1])
#weights2 = np.ones_like(inclination_distributions[2])/len(inclination_distributions[2])
#weights3 = np.ones_like(inclination_distributions[3])/len(inclination_distributions[3])
#
#n2, binsEdges2, patches2 = plt.hist(inclination_distributions[0], 50, weights=weights0, facecolor='yellow', alpha=0.5, label = 'RU')
#n2, binsEdges2, patches2 = plt.hist(inclination_distributions[1], 50, weights=weights1, facecolor='red', alpha=0.5, label = 'US')
#n2, binsEdges2, patches2 = plt.hist(inclination_distributions[2], 50, weights=weights2,  facecolor='green', alpha=0.5, label = 'CN')
#n2, binsEdges2, patches2 = plt.hist(inclination_distributions[3], 50, weights=weights3,  facecolor='blue', alpha=0.5, label = 'EU')
#
#plt.title('Inclination Distributions', fontsize=22)
#plt.legend(loc='upper left')
#plt.xlabel('degrees [DEG]', fontsize=22)
#plt.ylabel('normalized number of objects', fontsize=22)
#plt.rcParams.update({'font.size': 16})
#plt.show()

#oe_histograms = findHistograms()
#weights0 = np.ones_like(eccentricity_distributions[0])/len(eccentricity_distributions[0])
#weights1 = np.ones_like(eccentricity_distributions[1])/len(eccentricity_distributions[1])
#weights2 = np.ones_like(eccentricity_distributions[2])/len(eccentricity_distributions[2])
#weights3 = np.ones_like(eccentricity_distributions[3])/len(eccentricity_distributions[3])
##plt.hist(myarray, weights=weights)
#
#n2, binsEdges2, patches2 = plt.hist(eccentricity_distributions[0], 50, weights=weights0, facecolor='yellow', alpha=0.5, label = 'RU')
#n2, binsEdges2, patches2 = plt.hist(eccentricity_distributions[1], 50, weights=weights1, facecolor='red', alpha=0.5, label = 'US')
#n2, binsEdges2, patches2 = plt.hist(eccentricity_distributions[2], 50, weights=weights2,  facecolor='green', alpha=0.5, label = 'CN')
#n2, binsEdges2, patches2 = plt.hist(eccentricity_distributions[3], 50, weights=weights3,  facecolor='blue', alpha=0.5, label = 'EU')
#
#plt.title('Eccentricity Distributions', fontsize=22)
#plt.legend(loc='upper right')
#plt.xlabel('Eccentricity ratio', fontsize=22)
#plt.ylabel('normalized number of objects', fontsize=22)
#plt.rcParams.update({'font.size': 16})
#plt.show()


#oe_histograms = findHistograms()
#weights0 = np.ones_like(W_distributions[0])/len(W_distributions[0])
#weights1 = np.ones_like(W_distributions[1])/len(W_distributions[1])
#weights2 = np.ones_like(W_distributions[2])/len(W_distributions[2])
#weights3 = np.ones_like(W_distributions[3])/len(W_distributions[3])
#
#n2, binsEdges2, patches2 = plt.hist(W_distributions[0], 50, weights=weights0, facecolor='yellow', alpha=0.5, label = 'RU')
#n2, binsEdges2, patches2 = plt.hist(W_distributions[1], 50, weights=weights1, facecolor='red', alpha=0.5, label = 'US')
#n2, binsEdges2, patches2 = plt.hist(W_distributions[2], 50, weights=weights2,  facecolor='green', alpha=0.5, label = 'CN')
#n2, binsEdges2, patches2 = plt.hist(W_distributions[3], 50, weights=weights3,  facecolor='blue', alpha=0.5, label = 'EU')
#
#plt.title('Longitude Of The Ascending Node Distributions', fontsize=20)
#plt.legend(loc='upper left')
#plt.xlabel('degrees [DEG]', fontsize=22)
#plt.ylabel('normalized number of objects', fontsize=22)
#plt.rcParams.update({'font.size': 16})
#plt.show()

#oe_histograms = findHistograms()
#
#weights0 = np.ones_like(w_distributions[0])/len(w_distributions[0])
#weights1 = np.ones_like(w_distributions[1])/len(w_distributions[1])
#weights2 = np.ones_like(w_distributions[2])/len(w_distributions[2])
#weights3 = np.ones_like(w_distributions[3])/len(w_distributions[3])
##
#n2, binsEdges2, patches2 = plt.hist(w_distributions[0], 50, weights=weights0, facecolor='yellow', alpha=0.5, label = 'RU')
#n2, binsEdges2, patches2 = plt.hist(w_distributions[1], 50, weights=weights1, facecolor='red', alpha=0.5, label = 'US')
#n2, binsEdges2, patches2 = plt.hist(w_distributions[2], 50, weights=weights2,  facecolor='green', alpha=0.5, label = 'CN')
#n2, binsEdges2, patches2 = plt.hist(w_distributions[3], 50, weights=weights3,  facecolor='blue', alpha=0.5, label = 'EU')
##
#plt.title('Argument Of Periapsis Distributions', fontsize=22)
#plt.legend(loc='upper left')
#plt.xlabel('degrees [DEG]', fontsize=22)
#plt.ylabel('normalized number of objects', fontsize=22)
#plt.rcParams.update({'font.size': 16})
#plt.show()
#print importantAssets.items()[100][1][6][5:7]
