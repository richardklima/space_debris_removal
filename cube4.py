import time
import math
import logging
import PyKEP as kep
import numpy as np
import matplotlib.pyplot as plt
import re

from launch import get_new_tle_and_satcat, findHistograms, findImportantAssets, loadData_tle, loadData_satcat, getExtrapolationLinearFunction
from collections import namedtuple

from itertools import combinations
import random
from breakup import setmass, breakup, plot_orbits
from operator import itemgetter
from mpi4py import MPI
from cubeAlone import cube, collision_prob, setradii, update
from removeObjects import removeRiskyObjects

def decayObjects(planets, currentYear):
    for p in planets:
        yearOfLaunch = int(p.name[0:4])
        # print yearOfLaunch
        decayYear = currentYear - 10
        if yearOfLaunch > 2016 and yearOfLaunch < decayYear:
            planets.remove(p)
            # print "removing"

# MPI identification
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc_name = MPI.Get_processor_name()
mpi_instance_id = "proc_" + '{0:03d}'.format(rank) + "_of_" + '{0:03d}'.format(comm.size) + "_on_" + proc_name

# Create logger.
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Create file handler to log debug messages.
fh = logging.FileHandler(mpi_instance_id + '.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(logging.Formatter('[%(asctime)s] %(message)s'))

# Create console handler to display info messages.
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(logging.Formatter('[{}] %(message)s'.format(mpi_instance_id)))

# Add new handlers to the logger.
logger.addHandler(fh)
logger.addHandler(ch)


# MPI identification
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()

file_inputParameters = "paramInput.txt"
fParam = open(file_inputParameters)
numOfExperiments = int(re.sub("[^0-9]", "",fParam.next()))

# length of experiment in years, standard setting 200 years
timeHorizon = int(re.sub("[^0-9]", "",fParam.next()))

# Russia
removedObjectsPerYear1 = int(re.sub("[^0-9]", "",fParam.next()))

# USA
removedObjectsPerYear2 = int(re.sub("[^0-9]", "",fParam.next()))

# China
removedObjectsPerYear3 = int(re.sub("[^0-9]", "",fParam.next()))

# EU
removedObjectsPerYear4 = int(re.sub("[^0-9]", "",fParam.next()))

fParam.close()

removedObjectsPerYear = [removedObjectsPerYear1, removedObjectsPerYear2, removedObjectsPerYear3, removedObjectsPerYear4]

# for cycle for multiple experiments - output file is marked with number of cycle
for qq in range(1,numOfExperiments+1):
    satcat = loadData_satcat()
    debris = loadData_tle()

    # list of all EU states and ESA and its bodies
    EU_list = ["EU","ASRA","BEL","CZCH","DEN","ESA","ESRO","EST","EUME","EUTE","FGER","FR","FRIT","GER","GREC","HUN","IT","LTU","LUXE","NETH","NOR","POL","POR","SPN","SWED","SWTZ","UK"]
    agentList = ['CIS', 'US', 'PRC', 'EU']
    agentListWithEU = ['CIS', 'US', 'PRC', EU_list] 
#  filtering only objects in LEO
    debris_LEO = []
    ep = kep.epoch(16.*365.25)
    for p in debris:
        try:
            oe = p.osculating_elements(ep)
            if oe[0] * (1-oe[1]) < 6378000.0 + 2000000:
                debris_LEO.append(p)
        except:
            pass
    debris = debris_LEO

    # osculating elements distributions
    oe_histograms = []
    oe_histograms = findHistograms(agentListWithEU)
    yearLaunchesDistributions = oe_histograms[10]

    beginningOfSequence = 2005
    expectedLifespan = 10
    # initialization of number of new launches per month
    numberOfNewLaunches = [3,3,3,3]
    numberOfLaunchesOthers = 15
    # rate of increase
    #very conservative increase of new launches
    numberOfNewLaunchesRateOfIncrease = 0.1
    removeHorizon = 2 #how often we remove objects - 1 -> every year, 10 -> every 10 years

    # list of recent active assets used for making "list of important assets" for each of the agents
    importantAssets = findImportantAssets(satcat, 2016, expectedLifespan)

    setradii(debris, satcat)
    setmass(debris)

    sump = 0
    decadeYearMultiplier = 1
    numberOfCatastrophicCollisions = 0
    numberOfNonCatastrophicCollisions = 0
    totalNumberOfNewDebris = 0
    totalDebrisEvol = []
    spaceDensities = []

    collisionRiskCouplesVirtual = []
    agentThreats = map(list, [[]]*len(agentList))
    agentRisksProb = map(list, [[]]*len(agentList))
    agentYearRiskProb = [0] * len(agentList)
    agentTotalImportantAssets = map(list, [[]]*len(agentList))
    agentCollisions = []
    for agentIndex in range(0,len(agentList)):
        agentCollisions.append([0] * timeHorizon) 
    # list of removed objects and who removed it
    listOfRemovedObjects = []
    listOfRemovedObjectsFinal = []

    objectsWithPossitiveProbOfCollisionAll = []
    objectsProbsOfCollisionsAll = []
    commonRisksProb = []
    commonRisksProbYear = 0

    stringOfRemovals = str(removeHorizon)+"year_rem"+str(removedObjectsPerYear1) + str(removedObjectsPerYear2) + str(removedObjectsPerYear3) + str(removedObjectsPerYear4)

    virtualRunEnabled = 1
    if removedObjectsPerYear1 == 0 and removedObjectsPerYear2 == 0 and removedObjectsPerYear3 == 0 and removedObjectsPerYear4 == 0:
        virtualRunEnabled = 0
    yearToPlot = 0
    ep = int(math.floor(16*365.25))
    month = 0
    virtualRun = 0 #initial setting, if 1 the run is virtual for computing probs of collisions
    # start of removing objects
    yearPom = 2016

    #TODO check if cloning is done properly
    virtualDebris = []
    virtualDebris = debris[:]
    dataSaving = yearPom + 1
    # time horizon in days from year 2000
    timeHorizonInDays = int(math.ceil(16*365.25 + timeHorizon*365.25))

    # main loop, one step is 5 days
    while ep < timeHorizonInDays or virtualRun == 1:
        ep += 5
        year = int(math.floor((ep-16*365.25)/365.25)+2016)

        if virtualRun == 1 and year%removeHorizon == 0 and year != yearPom and virtualRunEnabled:
            ep = ep - int(removeHorizon*365.25)
            virtualRun = 0
            yearPom = year
            year = int(math.floor((ep-16*365.25)/365.25)+2016)

#            agentThreats = [agent1Threats, agent2Threats, agent3Threats, agent4Threats]
            #TODO do better!!
            listOfRemovedObjects, collisionRiskCouplesVirtual = removeRiskyObjects(importantAssets, collisionRiskCouplesVirtual, agentThreats, agentList, removedObjectsPerYear, debris)
            listOfRemovedObjectsFinal.extend(listOfRemovedObjects)
            listOfRemovedObjects = []
#            print "Switching to real run ",year
            logger.info("Switching to real run " + str(year)+", total objects: "+str(len(debris)))

        if virtualRun == 0 and virtualRunEnabled and dataSaving == year:
            dataSaving += 1
#            update(importantAssetsTLE, ep)
            for agentIndex in range(0,len(agentList)):
                agentRisksProb[agentIndex].append(agentYearRiskProb[agentIndex])
            commonRisksProb.append(commonRisksProbYear)
            agentYearRiskProb = [0] * len(agentList)
            commonRisksProbYear = 0
            noAssetsAgent = [0] * len(agentList)

            # statistics - file output
            for planet in importantAssets.items():
                for agentIndex in range(0,len(agentList)):
                    if planet[1][5].strip() in agentListWithEU[agentIndex]:
                        noAssetsAgent[agentIndex] += 1

            for agentIndex in range(0,len(agentList)):
                agentTotalImportantAssets[agentIndex].append(noAssetsAgent[agentIndex])

        if virtualRun == 0 and year%removeHorizon == 0 and year == yearPom and year != (2016+timeHorizon) and virtualRunEnabled:
            virtualRun = 1
            yearPom = year
            virtualDebris = debris[:]
            logger.info("Switching to virtual run " + str(year))
#            print "Switching to virtual run ",year
        if not virtualRunEnabled and year == yearPom:
            yearPom += 1

            for agentIndex in range(0,len(agentList)):
                agentRisksProb[agentIndex].append(agentYearRiskProb[agentIndex])

            commonRisksProb.append(commonRisksProbYear)
            agentYearRiskProb = [0] * len(agentList)
            commonRisksProbYear = 0
            noAssetsAgent = [0] * len(agentList)

            # statistics - file output
            for planet in importantAssets.items():
                for agentIndex in range(0,len(agentList)):
                    if planet[1][5].strip() in agentListWithEU[agentIndex]:
                        noAssetsAgent[agentIndex] += 1

            for agentIndex in range(0,len(agentList)):
                agentTotalImportantAssets[agentIndex].append(noAssetsAgent[agentIndex])

        if yearToPlot != int(math.ceil((ep-16*365.25)/365.25)) and virtualRun == 0:
            yearToPlot = int(math.ceil((ep-16*365.25)/365.25))
            # statistics - output file
            totalDebrisEvol.append(len(debris))

        # getting month for relaunching sequences
        monthPom = month
        month = math.ceil((((ep-16*365.25)/365.25 +0.000001) - math.floor((ep-16*365.25)/365.25))*12)

# Launching new satellites

        if virtualRun == 0 and monthPom != month :
            numberOfNewLaunches2 = [0] * len(agentList)
            launchIterator = 1
            for agentIndex in range(0,len(agentList)):
                # TODO linear extrapolation needs to be finished
                # TODO focus on reasonable launching scenario !!
                linFunction_a, linFunction_b = getExtrapolationLinearFunction(yearLaunchesDistributions, agentIndex)
                if linFunction_a < 0:
                    linFunction_a = 0
                currentYearLaunchesNumber = np.multiply(linFunction_a,year) + linFunction_b
                # TODO linear growth not yet implemented

                # simple launch number growth
                numberOfNewLaunches2[agentIndex] = numberOfNewLaunches[agentIndex]  + (numberOfNewLaunchesRateOfIncrease*(year-2015))
                numberOfNewLaunches2[agentIndex]  = int(round(numberOfNewLaunches2[agentIndex] ))
#                print currentYearLaunchesNumber, numberOfNewLaunches2[agentIndex]
                addingNewPlanet = 0
                for i in range(0,numberOfNewLaunches2[agentIndex] ):
                    newTLEobject, newSATCATobject = get_new_tle_and_satcat(agentIndex+1, oe_histograms, year, month, launchIterator)
                    debris.append(newTLEobject)
                    satcat[newSATCATobject[0]] = newSATCATobject[1]
                    addingNewPlanet = 1
                    launchIterator  += 1

                # launching of other agents
            for i in range(0,numberOfLaunchesOthers):
                newTLEobject, newSATCATobject = get_new_tle_and_satcat(len(agentList)+1, oe_histograms, year, month, launchIterator)
                debris.append(newTLEobject)
                satcat[newSATCATobject[0]] = newSATCATobject[1]
                addingNewPlanet = 1
                launchIterator  += 1
#
            if addingNewPlanet == 1:
                setmass(debris)
        #getting a list of active recent planets to form list of important assets for each agent
        importantAssets = findImportantAssets(satcat, year, expectedLifespan)

        if virtualRun == 1:
            # print "virtual", len(virtualDebris)
            update(virtualDebris, ep)
            volumes = cube(virtualDebris)
        else:
            # print "real", len(debris)
#            decayObjects(debris, year)
            update(debris, ep)
            volumes = cube(debris)

        maxp = 0
        for volume in volumes:
            for p1, p2 in combinations(volume, 2):
                if p1.name == p2.name and p1.mass ==  p2.mass:
                    continue
                # print p1.name,p2.name
                p1.name = p1.name.strip()
                p2.name = p2.name.strip()
                if tuple(p1._v) == tuple(p2._v):
                    pass
                else:
                    # probability of collision
                    Pij = collision_prob(p1, p2) * 2
                    # probability of collision over 5 days
                    P = Pij * 5. * kep.DAY2SEC
                    maxp = max(maxp, P)
                    sump += P
                    if P > 0 and virtualRun == 1:

                        if importantAssets.has_key(p1.name.strip()) and importantAssets.has_key(p2.name.strip()):
                            plnt1 = importantAssets[p1.name.strip()]
                            plnt2 = importantAssets[p2.name.strip()]
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(), plnt1[5].strip(), plnt2[5].strip(), P])
                        elif importantAssets.has_key(p1.name.strip()) and not importantAssets.has_key(p2.name.strip()):
                            plnt1 = importantAssets[p1.name.strip()]
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(), plnt1[5].strip(), "ROW", P])
                        elif not importantAssets.has_key(p1.name.strip()) and importantAssets.has_key(p2.name.strip()):
                            plnt2 = importantAssets[p2.name.strip()]
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(), "ROW", plnt2[5].strip(), P])
                        else:
                            collisionRiskCouplesVirtual.append([p1.name.strip(), p2.name.strip(),"ROW", "ROW", P])

                        for agentIndex in range(0,len(agentList)):
                            if importantAssets.has_key(p1.name.strip()):
                                plnt1 = importantAssets[p1.name.strip()]
                                if plnt1[5].strip() in agentListWithEU[agentIndex]:
                                    for nameOfObject in agentThreats[agentIndex]:
                                        if nameOfObject[0] == p2.name.strip():
                                            nameOfObject[1] += P
                                            break
                                    else:
                                        agentThreats[agentIndex].append([p2.name.strip(), P])

                        for agentIndex in range(0,len(agentList)):
                            if importantAssets.has_key(p2.name.strip()):
                                plnt2 = importantAssets[p2.name.strip()]
                                if plnt2[5].strip() in agentListWithEU[agentIndex]:
                                    for nameOfObject in agentThreats[agentIndex]:
                                        if nameOfObject[0] == p1.name.strip():
                                            nameOfObject[1] += P
                                            break
                                    else:
                                        agentThreats[agentIndex].append([p1.name.strip(), P])

                    if P > 0 and virtualRun == 0:
                        commonRisksProbYear += P
                        if importantAssets.has_key(p1.name.strip()):
                            planetA = importantAssets[p1.name.strip()]
                            for agentIndex in range(0,len(agentList)):
                                if planetA[5].strip() in agentListWithEU[agentIndex]:
                                # print "Russian object is threatened by ", p2.name
                                    agentYearRiskProb[agentIndex] += P

                        if importantAssets.has_key(p2.name.strip()):
                            planetB = importantAssets[p2.name.strip()] 
                            for agentIndex in range(0,len(agentList)):
                                if planetB[5].strip() in agentListWithEU[agentIndex]:
                                # print "Russian object is threatened by ", p2.name
                                    agentYearRiskProb[agentIndex] += P

                    # Monte-carlo simulation - there is a collision if random number is lower than prob of collision
                    if random.random() < P and virtualRun == 0:

                        for planet in importantAssets.items():
                            # if agent's assset is in risk of collision we save the object which threaten it
                            indexInArray = year-2016
                            for agentIndex in range(0,len(agentList)):
                                if (planet[0].strip() == p1.name.strip() or planet[0].strip() == p2.name.strip()) and planet[1][5].strip() in agentListWithEU[agentIndex]:
                                    agentCollisions[agentIndex][indexInArray] += 1


                        logger.info("!" * 25 + " BOOOOOM " + "!" * 25)
                        logger.info('planet {} with mass {}'.format(p1.name, p1.mass))
                        logger.info('planet {} with mass {}'.format(p2.name, p2.mass))
                        dv = np.linalg.norm(p1._v - p2._v)
                        logger.info('collision velocity ' +  str(dv))

                        if p2.mass < p1.mass:
                            catastrophRatio = (p2.mass*dv**2)/(2*p1.mass*1000)
                        else:
                            catastrophRatio = (p1.mass*dv**2)/(2*p2.mass*1000)
                        logger.info('catastrophic ratio: ' + str(catastrophRatio))
#                        print 'catastrophic ratio:', catastrophRatio
                        if catastrophRatio<40:
#                            print 'NON CATASTROPHIC COLLISION - NO DEBRIS'
                            logger.info('NON CATASTROPHIC COLLISION - NO DEBRIS')
                            numberOfNonCatastrophicCollisions += 1
                        else:
                            debris1, debris2 = breakup(ep, p1, p2)
                            logger.info(str(len(debris1) + len(debris2)) + ' new pieces from collision')
#                            print (len(debris1) + len(debris2)), 'new pieces from collision'
                            totalNumberOfNewDebris += (len(debris1) + len(debris2))
                            debris.extend(debris1)
                            debris.extend(debris2)
                            numberOfCatastrophicCollisions += 1
                            setmass(debris)
        logger.debug('%.2f %d %d %d %10.8f' % (float(ep)/365.25 - 16, len(debris), len(volumes), max(map(len, volumes)) if len(volumes) else 0, maxp))
    logger.info('There were ' + str(numberOfCatastrophicCollisions) + ' catastrophic collisions')
    logger.info('There were ' + str(numberOfNonCatastrophicCollisions) + ' NON-catastrophic collisions')
    logger.info('From collisions there were ' + str(totalNumberOfNewDebris) + ' of new debris')
    years = range(0,yearToPlot+1,1)

    totalDebrisEvol.append(len(debris))

    fileName = "exper"+stringOfRemovals+"thread"+str(rank)+"sim"+str(qq)+".txt"
    f = open(fileName,'w')
    f.write("%s\n" % years)
    f.write("%s\n" % totalDebrisEvol)
    f.write("%s\n" % numberOfCatastrophicCollisions)
    f.write("%s\n" % numberOfNonCatastrophicCollisions)
    f.write("%s\n" % totalNumberOfNewDebris)
    for agentIndex in range(0,len(agentList)):
        f.write("%s\n" % agentCollisions[agentIndex])
    for agentIndex in range(0,len(agentList)):
        f.write("%s\n" % agentRisksProb[agentIndex])
    f.write("%s\n" % commonRisksProb)
    for agentIndex in range(0,len(agentList)):
        f.write("%s\n" % agentTotalImportantAssets[agentIndex])
    for remObject in listOfRemovedObjectsFinal:
        f.write("%s\n" % remObject)
    f.close()

