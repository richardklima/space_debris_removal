import copy
import logging
from operator import itemgetter
import math

def updateCollisionTuples(tupleOfCollisions, removedObjects):
    newListOfCommonRiskPlanets = []
    for pln in tupleOfCollisions:
        if pln[0] in removedObjects:
            tupleOfCollisions.remove(pln)
        elif pln[1] in removedObjects:
            tupleOfCollisions.remove(pln)

    for couple in tupleOfCollisions:
        for planet in newListOfCommonRiskPlanets:
            if couple[0] == planet[0]:
                planet[1] += couple[4]
                break
        else:
            newListOfCommonRiskPlanets.append([couple[0],couple[4]])

        for planet in newListOfCommonRiskPlanets:
            if couple[1] == planet[0]:
                planet[1] += couple[4]
                break
        else:
            newListOfCommonRiskPlanets.append([couple[1],couple[4]])

    return newListOfCommonRiskPlanets, tupleOfCollisions

def removeRiskyObjects(importantAssets, collisionRiskCouplesVirtual, agentThreats, agentList, removedObjectsPerYear, debris):
    listOfRemovedObjects = []
    logger = logging.getLogger(__name__)
    for threats in agentThreats:
   # we do not remove any important assets
        for asset in threats:
            if importantAssets.has_key(asset[0]):
                threats.remove(asset)

    # need to use deep copy because of cloning list of lists of objects
    pomThreat = []
    pomThreats = copy.deepcopy(agentThreats)

    # remove objects until n objects are removed
    listOfRemoved = []
    iterator = [0] * len(agentList)
#   iterator1 = 0
    for agentIndex in range(0, len(agentList)):
        while iterator[agentIndex] < removedObjectsPerYear[agentIndex] and len(agentThreats[agentIndex]) != 0:
            highestThreat = max(agentThreats[agentIndex], key = itemgetter(1))
            agentThreats[agentIndex].remove(highestThreat)
            for d in debris:
                if highestThreat[0] == d.name:
#                    print agentList[agentIndex]+" removing", d.name, highestThreat[1]
                    logger.info(agentList[agentIndex]+" removing {} {}".format(d.name, highestThreat[1]))
                    debris.remove(d)
                    listOfRemovedObjects.append([agentList[agentIndex],"o",d.name,0,0,0,0])
                    listOfRemoved.append(d.name)
                    iterator[agentIndex] += 1
                    break
    
    listOfCommonRiskPlanets = []

    listOfCommonRiskPlanets, collisionRiskCouplesVirtual = updateCollisionTuples(collisionRiskCouplesVirtual, listOfRemoved)

    # REMOVING COMMON OBJECTS

    for agentIndex in range(0, len(agentList)):
        while iterator[agentIndex] < removedObjectsPerYear[agentIndex] and len(listOfCommonRiskPlanets) != 0:
            if iterator[agentIndex] < removedObjectsPerYear[agentIndex] and len(listOfCommonRiskPlanets) != 0:

                highestThreat = max(listOfCommonRiskPlanets, key = itemgetter(1))
                listOfCommonRiskPlanets.remove(highestThreat)
                # we do not remove important assets
                if importantAssets.has_key(highestThreat[0]):
                    continue
                for d in debris:
                    if highestThreat[0] == d.name:
#                        print agentList[agentIndex]," removing common", d.name, highestThreat[1]
                        logger.info(agentList[agentIndex]+" removing common {} {}".format(d.name, highestThreat[1]))
                        debris.remove(d)
                        listOfRemovedObjects.append([agentList[agentIndex],"c",d.name,0,0,0,0])
                        listOfCommonRiskPlanets, collisionRiskCouplesVirtual = updateCollisionTuples(collisionRiskCouplesVirtual, d.name)
                        iterator[agentIndex] += 1
                        break

    for remObj in listOfRemovedObjects:
        for agentIndex in range(0, len(agentList)):
            for threat in pomThreats[agentIndex]:
                if remObj[2] == threat[0]:
                    remObj[agentIndex+3] += threat[1]
#                    print remObj[2]

    return listOfRemovedObjects, collisionRiskCouplesVirtual 


