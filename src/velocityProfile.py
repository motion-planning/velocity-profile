from pylab import figure, show, rand
from matplotlib.patches import Ellipse
import matplotlib as mpl
import numpy as np
import math
from numpy import *


def perp(a):
    b = empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return


def seg_intersect(a1, a2, b1, b2):
    da = a2 - a1
    db = b2 - b1
    dp = a1 - b1
    dap = perp(da)
    denom = dot(dap, db)
    num = dot(dap, dp)
    return (num / denom) * db + b1


def ellipse_general_to_standard(A, B, C, D, E, F):
    """
    Convert an ellipse in general form:
      Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
    To standard form:
      ((x - x0) cos t - (y - y0) sin t)^2/a^2
      + ((x - x0) sin t + (y - y0) cos t)^2/b^2 = 1
    The ellipse has center (x0, y0), major and minor semi-axes a and b,
    and the angle to the semi-major axis is t.

    Parameters: A, B, C, D, E, F
    Returns: x0, y0, a, b, t
    """
    A, B, C, D, E, F = map(float, [A, B, C, D, E, F])
    # Matrix representation of conic section
    AQ = np.array([[A, B / 2, D / 2], [B / 2, C, E / 2], [D / 2, E / 2, F]])
    A33 = AQ[0:2, 0:2]
    # Formula for center
    x0, y0 = np.linalg.inv(A33).dot([-D / 2, -E / 2])
    # Each eigenvector of A33 lies along one of the axes
    evals, evecs = np.linalg.eigh(A33)
    # Semi-axes from reduced canonical equation
    a, b = np.sqrt(-np.linalg.det(AQ) / (np.linalg.det(A33) * evals))
    # Return major axis as "a" and angle to major axis between -pi/2 and pi/2
    t = np.arctan2(evecs[1, 0], evecs[0, 0])  # angle of axis "a"
    if b > a:
        a, b = b, a
        t += np.pi / 2
    if t < -np.pi / 2:
        t += np.pi
    elif t > np.pi / 2:
        t -= np.pi
    return (x0, y0, a, b, t)


def ellipse_parameters(a, b, x0, y0, xc, yc, vx, vy, r):
    dx = x0 - xc
    dy = y0 - yc
    A = a * a + b * b
    B = -2 * (a * vx + b * vy)
    C = vx * vx + vy * vy
    D = 2 * (a * dx + b * dy)
    E = -2 * (vx * dx + vy * dy)
    F = dx * dx + dy * dy - r * r
    return A, B, C, D, E, F


def calculateLineLength(line):
    base = math.sqrt((line[0] - line[2]) ** 2 + (line[1] - line[3]) ** 2)
    return base


def isPointonLine(ponitx, pointy, line):
    # if the point is on line then the slope of the point w.r.t. line
    # endpoints and slope of the line will be the same
    slope = (line[3] - line[1]) / (0.0001 + line[2] - line[0])
    slope1 = (pointy - line[1]) / (0.0001 + ponitx - line[0])
    if format(slope, '.2f') == format(slope1, '.2f'):
        return True
    else:
        return False


def rearrangeObstacleLines(obstacle, basetimeObstacle, velocityObstacle):
    i = 0
    length = 0
    tempobstacle = []
    startLine = 0
    # continue on the all obstacle lines untl to find where is the obstacle is
    # now for the given basetime
    for obstacleLine in obstacle:
        length += calculateLineLength(obstacleLine)
        if length / velocityObstacle > basetimeObstacle:
            startLine = i  # index of the line where the current obstacle is
            length -= calculateLineLength(obstacleLine)
            break
        i += 1
    # find the exact fractional line segment length of current line where the
    # obstacle is
    newLength = velocityObstacle * basetimeObstacle - length
    slope = (obstacleLine[3] - obstacleLine[1]) / \
        (0.0001 + obstacleLine[2] - obstacleLine[0])
    # coordinate of the obstacle on the current line it will be the starting
    # point of the obstacle now.
    xnew = math.sqrt(newLength ** 2 / (1 + slope ** 2)) + obstacleLine[0]
    ynew = slope * (xnew - obstacleLine[0]) + obstacleLine[1]
    # obstacle[startLine][0]=xnew
    # obstacle[startLine][1]=ynew
    # print obstacle
    i = 0
    # update line array.current line where the obstacle is will be the first
    # line.then other lines will be appended
    for obstacleLine in obstacle:
        # print i,(startLine+i)%len(obstacle)
        tempobstacle.append(obstacle[(startLine + i) % len(obstacle)])
        i += 1
    return tempobstacle, xnew, ynew


def findEllipse(lineObstacle, lineRobot, velocityObstacle, basetimeObstacle, basetimeRobot, obstacleStartX, obstacleStartY):
    # initialize for initial point checking of obstacle. after each robot line
    # solved,an obstacle has changed its position over time.
    tempx = 0
    tempy = 0
    # check if the initial position of obstacle is on this line. if so then
    # change the lines initial coordinates(x0,y0).
    if isPointonLine(obstacleStartX, obstacleStartY, lineObstacle):
        tempx = lineObstacle[0]
        tempy = lineObstacle[1]
        lineObstacle[0] = obstacleStartX
        lineObstacle[1] = obstacleStartY
    # print lineObstacle
    # calculate parameters for line equation. We need a,b,x0,y0.
    a = -(lineRobot[0] - lineRobot[2])
    b = -(lineRobot[1] - lineRobot[3])
    x0 = lineRobot[0]
    y0 = lineRobot[1]
    xc = lineObstacle[0]  # parameters for obstacle initial point.
    yc = lineObstacle[1]
    # calculate slope for the obstacle line
    slope = (lineObstacle[3] - lineObstacle[1]) / \
        (0.0001 + lineObstacle[2] - lineObstacle[0])
    # calculate Vx from the V for obstacle
    vx = math.sqrt(velocityObstacle ** 2 / (1 + slope ** 2))
    # calculate Vy from the V for obstacle
    vy = math.sqrt(velocityObstacle ** 2 - vx ** 2)
    # set the sign for Vy. If going downward then the sign is negative
    if lineObstacle[3] - lineObstacle[1] < 0:
        vy = -1 * vy
    # set the sign for Vx. If going leftward then the sign is negative
    if lineObstacle[2] - lineObstacle[0] < 0:
        vx = -1 * vx
    r = 5  # set the radius for the circular obstacle
    # print vx,vy
    # Calculate the conic equation Ax^2+Bxy+Cy^2+Dx+Ey+F=0
    A, B, C, D, E, F = ellipse_parameters(a, b, x0, y0, xc, yc, vx, vy, r)
    centerx, centery, h, w, t = ellipse_general_to_standard(
        A, B, C, D, E, F)  # Calculate the ellipse parameters

    # check if this is a valid rectangle.this means if the ellipse falls in
    # the line or outside the line length.the line starts from 0 and ends at 1
    if centerx < 0 or centerx > 1:
        centerx, centery, h, w, t = [0, 0, 0, 0, 0]
    # Re-Calculate the ellipse center for space coordinate.till now
    # space[0,1].we must convert it to real length in space by multiplying by
    # root(a^2+b^2).
    centerx = centerx * math.sqrt(a * a + b * b) + basetimeRobot
    # then add the current advancement of the robot by already traversed previous lines in the path.
    # Re-Calculate the ellipse center for time coordinate.Remember the time
    # coordinate is for obstacle time. The time always advances
    # .basetimeobstacle is the time of
    # #starting the obstacle for current obstacle line. basetimeobstacle
    # depends on how much the robot took for previous line.
    centery = centery + basetimeObstacle
    # then add the current advancement of the robot by already traversed previous lines in the path.
    # update the basetimeobstacle  after finishing current obstacle line
    basetimeObstacle += math.sqrt((lineObstacle[0] - lineObstacle[2]) ** 2 + (
        lineObstacle[1] - lineObstacle[3]) ** 2) / velocityObstacle
    # restore the original obstacle line which we might changed at the
    # starting of this function
    if tempx != 0:
        lineObstacle[0] = tempx
        lineObstacle[1] = tempy
    print centerx, centery, h, w, t
    return centerx, centery, h, w * math.sqrt(a * a + b * b), t, basetimeObstacle


def lineIntersection(line1, line2):
    value1 = (line1[2] - line1[0]) * (line2[1] - line1[3]) - \
        (line1[3] - line1[1]) * (line2[0] - line1[2])
    value2 = (line1[2] - line1[0]) * (line2[3] - line1[3]) - \
        (line1[3] - line1[1]) * (line2[2] - line1[2])
    value3 = (line2[2] - line2[0]) * (line1[1] - line2[3]) - \
        (line2[3] - line2[1]) * (line1[0] - line2[2])
    value4 = (line2[2] - line2[0]) * (line1[3] - line2[3]) - \
        (line2[3] - line2[1]) * (line1[2] - line2[2])
    if ((value1 > 0 and value2 < 0) or (value1 < 0 and value2 > 0)) and ((value3 > 0 and value4 < 0) or (value3 < 0 and value4 > 0)):
        print "intersection found"
        p1 = array([line1[0], line1[1]])
        p2 = array([line1[2], line1[3]])
        p3 = array([line2[2], 0])
        p4 = array([line2[2], line2[3]])
        print line1, line2
        return seg_intersect(p1, p2, p3, p4)
    return -1, -1


def velocityProfile(robotLineSpaceTimeX, robotLineSpaceTimeY, robotLine, rectList, velocityRobot, velocityLineList):
    # the endpoint of the line will be startpoint+the length of line
    lineEndX = calculateLineLength(robotLine) + robotLineSpaceTimeX
    # the endpoint of the line in time coordinate is different.it will be
    # starttime+the time reguired to travel along the line
    lineEndY = calculateLineLength(
        robotLine) / velocityRobot + robotLineSpaceTimeY
    # represents the line in space time coordinate system
    lineinSpaceTime = [
        robotLineSpaceTimeX, robotLineSpaceTimeY, lineEndX, lineEndY]
    # boolean variable to check if any intersection with any obstacle is found
    # or not. initialize it with false as no intersection at the beginning
    intersectionFound = False
    print 'linespacetime', lineinSpaceTime
    # for each rectangular obstacle we must check for possible intersection
    for rect in rectList:
        # represents one of the four lines of the rectangular obstacle.this one
        # is the vertical line on left side
        lineinSpaceTime1 = [rect[0], rect[1], rect[2], rect[3]]
        # this function returns the line intersection point. if no intersection
        # found it returns -1,-1
        x, y = lineIntersection(lineinSpaceTime, lineinSpaceTime1)
        if x != -1:
            print 'x,y', x, y
            # set the boolean variable true that an intersection has been found
            intersectionFound = True
            # append a line into the velocity list.the startpoint of this line
            # is the startpoint of the robot on this line.the endpoint of this
            # line is the intersection point
            velocityLineList.append(
                [robotLineSpaceTimeX, robotLineSpaceTimeY, x, y])
            # the startpoint is changed to the top left corner of the
            # rectangular obstacle
            robotLineSpaceTimeX = x
            robotLineSpaceTimeY = rect[1]
            lineEndY = ((lineEndY - y) / (lineEndX - x)) * \
                (lineEndX - rect[0]) + rect[1]
            lineinSpaceTime = [
                robotLineSpaceTimeX, robotLineSpaceTimeY, lineEndX, lineEndY]
            # print 'new linespacetime',lineinSpaceTime
            velocityLineList.append([x, y, x, rect[1]])
        # represents one of the four lines of the rectangular obstacle.this one
        # is the horizontal line on down side
        lineinSpaceTime1 = [rect[4], rect[5], rect[2], rect[3]]
        # this function returns the line intersection point. if no intersection
        # found it returns -1,-1
        x, y = lineIntersection(lineinSpaceTime, lineinSpaceTime1)
        if x != -1:
            print 'x,y', x, y
            intersectionFound = True
            # append a line into the velocity list.the startpoint of this line
            # is the startpoint of the robot on this line.the endpoint of this
            # line is the intersection point
            velocityLineList.append(
                [robotLineSpaceTimeX, robotLineSpaceTimeY, x, y])
            robotLineSpaceTimeX = x
            robotLineSpaceTimeY = rect[1]
            lineEndY = ((lineEndY - y) / (lineEndX - x)) * \
                (lineEndX - rect[0]) + rect[1]
            lineinSpaceTime = [
                robotLineSpaceTimeX, robotLineSpaceTimeY, lineEndX, lineEndY]
            velocityLineList.append([x, y, x, rect[1]])
    # if no intersection found then add the original line in the space time
    # velocity profile
    if intersectionFound == False:
        robotLineSpaceTimeX = lineEndX
        robotLineSpaceTimeY = lineEndY
        velocityLineList.append(lineinSpaceTime)
    if intersectionFound == True:
        robotLineSpaceTimeX = lineEndX
        robotLineSpaceTimeY = lineEndY
        # print 'new linespacetime',lineEndX,lineEndY
        velocityLineList.append(lineinSpaceTime)
    return robotLineSpaceTimeX, robotLineSpaceTimeY, velocityLineList


def main():
    obstacle = [[9, 45, 40, 42], [40, 42, 40, 23],
                [40, 23, 8, 7], [8, 7, 9, 45]]
    velocityObstacle = 3
    basetimeObstacle = 0
    obstacleStartX = 0
    obstacleStartY = 0
    # robot=[[1,5,15,20],[15,20,19,48],[19,48,33,37],[33,37,36,5]]
    robot = [
        [0.00, 0.00, 10.00, 10.00],
        [10.00, 10.00, 20.00, 20.00],
        [20.00, 20.00, 30.00, 30.00],
        [30.00, 30.00, 35.03, 35.03],
        [35.03, 35.03, 33.53, 40.00],
        [33.53, 40.00, 31.51, 50.00],
        [31.51, 50.00, 30.50, 60.00],
        [30.50, 60.00, 30.50, 69.50],
        [30.50, 69.50, 40.00, 77.42],
        [40.00, 77.42, 50.00, 83.90],
        [50.00, 83.90, 60.00, 88.53],
        [60.00, 88.53, 61.15, 88.85],
        [61.15, 88.85, 70.00, 92.46],
        [70.00, 92.46, 80.00, 94.50],
        [80.00, 94.50, 90.00, 94.50],
        [90.00, 94.50, 93.74, 93.74],
        [93.74, 93.74, 96.49, 90.00],
        [96.49, 90.00, 100.62, 80.62],
        [100.62, 80.62, 101.05, 80.00],
        [101.05, 80.00, 104.50, 70.00],
        [104.50, 70.00, 104.50, 60.00],
        [104.50, 60.00, 101.05, 50.00],
        [101.05, 50.00, 100.62, 49.38],
        [100.62, 49.38, 96.49, 40.00],
        [96.49, 40.00, 89.13, 30.00],
        [89.13, 30.00, 79.43, 20.57],
        [79.43, 20.57, 70.00, 21.77],
        [70.00, 21.77, 60.00, 24.30],
        [60.00, 24.30, 89.13, 30.00]]
    baseSpaceRobot = 0
    velocityRobot = 2
    fig = figure()
    ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")
    ellipseList = []
    robotLineSpaceTimeX, robotLineSpaceTimeY = 0, 0
    velocityLineList = []
    for robotLine in robot:
        print "####################################new run###################################"
        print obstacle, obstacleStartX, obstacleStartY, baseSpaceRobot, basetimeObstacle
        for obstacleLine in obstacle:
            centerx, centery, h, w, t, basetimeObstacle = findEllipse(
                obstacleLine, robotLine, velocityObstacle, basetimeObstacle, baseSpaceRobot, obstacleStartX, obstacleStartY)
            if h != 0:
                ellipseList.append([centerx - w / 2, centery + h / 2, centerx - w / 2, centery -
                                    h / 2, centerx + w / 2, centery - h / 2, centerx + w / 2, centery + h / 2])
            ellipse = mpl.patches.Rectangle(
                xy=(centerx - w / 2, centery - h / 2), width=w, height=h)
            fig.gca().add_artist(ellipse)
        robotLineSpaceTimeX, robotLineSpaceTimeY, velocityLineList = velocityProfile(
            robotLineSpaceTimeX, robotLineSpaceTimeY, robotLine, ellipseList, velocityRobot, velocityLineList)
        baseSpaceRobot += calculateLineLength(robotLine)
        basetimeObstacle = robotLineSpaceTimeY
        obstacle, obstacleStartX, obstacleStartY = rearrangeObstacleLines(
            obstacle, basetimeObstacle, velocityObstacle)
    #aList = [123, 'xyz', 'zara', 'abc', 'xyz'];

    velocityLineList.sort()
    tempvelocityLineList = []
    prevLineX = velocityLineList[0][0]
    prevLineY = velocityLineList[0][1]
    i = 0
    # velocityLineList=[]
    print tempvelocityLineList
    for line in velocityLineList:
        if line[0] < prevLineX and line[1] < prevLineY and line[0] == line[2]:
            tempvelocityLineList.append(line)
            tempvelocityLineList[i - 1][2] = line[0]
            tempvelocityLineList[i - 1][3] = line[1]
            prevLineX = line[2]
            prevLineY = line[3]
        if line[0] == prevLineX and line[1] == prevLineY:
            tempvelocityLineList.append(line)
            print line, prevLineX, prevLineY
            prevLineX = line[2]
            prevLineY = line[3]
        print prevLineX, prevLineY
        i = i + 1

    print tempvelocityLineList

    for line in tempvelocityLineList:
        fig.gca().plot(
            (line[0], line[2]), (line[1], line[3]), 'go-', label='line 1', linewidth=2)

    ax.set_xlim(0, 120)
    ax.set_ylim(0, 120)
    show()


if __name__ == '__main__':
    main()