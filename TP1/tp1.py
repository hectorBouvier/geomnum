#------------------------------------------------------
#
#  TP1 : Bezier curves, De Casteljau's algorithm
#  http://tiborstanko.sk/teaching/geo-num-2017/tp1.html
#  [03-Feb-2017]
#
#------------------------------------------------------
#
#  This file is a part of the course:
#    Geometrie numerique (spring 2017)
#    https://github.com/GeoNumTP/GeoNum2017
#    M1 Informatique
#    UFR IM2AG
#
#  Course lecturer:
#    Georges-Pierre.Bonneau at inria.fr
#
#  Practical part:
#    Tibor.Stanko at inria.fr
#
#------------------------------------------------------

import sys, os
import matplotlib.pyplot as plt
import numpy as np

TP = os.path.dirname(os.path.realpath(__file__)) + "/"
DATADIR = filename = TP+"data/"


#-------------------------------------------------
# READPOLYGON()
# Read Bezier control points from a file.
#
def ReadPolygon( filename ) :
    datafile = open(filename,'r');
    l = datafile.readline()
    degree = np.fromstring(l,sep=' ',dtype=int)
    BezierPts = np.fromfile(datafile,count=2*(degree+1),sep=' ',dtype=float)
    BezierPts = BezierPts.reshape(-1,2)
    return BezierPts


#-------------------------------------------------
# DECASTELJAU( ... )
# Perform the De Casteljau algorithm.
#
# Input
#    BezierPts :  (degree+1) x 2 matrix of Bezier control points
#    k         :  upper index of the computed point (depth of the algorithm)
#    i         :  lower index of the computed point
#    t         :  curve parameter in [0.0,1.0]
#
# Output
#    point b_i^k from the De Casteljau algorithm.
#
def DeCasteljau( BezierPts, k, i, t ) :
    if k==0:
        return BezierPts[i]
    else: 
        return t*DeCasteljau(BezierPts,k-1,i+1,t)+(1-t)*DeCasteljau(BezierPts,k-1,i,t)


#-------------------------------------------------
# BEZIERCURVE( ... )
# Compute points on the Bezier curve.
#
# Input
#    BezierPts :  (degree+1) x 2 matrix of Bezier control points
#    N         :  number of curve samples
#    
# Output
#    CurvePts  :  N x 2 matrix of curvepoints
#
def BezierCurve( BezierPts, N ) :
    
    # degree of the curve (one less than the number of control points)
    degree = BezierPts.shape[0]-1
    
    # initialize curvepoints as zeros
    CurvePts = np.zeros([N,2])
    T=np.linspace(0.0,1.0,num=N)
    for i in range(N):    
        CurvePts[i,:]=DeCasteljau(BezierPts,degree,0,T[i])
    return CurvePts

def AffichageDeCasteljau( BezierPts, k) :
    t=0.5
    if k!=0:
        B=np.zeros([k-1,2])
        for j in range(k-1):
            B[j]=(1-t)*BezierPts[j]+t*BezierPts[j+1]
        plt.plot(B[:,0],B[:,1])
        AffichageDeCasteljau(B,k-1)
#-------------------------------------------------
if __name__ == "__main__":
    
    # arg 1 : data name 
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
    else :
        dataname = "simple" # simple, infinity, spiral

    # arg 2 : sampling density
    if len(sys.argv) > 2 :
        density = int(sys.argv[2])
    else :
        density = 10

    # filename
    filename = DATADIR + dataname + ".bcv"
    
    # check if valid datafile
    if not os.path.isfile(filename) :
        print "error:  invalid dataname '" + dataname + "'"
        print "usage:  python tp1.py  [simple,infinity,spiral]  [sampling_density]"
        
    else :    
        # read control points
        BezierPts = ReadPolygon(filename)

        # compute curve points
        CurvePts = BezierCurve(BezierPts,density)

        # plot
        fig = plt.gcf()
        fig.canvas.set_window_title('TP1 Bezier curves')
        plt.title(dataname+', '+str(density)+" pts")
        
        # set axes with equal proportions
        plt.axis('equal')
        
        # plot the control polygon
        plt.plot( BezierPts[:,0], BezierPts[:,1], 'bo-' )
        
        # plot the curve
        plt.plot( CurvePts[:,0], CurvePts[:,1], 'r-' )
        #output the intermediate polygons of Bezier. Comment for a more clear output
        AffichageDeCasteljau( BezierPts,BezierPts.shape[0]-1) 
        #########
        ## TODO : Uncomment if you want to save the render as png image in the data/ dir
        #########
        plt.savefig( DATADIR + dataname + ".png" )
        
        #########
        ## TODO : Compute intermediate polygons b_i^k for k=1,...,degree-1 and i=0,...,degree-k
        #########
        
        #########
        ## TODO : Add plt.plot commands to plot the intermediate polygons
        #########
        
        plt.show()
