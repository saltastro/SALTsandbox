import math
import numpy as np

def imalign(xa, ya, xb, yb, maga=None, magb=None, tolerance=0, niter=5):
    """Register  a list of images by computing relative object
    shifts.  Given an input set of coordinates in xa and ya it will
    compute the shifts between the coordinates
    """
    pass

def immatch(xa, ya, xb, yb, maga=None, magb=None):
    """Given a set of coordinates, prodce an output list
       of matched coordinates 
    """ 
    pass

def imshift(xa, ya, xb, yb, fa=None, fb=None):
    """Assume a linear x-y shift between images,
    correlate the two lists of images to find the 
    best x-y shift

    Parameters
    -----------
    xa : `~numpy.ndarray`
        x position for sources in image a
    ya : `~numpy.ndarray`
        y position for sources in image a
    xb : `~numpy.ndarray`
        x position for sources in image b
    yb : `~numpy.ndarray`
        y position for sources in image b
    fa : `~numpy.ndarray`
        Flux  for sources in image a.  Can be set to None.
    fb : `~numpy.ndarray`
        Flux  for sources in image b.  Can be set to None.

    """
    if fa is None: fa=xa*0.0+1.0
    if fb is None: fb=xa*0.0+1.0

def measureangle(xa, ya, xb, yb, xc, yc):
    """Measure the angle from points b to c assuming point a is at the origin"""
    x1=xb-xa
    x2=xc-xa
    y1=yb-ya
    y2=yc-ya
    if y1==0:
       t1=0.5*math.pi
    else:
       t1=math.atan(x1/y1)
    if y2==0:
       t2=0.5*math.pi
    else:
       t2=math.atan(x2/y2)
    t=t1+t2
    if t>0.5*math.pi:  t=math.pi-t
    return math.cos(t)


def measuretriangle(xa, ya, xb, yb, xc, yc):
    """Measue the length of the perimeter, ratio of the
       longest to shortest side, and cosine of the angle
       between the longest and shortest side.

       Returns the list of these and the direction of 
       the arangement of vertices
    """
    s1=((xa-xb)**2+(ya-yb)**2)**0.5
    s2=((xa-xc)**2+(ya-yc)**2)**0.5
    s3=((xc-xb)**2+(yc-yb)**2)**0.5
    perimeter=s1+s2+s3

    long_side=max(s1,s2,s3)
    short_side=min(s1,s2,s3)
    ratio=short_side/long_side

    #determine the angle between the longest and shortest side
    if s1==long_side and s2==short_side: 
       cosang=measureangle(xa, ya, xb, yb, xc, yc)
       vert=0
       orient=1
    if s2==long_side and s1==short_side:
       cosang=measureangle(xa, ya, xc, yc, xb, yb)
       vert=0
       orient=-1

    if s1==long_side and s3==short_side: 
       cosang=measureangle(xb, yb, xa, ya, xc, yc)
       vert=1
       orient=-1
    if s3==long_side and s1==short_side:
       cosang=measureangle(xb, yb, xc, yc, xa, ya)
       vert=1
       orient=1

    if s3==long_side and s2==short_side: 
       cosang=measureangle(xc, yc, xb, yb, xa, ya)
       vert=2
       orient=-1
    if s2==long_side and s3==short_side:
       cosang=measureangle(xc, yc, xa, ya, xb, yb)
       vert=2
       orient=1

       
    return perimeter, ratio, cosang, vert, orient

    
def triangle(xa, ya, xb, yb, fa=None, fb=None, tolerance=3, maxratio=10):
    """Uses the triangle algorithm to match objects in the two coordinate
       to each other. From xyxymatch in iraf

    From iraf:

    The "triangles" algorithm  uses  a  sophisticated  pattern  matching
    technique  which  requires  no  tie point information from the user.
    It is expensive computationally and hence is restricted to a maximum
    of nmatch objects from the reference and input coordinate lists.
    
    The  "triangles"  algorithm  first  generates  a  list  of  all  the 
    possible triangles that can be formed from the points in each  list.
    For  a list of nmatch points this number is the combinatorial factor
    nmatch! / [(nmatch-3)! * 3!] or  nmatch * (nmatch-1) * (nmatch-2)  /
    6.   The length of the perimeter, ratio of longest to shortest side,
    cosine of the angle between  the  longest  and  shortest  side,  the
    tolerances  in  the  latter  two quantities and the direction of the
    arrangement of the  vertices  of  each  triangle  are  computed  and
    stored  in  a  table.   Triangles with vertices closer together than
    tolerance or with a ratio of the longest to  shortest  side  greater
    than  ratio  are  discarded.  The  remaining triangles are sorted in
    order of increasing ratio.  A sort merge algorithm is used to  match
    the   triangles   using   the  ratio  and  cosine  information,  the 
    tolerances in these quantities, and the maximum tolerances for  both
    lists.  Next  the  ratios of the perimeters of the matched triangles
    are  compared  to  the  average  ratio  for  the  entire  list,  and 
    triangles  which deviate too widely from the mean are discarded. The
    number of triangles remaining are  divided  into  the  number  which
    match  in  the  clockwise  sense  and the number which match int the
    counter-clockwise  sense.  Those  in  the  minority   category   are 
    eliminated.   The rejection step can be repeated up to nreject times
    or until no more rejections occur whichever comes first.   The  last
    step  in the algorithm is a voting procedure in which each remaining
    matched triangle casts three votes, one for  each  matched  pair  of
    vertices.   Points  which have fewer than half the maximum number of
    votes are discarded. The final set of matches  are  written  to  the
    output file.


    Parameters
    -----------
    xa : `~numpy.ndarray`
        x position for sources in image a
    ya : `~numpy.ndarray`
        y position for sources in image a
    xb : `~numpy.ndarray`
        x position for sources in image b
    yb : `~numpy.ndarray`
        y position for sources in image b
    fa : `~numpy.ndarray`
        Flux  for sources in image a.  Can be set to None.
    fb : `~numpy.ndarray`
        Flux  for sources in image b.  Can be set to None.

    """
    if fa is None: maga=xa*0.0+1.0
    if fb is None: magb=xa*0.0+1.0

    #create a list of all triangles in A-data set
    na=len(xa)
    perimeter_list=[]
    ratio_list=[]
    angle_list=[]
    vert_list=[]
    orient_list=[]
    cooa_list=[]
    for i in range(na):
      for j in range(i+1, na):
        for k in range(j+1, na):
           perimeter,ratio, angle, vert, orient=measuretriangle(xa[i], ya[i],  xa[j], ya[j], xa[k], ya[k])
           perimeter_list.append(perimeter)
           ratio_list.append(ratio)
           angle_list.append(angle)
           vert_list.append(vert)
           orient_list.append(orient)
           cooa_list.append([(xa[i], ya[i]),  (xa[j], ya[j]), (xa[k], ya[k])])
    pera_arr=np.array(perimeter_list)
    rata_arr=np.array(ratio_list)
    anga_arr=np.array(angle_list)
    vera_arr=np.array(vert_list)
    oria_arr=np.array(orient_list)

    #create a list of all triangles in A-data set
    nb=len(xb)
    perimeter_list=[]
    ratio_list=[]
    angle_list=[]
    vert_list=[]
    orient_list=[]
    coob_list=[]
    for i in range(nb):
      for j in range(i+1, nb):
        for k in range(j+1, nb):
           perimeter,ratio, angle, vert, orient=measuretriangle(xb[i], yb[i],  xb[j], yb[j], xb[k], yb[k])
           perimeter_list.append(perimeter)
           ratio_list.append(ratio)
           angle_list.append(angle)
           vert_list.append(vert)
           orient_list.append(orient)
           coob_list.append([(xb[i], yb[i]),  (xb[j], yb[j]), (xb[k], yb[k])])
    perb_arr=np.array(perimeter_list)
    ratb_arr=np.array(ratio_list)
    angb_arr=np.array(angle_list)
    verb_arr=np.array(vert_list)
    orib_arr=np.array(orient_list)

    #now match the two lists
    #TODO--make it be able to handle rotations/scale changes
    match_list=[]
    for i in range(len(cooa_list)):
        #if anga_arr[i]>tolerance and rata_arr[i]<maxratio:
        #   d=(rata_arr[i]-ratb_arr)**2/rata_arr[i]+(anga_arr[i]-angb_arr)**2/anga_arr[i]
        #   i_d=d.argsort()
        #   for j in range(len(coob_list)): 
        j=np.where((abs(perb_arr-pera_arr[i])<tolerance) )[0]
        if len(j)==1:
           j=j[0]
           inda_vert=get_vertices(vera_arr[i], oria_arr[i])
           indb_vert=get_vertices(verb_arr[i], orib_arr[i])
           for x,y in zip(inda_vert, indb_vert):
               match_list.append([cooa_list[i][x], coob_list[j][y]])
    return match_list
       
    
def get_vertices(vert, orient):
    vlist=[vert]
    for i in range(2):
        v=vlist[-1]+orient
        if v==3:  v=0
        if v==-1: v=2
        vlist.append(v)
    return vlist

