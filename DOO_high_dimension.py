from math import *
import numpy as np

from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
#from interval import interval, inf, imath

# function need to be optimized
#x = np.linspace(0,1,100) # 100 linearly spaced numbers
#y = (np.sin(13*x)*np.sin(27*x)+1)/2 # computing the values of sin(x
#x1=np.linspace(-5,10,100)
#x2=np.linspace(0,15,100)
#x1=9.42478
#x2=2.475
#print x1
#y=(x2-(5.1/(4*(np.pi)**2))*x1**2+(5/np.pi)*x1-6)**2+10*(1-1/(8*np.pi))*np.cos(x1)+10
#print y
def function(x1,x2):
    #y = (np.sin(13*x)*np.sin(27*x)+1)/2.0
    y1= (np.sin(13*x1)*np.sin(27*x1)+1)/2.0
    y2= (np.sin(13*x2)*np.sin(27*x2)+1)/2.0
    y=y1*y2


    #y=-((x2-(5.1/(4*(np.pi)**2))*x1**2+(5/np.pi)*x1-6)**2+10*(1-1/(8*np.pi))*np.cos(x1)+10)
    return y


#delta=14*2**(-h)
#cen=interval([0,1]).midpoint
#print cen
#test=np.sin(13*cen)
#print test
def split(x1left,x1right,x2left,x2right):
    left_mid=x1left+x1right/2.0
    right_mid=x2left+x2right/2.0
def midpoint(p1, p2):
    mid=(p1+p2)/2.0

    return mid

class Node:
    def __init__(self, x1center=None, x2center=None, depth=None,x1leftmin=None, x1rightmax=None, x2leftmin=None, x2rightmax=None,parent = None):

        self.parentNode = parent
        self.childNodes = []
        self.x1center=x1center
        self.x2center=x2center
        self.depth=depth
        self.x1leftmin=x1leftmin
        self.x1rightmax=x1rightmax
        self.x2leftmin=x1leftmin
        self.x2rightmax=x2rightmax

    def Selectnode(self,node,function_evalution):
        #for i in range(len(node)):
        #print node[0].center
        #print node[0].depth
        #print len(node)
        #s = sorted(self.childNodes, key = lambda c: c.wins/c.visits + 0.5*sqrt(2*log(self.visits)/c.visits))[-1]
        value=[]

        for i in range(len(node)):
            value.append(function(node[i].x1center,node[i].x2center)+11*2**(-2*node[i].depth))
            #print value[i]
        #child1=function(node[0].center)+14*2**(-node[0].depth)
        #child2=function(node[1].center)+14*2**(-node[1].depth)
        index=np.argmax(value)
        #print index
        #s = sorted(self.childNodes, key = lambda c: c.wins/c.visits + 0.5*sqrt(2*log(self.visits)/c.visits))[-1]
        #if child1>child2:
            #s=node[0]
        #if child1<child2:
            #s=node[1]
        #if child1==child2:
        hi=function(node[index].x1center,node[index].x2center)
            #s=random.choice(node[0],node[1])
        #print s.depth
        #function_evalution.append(node[index].x1center)
        print node[index].x1center,node[index].x2center
        print hi
        #print node[index].depth

        return node[index],index
    def Addnode(self,x1child_center,x2child_center, depth,x1leftmargin,x1rightmargin,x2leftmargin,x2rightmargin):
        n = Node(x1center=x1child_center,x2center=x2child_center, depth=depth, x1leftmin=x1leftmargin,x1rightmax=x1rightmargin, x2leftmin=x2leftmargin,x2rightmax=x2rightmargin,parent = self)
        self.childNodes.append(n)
        #print n.center
        return n


def DOO():
    rootnode = Node(depth=0,x1leftmin=0, x1rightmax=1,x2leftmin=0,x2rightmax=1)
    current_node=[]
    node=rootnode
    leaf=[]
    final=[]
    function_evalution=[]

    for i in range(150):


        if node.childNodes!=[]:

            #select the node
            node,index= node.Selectnode(current_node,function_evalution)
            current_node.pop(index)
        #spliting search space into 2-ary
        #calculate the center of the search space
        if node.depth%2==0:
            mid=midpoint(node.x1leftmin, node.x1rightmax)
            #print "vertival:", mid
            #x1left=midpoint(node.x1leftmin,mid)
            x1left_child_center=midpoint(node.x1leftmin,mid)
            x1right_child_center=midpoint(mid, node.x1rightmax)
            x2left_child_center=midpoint(node.x2leftmin, node.x2rightmax)
            x2right_child_center=x2left_child_center


            current_node.append(node.Addnode(x1left_child_center,x2left_child_center, node.depth+1, node.x1leftmin, mid, node.x2leftmin,node.x2rightmax))
            current_node.append(node.Addnode(x1right_child_center,x2right_child_center, node.depth+1, mid, node.x1rightmax, node.x2leftmin,node.x2rightmax))


        else:
            mid=midpoint(node.x2leftmin, node.x2rightmax)
            #print "horizaon:",mid
            x2left_child_center=midpoint(node.x2leftmin,mid)
            x2right_child_center=midpoint(mid, node.x2rightmax)
            x1left_child_center=midpoint(node.x1leftmin, node.x1rightmax)
            x1right_child_center=x1left_child_center


            current_node.append(node.Addnode(x1left_child_center,x2left_child_center, node.depth+1, node.x1leftmin, node.x1rightmax, node.x2leftmin,mid))
            current_node.append(node.Addnode(x1right_child_center,x2right_child_center, node.depth+1, node.x1leftmin, node.x1rightmax, mid,node.x2rightmax))

        #print mid


        #calculate two child center
        #left_center=midpoint(node.leftmin,mid)
        #right_center=midpoint(mid,node.rightmax)
        #print left_child
        #print right_child
        #print node.depth

        #expand two nodes at one time
        #current_node.append(node.Addnode(x1left_center,node.depth+1,node.leftmin,mid))
        #current_node.append(node.Addnode(right_center,node.depth+1,mid,node.rightmax))

        #print current_node[0]
    for i in range(len(current_node)):
        final.append(function(current_node[i].x1center,node.x2center)+11*2**(-2*current_node[i].depth))
    findex=np.argmax(final)
    #print function_evalution
    return current_node[findex].x1center,current_node[findex].x2center



b,d=DOO()
#print b,d
