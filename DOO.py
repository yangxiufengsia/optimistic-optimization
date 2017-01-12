from math import *
import numpy as np
import sys
#sys.setrecursionlimit(1000)
import matplotlib.pyplot as plt
#from interval import interval, inf, imath

# function need to be optimized
#x = np.linspace(0,1,100) # 100 linearly spaced numbers
#y = (np.sin(13*x)*np.sin(27*x)+1)/2 # computing the values of sin(x)/x
def function(x):
    y = (np.sin(13*x)*np.sin(27*x)+1)/2.0
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
    def __init__(self, center=None, depth=None,leftmin=None, rightmax=None, parent = None):

        self.parentNode = parent
        self.childNodes = []
        self.center=center
        self.depth=depth
        self.leftmin=leftmin
        self.rightmax=rightmax

    def Selectnode(self,node,function_evalution):
        #for i in range(len(node)):
        #print node[0].center
        #print node[0].depth
        #print len(node)
        #s = sorted(self.childNodes, key = lambda c: c.wins/c.visits + 0.5*sqrt(2*log(self.visits)/c.visits))[-1]
        value=[]

        for i in range(len(node)):
            value.append(function(node[i].center)+222*2**(-2*node[i].depth))
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
            #s=random.choice(node[0],node[1])
        #print s.depth
        function_evalution.append(node[index].center)

        return node[index],index
    def Addnode(self,center_value,depth,leftmargin,rightmargin):
        n = Node(center=center_value,depth=depth, leftmin=leftmargin,rightmax=rightmargin, parent = self)
        self.childNodes.append(n)
        #print n.center
        return n


def DOO():
    rootnode = Node(depth=0,leftmin=0, rightmax=1)
    current_node=[]
    node=rootnode
    leaf=[]
    final=[]
    function_evalution=[]

    for i in range(150):


        if node.childNodes!=[]:

            #select the node
            node ,index= node.Selectnode(current_node,function_evalution)
            current_node.pop(index)
        #spliting search space into 2-ary
        #calculate the center of the search space
        mid=midpoint(node.leftmin, node.rightmax)
        print mid
        #calculate two child center
        left_center=midpoint(node.leftmin,mid)
        right_center=midpoint(mid,node.rightmax)
        #print left_child
        #print right_child
        print node.depth

        #expand two nodes at one time
        current_node.append(node.Addnode(left_center,node.depth+1,node.leftmin,mid))
        current_node.append(node.Addnode(right_center,node.depth+1,mid,node.rightmax))

        #print current_node[0]
    for i in range(len(current_node)):
        final.append(function(current_node[i].center)+222*2**(-2*current_node[i].depth))
    findex=np.argmax(final)
    print function_evalution
    return current_node[findex].center



d=DOO()
print d


"""plot the final result"""
u=[]
x=[0.25, 0.75, 0.875, 0.375, 0.125, 0.625, 0.0625, 0.5625, 0.4375, 0.8125, 0.9375, 0.6875, 0.1875, 0.3125, 0.40625, 0.84375, 0.53125, 0.09375, 0.90625, 0.03125, 0.71875, 0.96875, 0.46875, 0.21875, 0.34375, 0.59375, 0.859375, 0.390625, 0.28125, 0.78125, 0.890625, 0.078125, 0.421875, 0.8671875, 0.65625, 0.546875, 0.046875, 0.3984375, 0.87109375, 0.86328125, 0.8828125, 0.869140625, 0.865234375, 0.8681640625, 0.8662109375, 0.86767578125, 0.86669921875, 0.867431640625, 0.867919921875, 0.8675537109375, 0.86749267578125, 0.8673095703125, 0.86761474609375, 0.867523193359375, 0.8675384521484375, 0.8675079345703125, 0.867584228515625, 0.8675308227539062, 0.8675155639648438, 0.8675270080566406, 0.8675251007080078, 0.8675193786621094, 0.8675289154052734, 0.8675260543823242, 0.8675265312194824, 0.867527961730957, 0.867525577545166, 0.8675262928009033, 0.8675258159637451, 0.8675261735916138, 0.8675264120101929, 0.8675262331962585, 0.867526113986969, 0.8675262033939362, 0.8675262182950974, 0.8675262629985809, 0.867526188492775, 0.8675262108445168, 0.8675261959433556, 0.8675262071192265, 0.867526214569807, 0.8675262089818716, 0.8675262052565813, 0.867526208050549, 0.8675262099131942, 0.8675262085162103, 0.8675262075848877, 0.8675262094475329, 0.8675262082833797, 0.867526208749041, 0.8675262078177184, 0.8675262081669644, 0.867526208399795, 0.8675262086326256, 0.8675262088654563, 0.8675262077013031, 0.8675262079341337, 0.8675262081087567, 0.867526208225172, 0.8675262083415873, 0.8675262084580027, 0.867526208574418, 0.8675262086908333, 0.8675262088072486, 0.8675262089236639, 0.8675262076430954, 0.8675262077595107, 0.867526207875926, 0.8675262079923414, 0.8675262080796529, 0.8675262081378605, 0.8675262081960682, 0.8675262082542758, 0.8675262083124835, 0.8675262083706912, 0.8675262084288988, 0.8675262084871065, 0.8675262085453141, 0.8675262086035218, 0.8675262086617295, 0.8675262087199371, 0.8675262087781448, 0.8675262088363525, 0.8675262088945601, 0.8675262076139916, 0.8675262076721992, 0.8675262077304069, 0.8675262077886146, 0.8675262078468222, 0.8675262079050299, 0.8675262079632375, 0.8675262080214452, 0.867526208065101, 0.8675262080942048, 0.8675262081233086, 0.8675262081524124, 0.8675262081815163, 0.8675262082106201, 0.8675262082397239, 0.8675262082688278, 0.8675262082979316, 0.8675262083270354, 0.8675262083561393, 0.8675262083852431, 0.8675262084143469, 0.8675262084434507, 0.8675262084725546, 0.8675262085016584, 0.8675262085307622]
print len(x)

for i in range(len(x)):
    f=(np.sin(13*x[i])*np.sin(27*x[i])+1)/2.0
    u.append(f)
print len(u)
plt.scatter(x,u)
t=np.linspace(0,1,200)
y=(np.sin(13*t)*np.sin(27*t)+1)/2.0
plt.plot(t,y,'r--')
plt.show()
