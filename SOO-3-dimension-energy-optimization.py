from math import *
import numpy as np
import sys
#sys.setrecursionlimit(1000)
import matplotlib.pyplot as plt
#from interval import interval, inf, imath

# function need to be optimized
#x = np.linspace(0,1,100) # 100 linearly spaced numbers
#y = (np.sin(13*x)*np.sin(27*x)+1)/2 # computing the values of sin(x)/x
def surface_function(x1,x2):
    y=-(x1**2+x2**2)
    return y
def Branin_function(x1,x2):
    y=-((x2-(5.1/(4*(np.pi)**2))*x1**2+(5/np.pi)*x1-6)**2+10*(1-1/(8*np.pi))*np.cos(x1)+10)
    return y

def function(x):
    y = (np.sin(13*x)*np.sin(27*x)+1)/2.0
    return y
def mini(x,y):
    if x<=y:
        return x
    else:
        return y


def split(x1left,x1right,x2left,x2right):
    left_mid=x1left+x1right/2.0
    right_mid=x2left+x2right/2.0
    return left_mid,right_mid
def midpoint(p1, p2):
    mid=(p1+p2)/2.0

    return mid

def tree_depth(node):
    dep=[]
    for i in range(len(node)):
        dep.append(node[i].depth)
    dindex=np.argmax(dep)

    return dep[dindex]
class Node:
    def __init__(self, x1center=None, x2center=None, x3center=None, depth=None,x1leftmin=None, x1rightmax=None,x2leftmin=None,
    x2rightmax=None, x3leftmin=None, x3rightmax=None,parent = None):

        self.parentNode = parent
        self.childNodes = []
        self.x1center=x1center
        self.x2center=x2center
        self.x3center=x3center
        self.depth=depth
        self.x1leftmin=x1leftmin
        self.x1rightmax=x1rightmax
        self.x2leftmin=x2leftmin
        self.x2rightmax=x2rightmax
        self.x3leftmin=x3leftmin
        self.x3rightmax=x3rightmax

    def Selectnode(self,node,h):
        leaves_of_depth=[]
        index=[]
        function_evalution=[]
        value=[]
        for i in range(len(node)):
            if node[i].depth==h:
                leaves_of_depth.append(node[i])
                index.append(i)

        for i in range(len(leaves_of_depth)):
            #value.append(Branin_function(leaves_of_depth[i].x1center, leaves_of_depth[i].x2center))
            value.append(surface_function(leaves_of_depth[i].x1center, leaves_of_depth[i].x2center))

            #print value
        new_index=np.argmax(value)
        #print new_index
        hl=node[index[new_index]].x1center
        hr=node[index[new_index]].x2center
        #fh=Branin_function(hl,hr)
        fh=surface_function(hl,hr)

        function_evalution.append(fh)
        #print hl
        #print hr
        print fh

        #function_evalution.append(hi)

        return node[index[new_index]],index[new_index],fh,function_evalution
    def Addnode(self,x1center,x2center,x3center, depth,x1leftmin,x1rightmax,x2leftmin,x2rightmax):
        n = Node(x1center=x1center,x2center=x2center,depth=depth, x1leftmin=x1leftmin,x1rightmax=x1rightmax, x2leftmin=x2leftmin,x2rightmax=x2rightmax, parent = self)
        self.childNodes.append(n)
        return n



def SOO():
    rootnode = Node(x1center=1.5,x2center=1.5,x3center, depth=0,x1leftmin=-1, x1rightmax=2,x2leftmin=-3,x2rightmax=1)
    current_node=[]
    current_node.append(rootnode)
    node=rootnode
    leaf=[]
    final=[]
    final_max=[]
    t=1
    h_tree=0
    ini_f=float("-inf")
    f_eva=0
    while t<=20000:
        v_max=float("-inf")
        h_max=t
        h_tree=tree_depth(current_node)
        loop=mini(h_max,h_tree)
        h=0
        while h<=loop:
            print h
            check_leaves=[]
            for i in range(len(current_node)):
                check_leaves.append(current_node[i].depth)
            if h in check_leaves:
                node,index,ini_f,ggg= node.Selectnode(current_node,h)
                final_max.append(ggg)
                f_eva=f_eva+1
                current_node.pop(index)
                if ini_f>=v_max:
                    x1mid,x2mid=split(node.x1leftmin, node.x1rightmax,node.x2leftmin,node.x2rightmax)
                    centerx11,centerx12=split(node.x1leftmin,x1mid,node.x2leftmin,x2mid)
                    centerx21,centerx22=split(node.x1leftmin,x1mid,x2mid,node.x2rightmax)
                    centerx31,centerx32=split(x1mid,node.x1rightmax,node.x2leftmin,x2mid)
                    centerx41,centerx42=split(x1mid,node.x1rightmax,x2mid,node.x2rightmax)

                    current_node.append(node.Addnode(centerx11,centerx12, node.depth+1,node.x1leftmin,x1mid,node.x2leftmin,x2mid))
                    current_node.append(node.Addnode(centerx21,centerx22, node.depth+1,node.x1leftmin,x1mid,x2mid,node.x2rightmax))
                    current_node.append(node.Addnode(centerx31,centerx32, node.depth+1,x1mid,node.x1rightmax,node.x2leftmin,x2mid))
                    current_node.append(node.Addnode(centerx41,centerx42, node.depth+1,x1mid,node.x1rightmax,x2mid,node.x2rightmax))

                    t=t+1
                    v_max=ini_f
            h=h+1
            if f_eva>=1000:
                break
        if f_eva>=1000:
            break

    kl=np.argmax(final_max)

    return final_max[kl]



d=SOO()
print d
