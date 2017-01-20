from math import *
import numpy as np
import sys
#sys.setrecursionlimit(1000)
import matplotlib.pyplot as plt
#from interval import interval, inf, imath

# function need to be optimized
#x = np.linspace(0,1,100) # 100 linearly spaced numbers
#y = (np.sin(13*x)*np.sin(27*x)+1)/2 # computing the values of sin(x)/x
def gland_function(x):
    y=x*(1-x)*(4-np.sqrt(abs(np.sin(60*x))))
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
    def __init__(self, center=None, depth=None,leftmin=None, rightmax=None, parent = None):

        self.parentNode = parent
        self.childNodes = []
        self.center=center
        self.depth=depth
        self.leftmin=leftmin
        self.rightmax=rightmax

    def Selectnode(self,node,function_evalution,h):
        leaves_of_depth=[]
        index=[]
        value=[]
        for i in range(len(node)):
            if node[i].depth==h:
                leaves_of_depth.append(node[i])
                index.append(i)

        for i in range(len(leaves_of_depth)):
            value.append(gland_function(leaves_of_depth[i].center))
            #print value
        new_index=np.argmax(value)
        #print new_index
        hi=node[index[new_index]].center
        fh=function(hi)
        #print hi
        #print fh

        function_evalution.append(fh)

        return node[index[new_index]],index[new_index],fh
    def Addnode(self,center_value,depth,leftmargin,rightmargin):
        n = Node(center=center_value,depth=depth, leftmin=leftmargin,rightmax=rightmargin, parent = self)
        self.childNodes.append(n)
        return n



def SOO():
    rootnode = Node(center=0.5, depth=0,leftmin=0, rightmax=1)
    current_node=[]
    current_node.append(rootnode)
    node=rootnode
    leaf=[]
    final=[]
    function_evalution=[]
    t=1
    h_tree=0
    ini_f=float("-inf")
    f_eva=0
    while t<=200:
        v_max=float("-inf")
        h_max=np.sqrt(t)
        h_tree=tree_depth(current_node)
        loop=mini(h_max,h_tree)
        h=0
        for h in range(0,t):
            #print h
            check_leaves=[]
            for i in range(len(current_node)):
                check_leaves.append(current_node[i].depth)
            if h in check_leaves:
                node,index,ini_f = node.Selectnode(current_node,function_evalution,h)
                print ini_f
                f_eva=f_eva+1
                current_node.pop(index)
                if ini_f>=v_max:
                    mid=midpoint(node.leftmin, node.rightmax)
                    left_center=midpoint(node.leftmin,mid)
                    right_center=midpoint(mid,node.rightmax)
                    current_node.append(node.Addnode(left_center,node.depth+1,node.leftmin,mid))
                    current_node.append(node.Addnode(right_center,node.depth+1,mid,node.rightmax))
                    t=t+1
                    v_max=ini_f
            h=h+1
            if f_eva>=150:
                break
        if f_eva>=150:
            break



    #for i in range(len(current_node)):
        #final.append(function(current_node[i].center)
    #findex=np.argmax(final)
    #print function_evalution
    #return current_node[findex].center



d=SOO()
#print d
