from math import *
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate
#import sys
#sys.setrecursionlimit(1000)
#import matplotlib.pyplot as plt
#from interval import interval, inf, imath

# function need to be optimized
#x = np.linspace(0,1,100) # 100 linearly spaced numbers
#y = (np.sin(13*x)*np.sin(27*x)+1)/2 # computing the values of sin(x)/x

def get_energy_ads(x,y,h):

    slab = fcc111('Cu', size=(4, 4, 2), vacuum=10.0)

    slab.set_calculator(EMT())
    e_slab = slab.get_potential_energy()

    amol = Atoms('N', positions=[(x, y, 0)])
    amol.set_calculator(EMT())
    e_N = amol.get_potential_energy()

    add_adsorbate(slab, amol, h, 'ontop')
    constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
    slab.set_constraint(constraint)
    #dyn = QuasiNewton(slab, trajectory='N2Cu.traj')
    #dyn.run(fmax=0.05)
    e_N_slab = slab.get_potential_energy()

    return  -(-e_slab - e_N +e_N_slab)
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


def split2(x1left,x1right,x2left,x2right):
    left_mid=x1left+x1right/2.0
    right_mid=x2left+x2right/2.0
    return left_mid,right_mid

def split3(x1left,x1right,x2left,x2right,x3left,x3right):
    x1mid=x1left+x1right/2.0
    x2mid=x2left+x2right/2.0
    x3mid=x3left+x3right/2.0

    return x1mid,x2mid,x3mid

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
    def __init__(self, x1center=None, x2center=None, x3center=None, depth=None,values=None, x1leftmin=None, x1rightmax=None,x2leftmin=None,
    x2rightmax=None, x3leftmin=None, x3rightmax=None,parent = None):

        self.parentNode = parent
        self.childNodes = []
        self.x1center=x1center
        self.x2center=x2center
        self.x3center=x3center
        self.depth=depth
        self.values=values
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
            #value.append(surface_function(leaves_of_depth[i].x1center, leaves_of_depth[i].x2center))
            value.append(get_energy_ads(leaves_of_depth[i].x1center,leaves_of_depth[i].x2center,leaves_of_depth[i].x3center))


            #print value
        new_index=np.argmax(value)
        #values.pop(new_index)
        #print new_index
        hl=node[index[new_index]].x1center
        hm=node[index[new_index]].x2center
        hr=node[index[new_index]].x3center
        #fh=Branin_function(hl,hr)
        #fh=surface_function(hl,hr)
        fh=get_energy_ads(hl,hm,hr)

        function_evalution.append(fh)
        #print hl
        #print hr
        print fh

        #function_evalution.append(hi)

        return node[index[new_index]],index[new_index],fh,function_evalution
    def Addnode(self,x1center,x2center,x3center, depth,x1leftmin,x1rightmax,x2leftmin,x2rightmax ,x3leftmin,x3rightmax):
        n = Node(x1center=x1center,x2center=x2center,x3center=x3center,depth=depth, x1leftmin=x1leftmin,x1rightmax=x1rightmax,
        x2leftmin=x2leftmin,x2rightmax=x2rightmax, x3leftmin=x3leftmin,x3rightmax=x3rightmax,parent = self)
        self.childNodes.append(n)
        return n



def SOO():
    rootnode = Node(x1center=1.5,x2center=1.5, x3center=2.5, depth=0,x1leftmin=0, x1rightmax=3,x2leftmin=0,x2rightmax=3,
    x3leftmin=1,x3rightmax=4)
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
    node_values=[]
    node_values.append([])
    while t<=20000:
        v_max=float("-inf")
        h_max=t
        h_tree=tree_depth(current_node)
        loop=mini(h_max,h_tree)
        h=0
        while h<=loop:
            #print h
            check_leaves=[]
            leaves_value=[]
            for i in range(len(current_node)):
                check_leaves.append(current_node[i].depth)
            for i in range(len(current_node)):
                if h in check_leaves:
                    leaves_value.append(current_node[i])
            if h in check_leaves:
                node,index,ini_f,ggg = node.Selectnode(leaves_value,h)
                final_max.append(ggg)

                current_node.pop(index)
                if ini_f>=v_max:
                    #x1mid,x2mid=split2(node.x1leftmin, node.x1rightmax,node.x2leftmin,node.x2rightmax)
                    x1mid,x2mid,x3mid=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax)
                    centerx11,centerx12,centerx13=split3(node.x1leftmin,x1mid,node.x2leftmin,x2mid,node.x3leftmin,x3mid)
                    centerx21,centerx22,centerx23=split3(node.x1leftmin,x1mid,node.x2leftmin,x2mid,x3mid,node.x3rightmax)
                    centerx31,centerx32,centerx33=split3(node.x1leftmin,x1mid,x2mid,node.x2rightmax,node.x3leftmin,x3mid)
                    centerx41,centerx42,centerx43=split3(node.x1leftmin,x1mid,x2mid,node.x2rightmax,x3mid,node.x3rightmax)
                    centerx51,centerx52,centerx53=split3(x1mid,node.x1leftmin,node.x2leftmin,x2mid,node.x3leftmin,x3mid)
                    centerx61,centerx62,centerx63=split3(x1mid,node.x1leftmin,node.x2leftmin,x2mid,x3mid,node.x3rightmax)
                    centerx71,centerx72,centerx73=split3(x1mid,node.x1leftmin,x2mid,node.x2rightmax,node.x3leftmin,x3mid)
                    centerx81,centerx82,centerx83=split3(x1mid,node.x1leftmin,x2mid,node.x2rightmax,x3mid,node.x3rightmax)
                    f_eva=f_eva+8




                    current_node.append(node.Addnode(centerx11,centerx12, centerx13, node.depth+1,node.x1leftmin,x1mid,node.x2leftmin,x2mid,node.x3leftmin,x3mid))
                    current_node.append(node.Addnode(centerx21,centerx22, centerx23, node.depth+1,node.x1leftmin,x1mid,node.x2leftmin,x2mid,x3mid,node.x3rightmax))
                    current_node.append(node.Addnode(centerx31,centerx32, centerx33, node.depth+1,node.x1leftmin,x1mid,x2mid,node.x2rightmax,node.x3leftmin,x3mid))
                    current_node.append(node.Addnode(centerx41,centerx42, centerx43, node.depth+1,node.x1leftmin,x1mid,x2mid,node.x2rightmax,x3mid,node.x3rightmax))
                    current_node.append(node.Addnode(centerx51,centerx52, centerx53, node.depth+1,x1mid,node.x1leftmin,node.x2leftmin,x2mid,node.x3leftmin,x3mid))
                    current_node.append(node.Addnode(centerx61,centerx62, centerx63, node.depth+1,x1mid,node.x1leftmin,node.x2leftmin,x2mid,x3mid,node.x3rightmax))
                    current_node.append(node.Addnode(centerx71,centerx72, centerx73, node.depth+1,x1mid,node.x1leftmin,x2mid,node.x2rightmax,node.x3leftmin,x3mid))
                    current_node.append(node.Addnode(centerx81,centerx82, centerx83, node.depth+1,x1mid,node.x1leftmin,x2mid,node.x2rightmax,x3mid,node.x3rightmax))



                    t=t+1
                    v_max=ini_f
            h=h+1
            if f_eva>=300:
                break
        if f_eva>=300:
            break

    kl=np.argmax(final_max)

    return final_max[kl]



d=SOO()
print d
