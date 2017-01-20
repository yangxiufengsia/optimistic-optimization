from math import *
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate
from math import floor
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

    add_adsorbate(slab, amol, h, position=(x,y))
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
    x1mid=(x1left+x1right)/2.0
    x2mid=(x2left+x2right)/2.0
    x3mid=(x3left+x3right)/2.0

    return x1mid,x2mid,x3mid
def split8(x,y):
    element=(y-x)/3.0
    return element
def small_interval(left,right):
    element=(left+right)/8.0
    #e=[]
    #e.append(left+element)
    #for i in range(8):
        #e.append(left+element)



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
    def __init__(self, values=None, x1center=None, x2center=None, x3center=None, depth=None, x1leftmin=None, x1rightmax=None,x2leftmin=None,
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
            value.append(leaves_of_depth[i].values)
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
        print "left_center:",hl
        print "mid_center:",hm
        print "right_center:",hr
        print fh
        #print node[index[new_index]].depth
        #function_evalution.append(hi)

        return node[index[new_index]],index[new_index],fh,function_evalution
    def Addnode(self,values,x1center,x2center,x3center, depth,x1leftmin,x1rightmax,x2leftmin,x2rightmax ,x3leftmin,x3rightmax):
        n = Node(values=values,x1center=x1center,x2center=x2center,x3center=x3center,depth=depth, x1leftmin=x1leftmin,x1rightmax=x1rightmax,
        x2leftmin=x2leftmin,x2rightmax=x2rightmax, x3leftmin=x3leftmin,x3rightmax=x3rightmax,parent = self)
        self.childNodes.append(n)
        return n



def SOO():
    #for i in range(1000)
    rootnode = Node(values=get_energy_ads(1.5,1.5,2.5),x1center=1.5,x2center=1.5, x3center=2.5, depth=1,x1leftmin=0.0, x1rightmax=3.0,x2leftmin=0.0,x2rightmax=3.0,
    x3leftmin=1.0,x3rightmax=4.0)
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
    """initialization of the function values in each depth"""
    for i in range(1000):
        node_values.append([])
    while t<=20000:
        v_max=float("-inf")
        h_max=floor(np.sqrt(t))
        h_tree=tree_depth(current_node)
        loop=mini(h_max,h_tree)
        h=1
        for h in range(0,int(floor(loop))):
            #print h
            check_leaves=[]
            leaves_value=[]
            for i in range(len(current_node)):
                check_leaves.append(current_node[i].depth)
            #for i in range(len(current_node)):
                #if h in check_leaves:
                    #leaves_value.append(current_node[i])
            #print check_leaves
            if h in check_leaves:
                node,index,ini_f,ggg = node.Selectnode(current_node,h)
                final_max.append(ggg)

                current_node.pop(index)
                if ini_f>=v_max:
                    most_width=[]
                    widthx1=node.x1rightmax-node.x1leftmin
                    most_width.append(widthx1)
                    widthx2=node.x2rightmax-node.x2leftmin
                    most_width.append(widthx2)
                    widthx3=node.x3rightmax-node.x3leftmin
                    #print widthx1,widthx2,widthx3
                    most_width.append(widthx3)
                    width_index=np.argmax(most_width)
                    if width_index==0:
                        ele=split8(node.x1leftmin,node.x1rightmax)
                        e2=node.x1leftmin+ele
                        e3=e2+ele

                        centerx11,centerx12,centerx13=split3(node.x1leftmin,e2,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax)
                        centerx21,centerx22,centerx23=split3(e2,e3,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax)
                        centerx31,centerx32,centerx33=split3(e3,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax)

                        value1=get_energy_ads(centerx11,centerx12,centerx13)
                        value2=get_energy_ads(centerx21,centerx22,centerx23)
                        value3=get_energy_ads(centerx31,centerx32,centerx33)
                        f_eva=f_eva+3

                        current_node.append(node.Addnode(value1,centerx11,centerx12, centerx13, node.depth+1,node.x1leftmin,e2,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                        current_node.append(node.Addnode(value2,centerx21,centerx22, centerx23, node.depth+1,e2,e3,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                        current_node.append(node.Addnode(value3,centerx31,centerx32, centerx33, node.depth+1,e3,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))


                    if width_index==2:
                        ele=split8(node.x3leftmin,node.x3rightmax)
                        e2=node.x3leftmin+ele
                        e3=e2+ele

                        centerx11,centerx12,centerx13=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,e2)
                        centerx21,centerx22,centerx23=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e2,e3)
                        centerx31,centerx32,centerx33=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e3,node.x3rightmax)


                        value1=get_energy_ads(centerx11,centerx12,centerx13)
                        value2=get_energy_ads(centerx21,centerx22,centerx23)
                        value3=get_energy_ads(centerx31,centerx32,centerx33)

                        f_eva=f_eva+3


                        current_node.append(node.Addnode(value1,centerx11,centerx12, centerx13, node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,e2))
                        current_node.append(node.Addnode(value2,centerx21,centerx22, centerx23, node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e2,e3))
                        current_node.append(node.Addnode(value3,centerx31,centerx32, centerx33, node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e3,node.x3rightmax))





                    if width_index==1:
                        ele=split8(node.x2leftmin,node.x2rightmax)
                        e2=node.x2leftmin+ele
                        e3=e2+ele

                        centerx11,centerx12,centerx13=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,e2,node.x3leftmin,node.x3rightmax)
                        centerx21,centerx22,centerx23=split3(node.x1leftmin,node.x1rightmax,e2,e3,node.x3leftmin,node.x3rightmax)
                        centerx31,centerx32,centerx33=split3(node.x1leftmin,node.x1rightmax,e3,node.x2rightmax,node.x3leftmin,node.x3rightmax)


                        value1=get_energy_ads(centerx11,centerx12,centerx13)
                        value2=get_energy_ads(centerx21,centerx22,centerx23)
                        value3=get_energy_ads(centerx31,centerx32,centerx33)

                        f_eva=f_eva+3


                        current_node.append(node.Addnode(value1,centerx11,centerx12, centerx13, node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,e2,node.x3leftmin,node.x3rightmax))
                        current_node.append(node.Addnode(value2,centerx21,centerx22, centerx23, node.depth+1,node.x1leftmin,node.x1rightmax,e2,e3,node.x3leftmin,node.x3rightmax))
                        current_node.append(node.Addnode(value3,centerx31,centerx32, centerx33, node.depth+1,node.x1leftmin,node.x1rightmax,e3,node.x2rightmax,node.x3leftmin,node.x3rightmax))

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
