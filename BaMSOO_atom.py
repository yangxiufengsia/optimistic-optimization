from math import *
import numpy as np
import sys
#sys.setrecursionlimit(1000)
#import matplotlib.pyplot as plt
from numpy import linalg as LA
from numpy.linalg import inv
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate
from math import floor


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

ggfh=get_energy_ads(2.5,1.5,1.5)
print ggfh
#from interval import interval, inf, imath

# function need to be optimized
#x = np.linspace(0,1,100) # 100 linearly spaced numbers
#y = (np.sin(13*x)*np.sin(27*x)+1)/2 # computing the values of sin(x)/x
def energy_calculation():
    pass

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

def split8(x,y):
    element=(y-x)/3.0
    return element
def split(x1left,x1right,x2left,x2right):
    left_mid=(x1left+x1right)/2.0
    right_mid=(x2left+x2right)/2.0
    return left_mid,right_mid
def midpoint(p1, p2):
    mid=(p1+p2)/2.0

    return mid

def split3(x1left,x1right,x2left,x2right,x3left,x3right):
    x1mid=(x1left+x1right)/2.0
    x2mid=(x2left+x2right)/2.0
    x3mid=(x3left+x3right)/2.0

    return x1mid,x2mid,x3mid

def tree_depth(node):
    dep=[]
    #print "node len:",len(node)
    for i in range(len(node)):
        dep.append(node[i].depth)
    dindex=np.argmax(dep)

    return dep[dindex]

def Branin_function(x1,x2):
    y=-((x2-(5.1/(4*(np.pi)**2))*x1**2+(5/np.pi)*x1-6)**2+10*(1-1/(8*np.pi))*np.cos(x1)+10)
    return y



def kernel_calculation(idtr,xtr,xte):
    ktr=np.zeros((idtr,idtr))
    ktr_new=np.matrix(ktr)
    kte=np.zeros((1,idtr))
    kte_new=np.matrix(kte)
    #kte_train=np.zeros((idtr,idtr))
    #kte_test=np.zeros((idtr_l,idtr_l))
    #kte_kte=np.zeros((idtr,idtr_l))
    ### training kernel matrix calculation
    for i in range(0,idtr):
        for j in range(i,idtr):
            #ktr_new[i,j]=exp(-0.5*0.03*LA.norm(xtr[:,i]-xtr[:,j])**2)
            ktr_new[i,j]=0.5*(1+np.sqrt((5*LA.norm(xtr[:,i]-xtr[:,j])**2)/0.05)+
            5*LA.norm(xtr[:,i]-xtr[:,j])**2/(3*0.05))*np.exp(-(np.sqrt(5*LA.norm(xtr[:,i]-xtr[:,j])**2)/0.05))
            #print ktr[i,j]
            #ktr_new[j,i]=ktr_new[i,j]

    for i in range(1):
        for j in range(0,idtr):
            #kte_new[i,j]=exp(-0.5*0.03*(xte[:,i]-xtr[:,j])**2)
            kte_new[i,j]=0.5*(1+np.sqrt((5*LA.norm(xte[:,i]-xtr[:,j])**2)/0.05)+
            5*LA.norm(xte[:,i]-xtr[:,j])**2/(3*0.05))*np.exp(-np.sqrt((5*LA.norm(xte[:,i]-xtr[:,j])**2)/0.05))

    #print kte_new.shape

    return ktr_new,kte_new




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
        value=[]
        for i in range(len(node)):
            if node[i].depth==h:
                leaves_of_depth.append(node[i])
                index.append(i)

        for i in range(len(leaves_of_depth)):
            value.append(leaves_of_depth[i].values)
            #print value
        new_index=np.argmax(value)
        hl=node[index[new_index]].x1center
        hm=node[index[new_index]].x2center
        hr=node[index[new_index]].x3center

        #print "left_center:",hl
        #print "mid_center:,",hm
        #print "right_center:",hr
        fh=node[index[new_index]].values
        #print "function value:",fh
        #gl=function(node[index[new_index]].center)
        #print gl


        return node[index[new_index]],index[new_index],fh
    def Addnode(self,values,x1center,x2center,x3center, depth,x1leftmin,x1rightmax,x2leftmin,x2rightmax ,x3leftmin,x3rightmax):
        n = Node(values=values,x1center=x1center,x2center=x2center,x3center=x3center,depth=depth, x1leftmin=x1leftmin,x1rightmax=x1rightmax,
        x2leftmin=x2leftmin,x2rightmax=x2rightmax, x3leftmin=x3leftmin,x3rightmax=x3rightmax,parent = self)
        self.childNodes.append(n)
        return n



def SOO():
    """set g(0,0)=f(x(0,0))"""
    ini_f=float("-inf")
    g_function=get_energy_ads(1.5,1.5,2.5)

    """set f^=g(0,0)"""
    f_t=g_function

    """initialize the tree"""
    rootnode = Node(values=get_energy_ads(1.5,1.5,2.5),x1center=1.5,x2center=1.5, x3center=2.5, depth=0,x1leftmin=0.0, x1rightmax=3.0,x2leftmin=0.0,x2rightmax=3.0,
    x3leftmin=1.0,x3rightmax=4.0)


    """set t=1,n=1,N=1 and D_t={(x00,g00)}"""
    t=1 ## represent the number of sampled data
    n=1 ## represent the iteration of the algorithm
    N=1 ## represent the number of GP calculation
    D_x_train=np.array([1.5,1.5,2.5])
    D_y_train=np.array([g_function])

    current_node=[]
    current_node.append(rootnode)
    #print current_node[0].center
    node=rootnode
    leaf=[]
    final=[]
    function_evalution=[]
    final_result=[]
    t=1
    h_tree=0
    #ini_f=float("-inf")
    f_eva=0
    while n<=2000:
        print n
        v_max=float("-inf")
        h_max=np.sqrt(n)
        h_tree=tree_depth(current_node)
        loop=mini(h_max,h_tree)
        #print loop
        h=0
        for h in range(0,n):
            #print h
            #print len(current_node)
            #print "tree depth:", h
            #print len(current_node)
            check_leaves=[]
            for i in range(len(current_node)):
                check_leaves.append(current_node[i].depth)
            #print len(check_leaves)
            if h in check_leaves:
                node,index,g_value = node.Selectnode(current_node,h)
                #print "g_value:",g_value
                #f_eva=f_eva+1
                current_node.pop(index)
                #print len(current_node)
                if g_value>=v_max:
                    center=[]
                    most_width=[]
                    widthx1=node.x1rightmax-node.x1leftmin
                    most_width.append(widthx1)
                    widthx2=node.x2rightmax-node.x2leftmin
                    most_width.append(widthx2)
                    widthx3=node.x3rightmax-node.x3leftmin
                    #print widthx1,widthx2,widthx3
                    most_width.append(widthx3)
                    width_index=np.argmax(most_width)
                    #print width_index
                    if width_index==2:
                        ele=split8(node.x3leftmin,node.x3rightmax)
                        e2=node.x3leftmin+ele
                        e3=e2+ele

                        centerx11,centerx12,centerx13=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,e2)
                        centerx21,centerx22,centerx23=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e2,e3)
                        centerx31,centerx32,centerx33=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e3,node.x3rightmax)


                        leftnode=[centerx11,centerx12,centerx13]
                        midnode=[centerx21,centerx22,centerx23]
                        rightnode=[centerx31,centerx32,centerx33]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)
                        #print "i=2:",center[0][0],center[0][1],center[0][2]

                        #print "i=2:",center[1][0],center[1][1],center[1][2]

                        #print "i=2:",center[2][0],center[2][1],center[2][2]




                        #node_function_value=[]

                        for i in range(3):


                            N=N+1
                            D_x_test=np.array([center[i]])
                            ktrain,ktest=kernel_calculation(t,np.matrix(D_x_train),np.matrix(D_x_test))
                            u_N=ktest*inv(ktrain)*np.matrix(D_y_train).T
                            #print "u_N:",u_N
                            delta_N=abs(1.0-ktest*inv(ktrain)*ktest.T)
                            #print "delta_N:",delta_N
                            beta_N=np.sqrt((2*np.log((np.pi)**2*(N**2)))/(6*0.05))
                            #print "beta_N:",beta_N
                            U_N=u_N+beta_N*delta_N
                            #print "U_N:",U_N
                            #print "f_t:",f_t
                            if U_N>=f_t:
                                value=get_energy_ads(center[i][0],center[i][1],center[i][2])
                                function_evalution.append(value)
                                kl=np.argmax(function_evalution)
                                final_result.append(function_evalution[kl])
                                print "centers:",center[i][0],center[i][1],center[i][2]
                                print "function value:",value
                                f_eva=f_eva+1
                                g_value=value
                                #node_function_value.append(g_value)
                                if i==0:
                                    current_node.append(node.Addnode(g_value,center[0][0],center[0][1],center[0][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,e2))
                                if i==1:
                                    current_node.append(node.Addnode(g_value,center[1][0],center[1][1],center[1][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e2,e3))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                                if i==2:
                                    current_node.append(node.Addnode(g_value,center[2][0],center[2][1],center[2][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e3,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                                t=t+1
                                D_x_train=np.c_[D_x_train,center[i]]
                                D_y_train=np.c_[D_y_train,value]
                            else:
                                L_N=u_N-beta_N*delta_N
                                #print "L_N:",L_N
                                g_value=L_N
                                node_function_value.append(g_value)
                                if i==0:
                                    current_node.append(node.Addnode(g_value,center[0][0],center[0][1],center[0][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,e2))
                                if i==1:
                                    current_node.append(node.Addnode(g_value,center[1][0],center[1][1],center[1][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e2,e3))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                                if i==2:
                                    current_node.append(node.Addnode(g_value,center[2][0],center[2][1],center[2][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,node.x2rightmax,e3,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                            if g_value>f_t:
                                f_t=g_value


                    if width_index==1:
                        ele=split8(node.x2leftmin,node.x2rightmax)
                        e2=node.x2leftmin+ele
                        e3=e2+ele

                        centerx11,centerx12,centerx13=split3(node.x1leftmin,node.x1rightmax,node.x2leftmin,e2,node.x3leftmin,node.x3rightmax)
                        centerx21,centerx22,centerx23=split3(node.x1leftmin,node.x1rightmax,e2,e3,node.x3leftmin,node.x3rightmax)
                        centerx31,centerx32,centerx33=split3(node.x1leftmin,node.x1rightmax,e3,node.x2rightmax,node.x3leftmin,node.x3rightmax)


                        leftnode=[centerx11,centerx12,centerx13]
                        midnode=[centerx21,centerx22,centerx23]
                        rightnode=[centerx31,centerx32,centerx33]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)
                        #print "i=1:",center[0][0],center[0][1],center[0][2]

                        #print "i=1:",center[1][0],center[1][1],center[1][2]

                        #print "i=1:",center[2][0],center[2][1],center[2][2]


                        node_function_value=[]

                        for i in range(3):


                            N=N+1
                            D_x_test=np.array([center[i]])
                            ktrain,ktest=kernel_calculation(t,np.matrix(D_x_train),np.matrix(D_x_test))
                            u_N=ktest*inv(ktrain)*np.matrix(D_y_train).T
                            #print "u_N:",u_N
                            delta_N=abs(1.0-ktest*inv(ktrain)*ktest.T)
                            #print "delta_N:",delta_N
                            beta_N=np.sqrt((2*np.log((np.pi)**2*(N**2)))/(6*0.05))
                            #print "beta_N:",beta_N
                            U_N=u_N+beta_N*delta_N
                            #print "U_N:",U_N
                            #print "f_t:",f_t
                            if U_N>=f_t:
                                value=get_energy_ads(center[i][0],center[i][1],center[i][2])
                                function_evalution.append(value)
                                kl=np.argmax(function_evalution)
                                final_result.append(function_evalution[kl])
                                print "centers:",center[i][0],center[i][1],center[i][2]
                                print "function value:",value
                                f_eva=f_eva+1
                                g_value=value
                                #node_function_value.append(g_value)
                                if i==0:
                                    current_node.append(node.Addnode(g_value,center[0][0],center[0][1],center[0][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,e2,node.x3leftmin,node.x3rightmax))
                                if i==1:
                                    current_node.append(node.Addnode(g_value,center[1][0],center[1][1],center[1][2],node.depth+1,node.x1leftmin,node.x1rightmax,e2,e3,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                                if i==2:
                                    current_node.append(node.Addnode(g_value,center[2][0],center[2][1],center[2][2],node.depth+1,node.x1leftmin,node.x1rightmax,e3,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                                t=t+1
                                D_x_train=np.c_[D_x_train,center[i]]
                                D_y_train=np.c_[D_y_train,value]
                            else:
                                L_N=u_N-beta_N*delta_N
                                #print "L_N:",L_N
                                g_value=L_N
                                node_function_value.append(g_value)
                                if i==0:
                                    current_node.append(node.Addnode(g_value,center[0][0],center[0][1],center[0][2],node.depth+1,node.x1leftmin,node.x1rightmax,node.x2leftmin,e2,node.x3leftmin,node.x3rightmax))
                                if i==1:
                                    current_node.append(node.Addnode(g_value,center[1][0],center[1][1],center[1][2],node.depth+1,node.x1leftmin,node.x1rightmax,e2,e3,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                                if i==2:
                                    current_node.append(node.Addnode(g_value,center[2][0],center[2][1],center[2][2],node.depth+1,node.x1leftmin,node.x1rightmax,e3,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                            if g_value>f_t:
                                f_t=g_value
                    #print left_center,right_center
                    if width_index==0:
                        ele=split8(node.x1leftmin,node.x1rightmax)
                        e2=node.x1leftmin+ele
                        e3=e2+ele
                        centerx11,centerx12,centerx13=split3(node.x1leftmin,e2,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax)
                        centerx21,centerx22,centerx23=split3(e2,e3,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax)
                        centerx31,centerx32,centerx33=split3(e3,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax)
                        leftnode=[centerx11,centerx12,centerx13]
                        midnode=[centerx21,centerx22,centerx23]
                        rightnode=[centerx31,centerx32,centerx33]
                        center.append(leftnode)
                        center.append(midnode)
                        center.append(rightnode)
                        #print "i=0:",center[0][0],center[0][1],center[0][2]

                        #print "i=0:",center[1][0],center[1][1],center[1][2]

                        #print "i=0:",center[2][0],center[2][1],center[2][2]

                        node_function_value=[]

                        for i in range(3):
                            N=N+1
                            D_x_test=np.array([center[i]])
                            ktrain,ktest=kernel_calculation(t,np.matrix(D_x_train),np.matrix(D_x_test))
                            u_N=ktest*inv(ktrain)*np.matrix(D_y_train).T
                            #print "u_N:",u_N
                            delta_N=abs(1.0-ktest*inv(ktrain)*ktest.T)
                            #print "delta_N:",delta_N
                            beta_N=np.sqrt((2*np.log((np.pi)**2*(N**2)))/(6*0.05))
                            #print "beta_N:",beta_N
                            U_N=u_N+beta_N*delta_N
                            #print "U_N:",U_N
                            #print "f_t:",f_t
                            if U_N>=f_t:
                                value=get_energy_ads(center[i][0],center[i][1],center[i][2])
                                function_evalution.append(value)
                                kl=np.argmax(function_evalution)
                                final_result.append(function_evalution[kl])
                                print "centers:",center[i][0],center[i][1],center[i][2]
                                print "function value:",value
                                f_eva=f_eva+1
                                g_value=value
                                #node_function_value.append(g_value)
                                if i==0:
                                    current_node.append(node.Addnode(g_value,center[0][0],center[0][1],center[0][2],node.depth+1,node.x1leftmin,e2,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                if i==1:
                                    current_node.append(node.Addnode(g_value,center[1][0],center[1][1],center[1][2],node.depth+1,e2,e3,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                                if i==2:
                                    current_node.append(node.Addnode(g_value,center[2][0],center[2][1],center[2][2],node.depth+1,e3,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                                t=t+1
                                D_x_train=np.c_[D_x_train,center[i]]
                                D_y_train=np.c_[D_y_train,value]
                            else:
                                L_N=u_N-beta_N*delta_N
                                #print "L_N:",L_N
                                g_value=L_N
                                node_function_value.append(g_value)
                                if i==0:
                                    current_node.append(node.Addnode(g_value,center[0][0],center[0][1],center[0][2],node.depth+1,node.x1leftmin,e2,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                if i==1:
                                    current_node.append(node.Addnode(g_value,center[1][0],center[1][1],center[1][2],node.depth+1,e2,e3,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                                if i==2:
                                    current_node.append(node.Addnode(g_value,center[2][0],center[2][1],center[2][2],node.depth+1,e3,node.x1rightmax,node.x2leftmin,node.x2rightmax,node.x3leftmin,node.x3rightmax))
                                    #current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                            if g_value>f_t:
                                f_t=g_value





                    n=n+1
                    v_max=g_value

            #h=h+1
            if f_eva>=300:
                break
        if f_eva>=300:
            break

    return final_result



    #for i in range(len(current_node)):
        #final.append(function(current_node[i].center)
    #findex=np.argmax(final)
    #print function_evalution
    #return current_node[findex].center



d=SOO()
print d
