from math import *
import numpy as np
import sys
#sys.setrecursionlimit(1000)
import matplotlib.pyplot as plt
from numpy import linalg as LA
from numpy.linalg import inv


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
            ktr_new[i,j]=(1+np.sqrt((5*LA.norm(xtr[:,i]-xtr[:,j])**2)/0.15)+
            5*LA.norm(xtr[:,i]-xtr[:,j])**2/(3*0.15))*np.exp(-(np.sqrt(5*LA.norm(xtr[:,i]-xtr[:,j])**2)/0.15))
            #print ktr[i,j]
            #ktr_new[j,i]=ktr_new[i,j]

    for i in range(1):
        for j in range(0,idtr):
            #kte_new[i,j]=exp(-0.5*0.03*(xte[:,i]-xtr[:,j])**2)
            kte_new[i,j]=(1+np.sqrt((5*LA.norm(xte[:,i]-xtr[:,j])**2)/0.15)+
            5*LA.norm(xte[:,i]-xtr[:,j])**2/(3*0.15))*np.exp(-np.sqrt((5*LA.norm(xte[:,i]-xtr[:,j])**2)/0.15))

    #print kte_new.shape

    return ktr_new,kte_new




class Node:
    def __init__(self, center=None, value=None, depth=None,leftmin=None, rightmax=None, parent = None):

        self.parentNode = parent
        self.childNodes = []
        self.center=center
        self.value=value
        self.depth=depth
        self.leftmin=leftmin
        self.rightmax=rightmax

    def Selectnode(self,node,h):
        leaves_of_depth=[]
        index=[]
        value=[]
        for i in range(len(node)):
            if node[i].depth==h:
                leaves_of_depth.append(node[i])
                index.append(i)

        for i in range(len(leaves_of_depth)):
            value.append(leaves_of_depth[i].value)
            #print value
        new_index=np.argmax(value)
        #print new_index
        #hi=node[index[new_index]].center
        #fh=function(hi)
        #print hi
        #print fh
        fh=node[index[new_index]].value
        #print "function value:",fh
        #gl=function(node[index[new_index]].center)
        #print gl


        #function_evalution.append(fh)

        return node[index[new_index]],index[new_index],fh
    def Addnode(self,center_value,value, depth,leftmargin,rightmargin):
        n = Node(center=center_value,value=value, depth=depth, leftmin=leftmargin,rightmax=rightmargin, parent = self)
        self.childNodes.append(n)
        return n



def SOO():
    """set g(0,0)=f(x(0,0))"""
    ini_f=float("-inf")
    g_function=function(0.5)

    """set f^=g(0,0)"""
    f_t=g_function

    """initialize the tree"""
    rootnode = Node(center=0.5, value=g_function,depth=0,leftmin=0, rightmax=1)


    """set t=1,n=1,N=1 and D_t={(x00,g00)}"""
    t=1 ## represent the number of sampled data
    n=1 ## represent the iteration of the algorithm
    N=1 ## represent the number of GP calculation
    D_x_train=np.array([0.5])
    D_y_train=np.array([g_function])

    current_node=[]
    current_node.append(rootnode)
    #print current_node[0].center
    node=rootnode
    leaf=[]
    final=[]
    function_evalution=[]
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
        h=0
        for h in range(0,n):
            #print len(current_node)
            #print h
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
                    mid=midpoint(node.leftmin, node.rightmax)
                    left_center=midpoint(node.leftmin,mid)
                    right_center=midpoint(mid,node.rightmax)
                    #print left_center,right_center
                    center.append(left_center)
                    center.append(right_center)
                    node_function_value=[]

                    for i in range(2):


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
                            value=function(center[i])
                            print "function value:",value
                            f_eva=f_eva+1
                            g_value=value
                            #node_function_value.append(g_value)
                            if i==0:
                                current_node.append(node.Addnode(center[i],g_value,node.depth+1,node.leftmin,mid))
                            if i==1:
                                current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))
                            t=t+1
                            D_x_train=np.c_[D_x_train,center[i]]
                            D_y_train=np.c_[D_y_train,value]
                        else:
                            L_N=u_N-beta_N*delta_N
                            #print "L_N:",L_N
                            g_value=L_N
                            node_function_value.append(g_value)
                            if i==0:
                                current_node.append(node.Addnode(center[i],g_value,node.depth+1,node.leftmin,mid))
                            if i==1:
                                current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                            #if i==0:
                            #    current_node.append(node.Addnode(center[i],g_value,node.depth+1,node.leftmin,mid))
                            #if i==1:
                            #    current_node.append(node.Addnode(center[i],g_value,node.depth+1,mid,node.rightmax))

                        if g_value>f_t:
                            f_t=g_value

                    n=n+1
                    v_max=g_value

            h=h+1
            if f_eva>=200:
                break
        if f_eva>=200:
            break



    #for i in range(len(current_node)):
        #final.append(function(current_node[i].center)
    #findex=np.argmax(final)
    #print function_evalution
    #return current_node[findex].center



d=SOO()
#print d
