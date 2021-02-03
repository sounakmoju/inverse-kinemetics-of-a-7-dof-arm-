import sympy as sp
from sympy.physics.vector import init_vprinting
init_vprinting(use_latex='mathjax', pretty_print=False)
from IPython.display import Image
from sympy.physics.mechanics import dynamicsymbols
import math
from math import pi
import numpy as np
import scipy
import scipy.linalg
#from mat_mul import mat_mul1
import random
#from sympy.physics.mechanics import dynamicsymbols
from sympy import *
import intersection from intersection1
import selectphi from selecting_phi
#jointlimits
#theta_1_u=
#theta_1_l=
#theta_2_u=
#theta_2_l=
#theta_3_u=
#theta_3_l=
#theta_4_u=
#theta_4_l=
#theta_5_u=
#theta_5_l=
#theta_6_u=
#theta_6_l=
alpha,beta,gamma,p_x,p_y,p_z,l_bs,l_se,l_we,l_wt,theta_1,theta_2,theta_3,theta_4,theta_5,theta_6,u_x,u_y,u_z,theta_11,theta_22,theta_33,theta_44,alpha_1,beta_2,phi=dynamicsymbols('alpha,beta,gamma,p_x,p_y,p_z,l_bs,l_se,l_we,l_wt,theta_1,theta_2,theta_3,theta_4,theta_5,theta_6,u_x,u_y,u_z,theta_11,theta_22,theta_33,theta_44,alpha_1,beta_2,phi')
rot_end=sp.Matrix([[sp.cos(beta), -sp.sin(beta)*sp.cos(gamma), sp.sin(beta)*sp.sin(gamma)],
                 [sp.sin(beta)*sp.cos(alpha), (sp.cos(alpha)*sp.cos(beta)*sp.cos(gamma)-sp.sin(alpha)*sp.sin(gamma)), (-sp.cos(alpha)*sp.cos(beta)*sp.sin(gamma)-sp.sin(alpha)*sp.cos(gamma))],
                 [sp.sin(alpha)*sp.sin(beta), (sp.cos(alpha)*sp.sin(gamma)+ sp.cos(beta)*sp.cos(gamma)*sin(alpha)), (sp.cos(alpha)*sp.cos(gamma)-sp.cos(beta)*sp.sin(alpha)*sin(gamma))]])
#Euler angle XZX for orientation#
#print(rot_end[0,2])
pos_end=sp.Matrix([[p_x],[p_y],[p_z]])
#print(pos_end.shape)
w=np.zeros((3,1))
l_wt=sp.Matrix([[0],[0],[l_wt]])

w=pos_end-(rot_end*(l_wt))
w_x=w[0,0]
w_y=w[1,0]
w_z=w[2,0]
L_bs=sp.Matrix([[0],[0],[l_bs]])
#wrist position

l_sw=w-L_bs
#print(l_sw[1,0])
#shoulder to wrist position
k=((l_sw[0,0])*(l_sw[0,0])+(l_sw[1,0])*(l_sw[1,0])+(l_sw[2,0])*(l_sw[2,0]))
mag=math.sqrt((l_sw[0,0])*(l_sw[0,0])+(l_sw[1,0])*(l_sw[1,0])+(l_sw[2,0])*(l_sw[2,0]))
u_norm=l_sw/mag
theta_41=math.acos(l_se**2+l_we**2-k)/(2*l_se*l_we)
theta_4=(pi-theta_41)
u_x=u_norm[0,0]
u_y=u_norm[1,0]
u_z=u_norm[2,0]
if theta_4_l<theta_4<theta_4_u:
    print('no solution')
else:
    skw=sp.Matrix([[0,-u_z,u_y],
                 [u_z,0,-u_x],
                 [-u_y,u_x,0]])
    theta_11=math.atan(w_y,w_x)
    alpha_1=math.asin((w_z-l_bs)/(l_sw))
    beta_2=math.acos((l_se**2+l_we**2 -l_we**2)/(2*l_se*l_we))
    theta_22=(pi/2)-alpha_1-beta_2
    R_01=sp.Matrix([[sp.cos(theta_11),0,-sp.sin(theta_11)],[sp.sin(theta_11),0,sp.cos(theta_11)],[0,-1,0]])
    R_12=sp.Matrix([[sp.cos(theta_22),0,sp.sin(theta_22)],[sp.sin(theta_22),0,-sp.cos(theta_22)],[0,1,0]])
    R_23=sp.Matrix([[1,0,0],[0,0,1],[0,-1,0]])
    R_03=R_01*R_02*R_23
    I=np.identity(3)
    skw2=sp.Matrix([[(u_z**2+u_y**2),(-u_x*u_y),(-u_x*u_z)],
                   [(-u_x*u_y),(u_x**2+u_z**2),(-u_z*u_y)],
                 [(-u_x*u_z),(-u_z*u_y),(u_y**2+u_x**2)]])
    X_S=skw*R_03
    Y_S=-(skw2)*R_03
    Z_S=(I+skw2)*R_03

    theta_1=sp.atan(-(sp.sin(phi)*X_S[1,1]+sp.cos(phi)*Y_S[1,1]+Z_S[1,1]),-(sp.sin(phi)*X_S[0,1]+sp.cos(phi)*Y_S[0,1]+Z_S[0,1]))
    theta_2=sp.acos(-(sp.sin(phi)*X_S[2,1]+sp.cos(phi)*Y_S[2,1]+Z_S[2,1])) 
    d={}
    #all calculations for joints 1,3,5,7
    x_1=Y_S[0,1]*Z_S[1,1]-Y_S[1,1]*Z_S[0,1]
    y_1=X_S[1,1]*Z_S[0,1]-X_S[0,1]*Z_S[1,1]
    z_1=X_S[1,1]*Y_S[0,1]-X_S[0,1]*Y_S[1,1]
    x_3=-Y_S[2,0]*Z_S[2,2]+Y_S[2,2]*Z_S[2,0]
    y_3=-X_S[2,2]*Z_S[2,0]+X_S[2,0]*Z_S[2,2]
    z_3=-X_S[2,2]*Y_S[2,0]+X_S[2,0]*Y_S[2,2]
    x_5=Y_w[0,2]*Z_w[1,2]-Y_w[1,2]*Z_w[0,2]
    y_5=X_w[1,2]*Z_w[0,2]-X_w[0,2]*Z_w[1,2]
    z_5=X_w[1,2]*Y_w[0,2]-X_w[0,2]*Y_w[1,2]
    x_7=-Y_w[2,0]*Z_w[2,1]+Y_w[2,1]*Z_w[2,0]
    y_7=-X_w[2,1]*Z_w[2,0]+X_w[2,0]*Z_w[2,1]
    z_7=-X_w[2,1]*Y_w[2,0]+X_w[2,0]*Y_w[2,1]
    v_1=x_1**2+y_1**2-z_1**2
    v_3=x_3**2+y_3**2-z_3**2
    v_5=x_4**2+y_4**2-z_4**2
    v_7=x_7**2+y_7**2-z_7**2
    a_11_u=X_S[0,1]*(math.tan(theta_1_u))-X_S[1,1]
    b_11_u=Y_S[0,1]*(math.tan(theta_1_u))-Y_S[1,1]
    c_11_u=Z_S[0,1]*(math.tan(theta_1_u))-Z_S[1,1]
    a_11_l=X_S[0,1]*(math.tan(theta_1_l))-X_S[1,1]
    b_11_l=Y_S[0,1]*(math.tan(theta_1_l))-Y_S[1,1]
    c_11_l=Z_S[0,1]*(math.tan(theta_1_l))-Z_S[1,1]
    a_33_u=-X_S[2,0]*(math.tan(theta_3_u))-X_S[2,2]
    b_33_u=-Y_S[2,0]*(math.tan(theta_3_u))-Y_S[2,2]
    c_33_u=-Z_S[2,0]*(math.tan(theta_3_u))-Z_S[2,2]
    a_33_l=-X_S[2,0]*(math.tan(theta_3_l))-X_S[2,2]
    b_33_l=-Y_S[2,0]*(math.tan(theta_3_l))-Y_S[2,2]
    c_33_l=-Z_S[2,0]*(math.tan(theta_3_l))-Z_S[2,2]
    a_55_u=X_w[0,2]*(math.tan(theta_5_u))-X_w[1,2]
    b_55_u=Y_w[0,2]*(math.tan(theta_5_u))-Y_w[1,2]
    c_55_u=Z_w[0,2]*(math.tan(theta_5_u))-Z_w[1,2]
    a_55_l=X_w[0,2]*(math.tan(theta_5_l))-X_w[1,2]
    b_55_l=Y_w[0,2]*(math.tan(theta_5_l))-Y_w[1,2]
    c_55_l=Z_w[0,2]*(math.tan(theta_5_l))-Z_w[1,2]
    a_55_u=X_w[0,2]*(math.tan(theta_5_u))-X_w[1,2]
    b_55_u=Y_w[0,2]*(math.tan(theta_5_u))-Y_w[1,2]
    c_55_u=Z_w[0,2]*(math.tan(theta_5_u))-Z_w[1,2]
    a_55_l=X_w[0,2]*(math.tan(theta_5_l))-X_w[1,2]
    b_55_l=Y_w[0,2]*(math.tan(theta_5_l))-Y_w[1,2]
    c_55_l=Z_w[0,2]*(math.tan(theta_5_l))-Z_w[1,2]
    a_77_u=-X_w[2,0]*(math.tan(theta_7_u))-X_w[2,1]
    b_77_u=-Y_w[2,0]*(math.tan(theta_7_u))-Y_w[2,1]
    c_77_u=-Z_w[2,0]*(math.tan(theta_7_u))-Z_w[2,1]
    a_77_l=-X_w[2,0]*(math.tan(theta_7_l))-X_w[2,1]
    b_77_l=-Y_w[2,0]*(math.tan(theta_7_l))-Y_w[2,1]
    c_77_l=-Z_w[2,0]*(math.tan(theta_7_l))-Z_w[2,1]
    phi_1_u_ma=math.atan(a_11_u,b_11_u)+math.acos(c_11_u,math.sqrt((a_11_u**2,b_11_u**2)))
    phi_1_u_mi=math.atan(a_11_u,b_11_u)-math.acos(c_11_u,math.sqrt((a_11_u**2,b_11_u**2)))
    phi_1_l_ma=math.atan(a_11_l,b_11_l)+math.acos(c_11_l,math.sqrt((a_11_l**2,b_11_l**2)))
    phi_1_l_mi=math.atan(a_11_l,b_11_l)-math.acos(c_11_l,math.sqrt((a_11_l**2,b_11_l**2)))

    phi_3_u_ma=math.atan(a_33_u,b_33_u)+math.acos(c_33_u,math.sqrt((a_33_u**2,b_33_u**2)))
    phi_3_u_mi=math.atan(a_33_u,b_33_u)-math.acos(c_33_u,math.sqrt((a_33_u**2,b_33_u**2)))
    phi_3_l_ma=math.atan(a_33_l,b_33_l)+math.acos(c_33_l,math.sqrt((a_33_l**2,b_33_l**2)))
    phi_3_l_mi=math.atan(a_33_l,b_33_l)-math.acos(c_33_l,math.sqrt((a_33_l**2,b_33_l**2)))

    phi_5_u_ma=math.atan(a_55_u,b_55_u)+math.acos(c_55_u,math.sqrt((a_55_u**2,b_55_u**2)))
    phi_5_u_mi=math.atan(a_55_u,b_55_u)-math.acos(c_55_u,math.sqrt((a_55_u**2,b_55_u**2)))
    phi_5_l_ma=math.atan(a_55_l,b_55_l)+math.acos(c_55_l,math.sqrt((a_55_l**2,b_55_l**2)))
    phi_5_l_mi=math.atan(a_33_l,b_33_l)-math.acos(c_33_l,math.sqrt((a_33_l**2,b_33_l**2)))

    phi_7_u_ma=math.atan(a_77_u,b_77_u)+math.acos(c_77_u,math.sqrt((a_77_u**2,b_77_u**2)))
    phi_7_u_mi=math.atan(a_77_u,b_55_u)-math.acos(c_77_u,math.sqrt((a_77_u**2,b_77_u**2)))
    phi_7_l_ma=math.atan(a_77_l,b_77_l)+math.acos(c_77_l,math.sqrt((a_77_l**2,b_77_l**2)))
    phi_7_l_mi=math.atan(a_77_l,b_77_l)-math.acos(c_33_l,math.sqrt((a_77_l**2,b_77_l**2)))
    #following maximum minimum has to be decided by the second derivative test 

    phi_cr_1_ma=2*math.atan((x_1-v_1)/(y_1-z_1))
    phi_cr_1_mi=2*math.atan((x_1+v_1)/(y_1-z_1))
    phi_cr_3_ma=2*math.atan((x_3-v_3)/(y_3-z_3))
    phi_cr_3_mi=2*math.atan((x_3+v_3)/(y_3-z_3))
    phi_cr_5_ma=2*math.atan((x_5-v_5)/(y_5-z_5))
    phi_cr_5_mi=2*math.atan((x_5+v_5)/(y_5-z_5))
    phi_cr_7_ma=2*math.atan((x_5-v_5)/(y_5-z_5))
    phi_cr_7_mi=2*math.atan((x_5+v_5)/(y_5-z_5))
    phi_singu_1=2*math.atan(x_1/(y_1-z_1))
    phi_singu_3=2*math.atan(x_3/(y_3-z_3))
    phi_singu_5=2*math.atan(x_5/(y_5-z_5))
    phi_singu_7=2*math.atan(x_7/(y_7-z_7))
    d_singu1={1:phi_singu_1,3:phi_singu_2,5:phi_singu_5,7: phi_singu_7}
    
    
    

    d_2={1:[phi_cr_1_ma,phi_cr_1_mi],3:[phi_cr_3_ma,phi_cr_3_mi],5:[phi_cr_5_ma,phi_cr_5_mi],7:[phi_cr_7_ma,phi_cr_7_mi]}
    d_1={1:[phi_1_u_ma,phi_1_u_mi,phi_1_l_ma,phi_1_l_mi],3:[phi_3_u_ma,phi_3_u_mi,phi_3_l_ma,phi_3_l_mi],5:[phi_5_u_ma,phi_5_u_mi,phi_5_l_ma,phi_5_l_mi],7:[phi_7_u_ma,phi_7_u_mi,phi_7_l_ma,phi_7_l_mi]}
    d={1:v_1,3:v_3,5:v_5,7:v_7}
    d_4={1:theta_1,3:theta_3,5:theta_5,7:theta_5]
    d_5={1:[theta_1_u,theta_1_l],3:[theta_2_u,theta_2_l],5:[theta_5_u,theta_5_l],7:[theta_7_u,theta_7_l]
    #d_3=dict.fromkeys([1,2,3,4,5,6], list())
    list_1=[]
    list_1_singular=[]

    for key in d:
         
         if d[key]<0:
         
            list1=[]
            list1.append(d_1[key][3])
            list1.append(d_1[key][0])
            list_1.append(list1)
            list1.clear()
         elif d[key]>0:
            list_2=[]
            phi_max=d_2[key][0]
            phi_min=d_2[key][1]
            theta_max=d_4[key('phi_max')]
            theta_min=d_4[key('phi_max')]
            if theta_max<d_5[key][0] and theta_min<d_5[key][0] and theta_max>d_5[key][1] and theta_min>d_5[key][1]:
                list1.append(-180)
                list1.append(180)
                list_1.append(list1)
                list1.clear()
            elif theta_max>d_5[key][0] and theta_min<d_5[key][0] and theta_min>d_5[key][1]:
                list1.append(-180)
                list1.append(d_1[key][1])
                list_1.append(list1)
                list1.clear()
                list1.append(d_1[key][0])
                list1.append(180)
                list_1.append(list1)
                list1.clear()
            elif theta_max<d_5[key][0] and theta_max>d_5[key][1] and theta_min<d_5[key][1]:
                list1.append(-180)
                list1.append(d_1[key][3])
                list_1.append(list1)
                list_1.clear()
                list1.append(d_1[key][2])
                list1.append(180)
                list_1.append(list1)
                list_1.clear()
            elif theta_max>d_5[key][0] and theta_min<d_5[key][1]:
                list1.append(-180)
                list1.append(d_1[key][1])
                list_1.append(list1)
                list_1.clear()
                list1.append(d_1[key][0])
                list1.append(d_1[key][3])
                list_1.append(list1)
                list1.clear()
                list1.append(d_1[key][2])
                list1.append(180)
                list_1.append(list1)
                list1.clear
            else:
                print("no solution exist")
                
        else:
            list_1_singular.append(d_singu1[key])
        #all calculations for joints 1,3,5,7
        #all calculations for joints 2,6
        #theta_2=sp.acos(-(sp.sin(phi)*X_S[2,1]+sp.cos(phi)*Y_S[2,1]+Z_S[2,1]))
        #theta_6=sp.asin(-(sp.sin(phi)*X_w[2,2]+sp.cos(phi)*Y_w[2,2]+Y_w[2,2]))
    phi_2_cri=math.atan(X_S[2,1],Y_S[2,1])
    phi_6_cri=math.atan(X_w[2,2],Y_w[2,2])
    phi_2_cri2=theta_2_cri+pi
    phi_2_cri21=theta_2_cri-pi
    phi_6_cri2=theta_6_cri+pi
    phi_6_cri21=theta_6_cri-pi
        #after checking the sign and the range 
    phi_2_max=phi_2_cri
    phi_2_min=phi_2_cri21
    phi_6_max=phi_6_cri
    phi_6_min=phi_6_cri21

    a_22_u=X_S[2,1]
    b_22_u=Y_S[2,1]
    c_22_u=(math.cos(theta_2_u))-Z_S[2,1]
    a_22_l=X_S[2,1]
    b_22_l=Y_S[2,1]
    c_22_l=(math.cos(theta_2_l))-Z_S[2,1]
    a_66_u=X_w[2,2]
    b_66_u=Y_w[2,2]
    c_66_u=(math.sin(theta_6_u))-Z_S[2,2]
    a_66_l=X_w[2,2]
    b_66_l=Y_w[2,2]
    c_66_l=(math.sin(theta_6_l))-Z_S[2,2]
    phi_2_u_ma=math.atan(a_22_u,b_22_u)+math.acos(c_22_u,math.sqrt((a_22_u**2,b_22_u**2)))
    phi_2_u_mi=math.atan(a_22_u,b_22_u)-math.acos(c_22_u,math.sqrt((a_22_u**2,b_22_u**2)))
    phi_2_l_ma=math.atan(a_22_l,b_22_l)+math.acos(c_22_l,math.sqrt((a_22_l**2,b_22_l**2)))
    phi_2_l_mi=math.atan(a_11_l,b_11_l)-math.acos(c_11_l,math.sqrt((a_11_l**2,b_11_l**2)))

    phi_6_u_ma=math.atan(a_66_u,b_66_u)+math.acos(c_66_u,math.sqrt((a_66_u**2,b_66_u**2)))
    phi_6_u_mi=math.atan(a_66_u,b_33_u)-math.acos(c_66_u,math.sqrt((a_66_u**2,b_66_u**2)))
    phi_6_l_ma=math.atan(a_66_l,b_66_l)+math.acos(c_66_l,math.sqrt((a_66_l**2,b_66_l**2)))
    phi_6_l_mi=math.atan(a_66_l,b_33_l)-math.acos(c_66_l,math.sqrt((a_66_l**2,b_66_l**2)))
    #2 tan−1xcc − zc + 1
    phi_singu_2_1=math.atan(X_S[2,1],(Y_S[2,1]-Z_S[2,1]+1))
    phi_singu_2_1=math.atan(X_S[2,1],(Y_S[2,1]-Z_S[2,1]-1))
    phi_singu_6_1=math.atan(X_w[2,1],(Y_w[2,1]-Z_w[2,1]+1)
    phi_singu_6_1=math.atan(X_S[2,1],(Y_S[2,1]-Z_S[2,1]-1))
    
    
    

    d_22={2:[ phi_2_max,phi_2_min],6:[phi_6_max,phi_6_min]}
    d_11={2:[phi_2_u_ma,phi_2_u_mi,phi_2_l_ma,phi_2_l_mi],6:[phi_6_u_ma,phi_6_u_mi,phi_6_l_ma,phi_6_l_mi]}
    #d={1:v_1,3:v_3,5:v_5,7:v_7}
    d_44={2:theta_2,6:theta_6]
    d_55={2:[theta_2_u,theta_2_l],6:[theta_6_u,theta_6_l]}
    d_singu1={2:phi_singu_2,6:phi_singu_6}
    
    for key in d_22:
        phi_max=d_22[key][0]
        phi_min=d_22[key][1]
        theta_max=d_44[key('phi_max')]
        theta_min=d_44[key('phi_max')]
        if theta_max<d_55[key][0] and theta_min<d_55[key][0] and theta_max>d_55[key][1] and theta_min>d_55[key][1]:
            list1.append(-180)
            list1.append(180)
            list_1.append(list1)
            list1.clear()
        elif theta_max>d_55[key][0] and theta_min<d_55[key][0] and theta_min>d_55[key][1]:
            list1.append(-180)
            list1.append(d_11[key][1])
            list_1.append(list1)
            list1.clear()
            list1.append(d_11[key][0])
            list1.append(180)
            list_1.append(list1)
            list1.clear()
        elif theta_max<d_55[key][0] and theta_max>d_55[key][1] and theta_min<d_55[key][1]:
            list1.append(-180)
            list1.append(d_11[key][3])
            list_1.append(list1)
            list_1.clear()
            list1.append(d_11[key][2])
            list1.append(180)
            list_1.append(list1)
            list_1.clear()
        elif theta_max>d_55[key][0] and theta_min<d_55[key][1]:
            list1.append(-180)
            list1.append(d_11[key][1])
            list_1.append(list1)
            list_1.clear()
            list1.append(d_11[key][0])
            list1.append(d_11[key][3])
            list_1.append(list1)
            list1.clear()
            list1.append(d_11[key][2])
            list1.append(180)
            list_1.append(list1)
            list1.clear
        else:
            print("no solution exist")
                
        #else:
            #list_1_singular.append(d_singu1[key])
    list_12=intersection(list_1)
    phi_selection=selectphi(list_12,list_1_singular)
        


    
    
    
            
            
            
               

