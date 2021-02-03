import random
def selectphi(lost,singular):
    phi_1=[]
    phi_1=list((lost)-set(singular))
    
    phi=random.choice(phi_1)

    return phi
##d={1,2,8,9,0,6,45}
##d_1=[2,9,6]
##print(selectphi(d,d_1))
    
    
    
