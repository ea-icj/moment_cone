#Attention, ne fonctionne que pour d_i<=9 et s<=9

lim_st=input('Enter time limit per Grobner computation, (ex "3" for 3 seconds): ')
lim=float(lim_st)

restr_st=input('Probabilist? Y/N ')
if restr_st=="Y":
   restr=True
else:
   restr=False
   print("Formal computations chosen")


#reference list of tw
#in 19/11/24, List_Tau_W is the last list in kron.py still containing all the informations (tau,w)

List_tw=List_Tau_W

def dict_Vtau_neg_ou_nul(tau,Poids,d):
    Res={}
    for chi in Poids:
        p=ScalarProd(chi,tau,d)
        if p<=0:
           if p in Res.keys():
              Res[p].append(chi)
           else :
              Res[p]=[chi]
    return(Res)

def Dict_wUw(tw,d):
    Inv_w=[]
    for k,w in enumerate(tw[3:3+len(d)]): # Inversions of (w1,...,ws)
        Inv_w+=[[k+1,inv[0]-1,inv[1]-1] for inv in w[2]]
    # now, as a dictionary
    dic_Inv_w={}
    for beta in Inv_w:
        p=tw[0][1+sum(d[:beta[0]-1])+beta[1]]-tw[0][1+sum(d[:beta[0]-1])+beta[2]]
        if p not in dic_Inv_w.keys() :
            dic_Inv_w[p]=[beta]
        else :
            dic_Inv_w[p].append(beta)
    return dic_Inv_w

def var_index(index):
    st=""
    for i in index:
       st+=str(i)
    return(st)

def str_of_var(l,seed):
      l_gen=[]
      for xhi in l:
        st=seed+var_index(xhi)
        l_gen.append(st)
      return l_gen

def flat_dict(dic):
   Res=[]
   for x in dic.keys():
      Res+=dic[x]
   return Res

def polyring(Poids,tw,d,restr=False):
   dict_Vneg=dict_Vtau_neg_ou_nul(tw[0],Poids,d)
   fdicv=flat_dict(dict_Vneg)
   gens_V=str_of_var(fdicv,"v")
   dict_Inv_w=Dict_wUw(tw,d)
   fdicu=flat_dict(dict_Inv_w)
   gens_U=str_of_var(fdicu,"u")
   if restr:
      R1=QQ
   else:
      R1=PolynomialRing(QQ,gens_V)
   F1=FractionField(R1)
   R=PolynomialRing(F1,gens_U, order='deglex')
   return R
   
def list_gi(dic_gi):
   dk=list(dic_gi.keys())
   if len(dk)==1:
      return [[x] for x in dic_gi[dk[0]]]
   else:
      dicp=dic_gi.copy()
      del dicp[dk[0]]
      giplus=list_gi(dicp)
      Res=[]
      for g1 in dic_gi[dk[0]]:
         for gi in giplus:
            Res.append([g1]+gi)
      return(Res)



def map_pi(Poids,tw,d,restr=False): 
   """ computes the map pi formally in terms of variables v[i_1,i_2,i_3] in V^{tau<=0} and u[k,i,j] (U_k=1+ sum u_[k,i,j])"""
   R=polyring(Poids,tw,d,restr)
   F1=R.base_ring()
   v=vector(R,prod(d))
   dict_Inv_w=Dict_wUw(tw,d)
   fdicu=flat_dict(dict_Inv_w)
   gi_pos={}
   for s in range(len(d)):
       gi_pos[s+1]=[[s+1,i,i] for i in range(d[s])]
   for u in fdicu:
       if u[0] in gi_pos.keys():
          gi_pos[u[0]].append(u)
       else :
          print("error")
          gi_pos[u[0]]=[u]
   lgi=list_gi(gi_pos)
   for p in flat_dict(dict_Vtau_neg_ou_nul(tw[0],Poids,d)):
      if restr:
          randp=randint(1,1000)
          cp=randp
      else:
          randp=F1.gens_dict()["v"+var_index(p)]
          cp=randp
      for gis in lgi:
         c=cp
         pp=p.copy()
         for k,gi in enumerate(gis):
            if gi[2]+1==p[k]:
               pp[k]=gi[1]+1
               if (gi[1]<gi[2]):  #If, in [k,i,j], we have i=j then the coefficient on the diagonal is 1
                  c*=R.gens_dict()["u"+var_index(gi)]
            else:
              c=0
         v[indexinPoids(pp,d)]+=c
      v[indexinPoids(p,d)]-=randp
   return(v)

def eq_fibre_pi(Poids,tw,d,restr=False):
   eq=[]
   vv=map_pi(Poids,tw,d,restr)
   dict_Vtau_pos=tw[2].copy()
   if 0 in dict_Vtau_pos.keys():
      del dict_Vtau_pos[0]
   for p in flat_dict(dict_Vtau_pos):
      eq.append(vv[indexinPoids(p,d)])
   return eq

def Groebner_basis_brut(tw,restr=False):
   eq=eq_fibre_pi(Poids0,tw,d0,restr)
   if len(eq)>0:
      R=eq[0].parent()
      II=R.ideal(eq) 
      return(II.groebner_basis())
   else:
      return("0 equations")

def Test_deg_Groebner(tw,restr=False):
   GB=Groebner_basis_brut(tw,restr)
   if isinstance(GB, str):
      return(GB)
   print(GB)
   dic={-1:len(GB.ring().gens())}
   for eq in GB:
      dg=eq.degree()
      if dg in dic.keys():
         dic[dg]+=1
      else:
         dic[dg]=1
   dgs=dic.keys()
   t=(len(dgs)==2) and (1 in dgs) and (dic[1]==dic[-1])
   if not(t):
      print(dic,GB)
   return(t)

def Test_deg_Groebner_restr(tw):
   return Test_deg_Groebner(tw,True)

def long_calculation(Liste, fonction, lim, restr=False):
    Res=[]
    for i,l in enumerate(Liste):
       alarm(lim)
       succ=True
       try:
           print('starting calculation..', i)
           resl=fonction(l,restr)
       except:
           print('did not complete! in ', lim, ' seconds')
           succ=False
       # if the computation finished early, though, the alarm is still ticking!
       # so let's turn it off..
       cancel_alarm()
       if succ:
          print(i, "success", resl)
          Res.append([i,l,resl])
       else: 
          print(i, "fail")
    return(Res)


################ Calculs explicites

t1=time.time()
quick_res=long_calculation(List_tw,Test_deg_Groebner,lim,restr)
t2=time.time()
print("truncated Groebner in ", t2-t1, " seconds")

label_bir_quick_res=[]
for m in quick_res:
  if m[2]:
    label_bir_quick_res.append(m[0])

label_bir_quick_res

print("over ", len(List_tw), " initial equalities")
print(len(label_bir_quick_res), " are indeed necessary (birationnality condition)")
print(len(quick_res)-len(label_bir_quick_res), " are non-necessary (non-birat)")
print("For ", len(List_tw)-len(quick_res), " of them, the Grobner computation did not ended within ", lim, " seconds")


"""
ex0=List_Tau_W[7]
vv=map_pi(Poids0,ex0,d0)
eq=eq_fibre_pi(Poids0,ex0,d0)
Test_deg_Groebner(ex0)
R=vv.base_ring()
II=R.ideal([u for u in vv]) 
II.groebner_basis()  


R.inject_variables() 
R.gens_dict()['v113']

import signal

def handler(signum, frame):
     print("Forever is over!")
     raise Exception("end of time")

signal.signal(signal.SIGALRM, handler)


for i,tw in enumerate(List_Tau_W):
    signal.alarm(10)
    try:
       print(i,Test_deg_Groebner(tw))
    except: 
       print("timeout")

 
for i,tw in enumerate(List_Tau_W[29:):
    print(i+29,Test_deg_Groebner(tw))

#Comprendre le 29 (simple à conceptualiser, compliqué pour la base de Grobner)
"""




"""
def polyring(tw,d):
   gens_V=[]
   for xhi in tw[2][0]:
      st="v"
      for i in xhi:
         st+=str(i)
      gens_V.append(st)
   dict_Inv_w=Dict_wUw(tw,d)
   R=PolynomialRing(QQ,gens_V)
   return R
"""
#R = PolynomialRing(QQ, ['u%s'%p for p in primes(100)]+['v%s' for ])




"""
def add_time(cur_times, message)
   tps=time.time()
   cur_times.append(tps)
   stream(message)
   print(message,tps)
"""

"""
def long_calculation(Liste, fonction, lim):
    Res=[]
    for i,l in enumerate(Liste):
       alarm(lim)
       succ=True
       try:
           print('starting calculation..', i)
           resl=fonction(l)
       except AlarmInterrupt:
           print('did not complete!')
           succ=False
       # if the computation finished early, though, the alarm is still ticking!
       # so let's turn it off..
       cancel_alarm()
       if succ:
          print(i, "success", resl)
          Res.append([i,l,resl])
       else: 
          print(i, "fail")
    return(Res)

long_calculation(List_Tau_W,Test_deg_Groebner,3)

def long_calculation2(Liste, fonction):
    Res=[]
    for i,l in enumerate(Liste):
       succ=True
       try:
           print('starting calculation..', i)
           resl=fonction(l)
       except: #KeyboardInterrupt:
           print('did not complete!')
           succ=False
       if succ:
          print(i, "success", resl)
          Res.append([i,l,resl])
       else: 
          print(i, "fail")
    return(Res)





"""





