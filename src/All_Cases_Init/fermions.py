def hyp_fermions(V):
    resnd=[]
    res=[]
    weights=[chi for chi in V.all_weights]
    #print(weights)
    for S in itertools.combinations(weights,V.G.rank-1):
        if check_hyperplane_dim(S, V.G.rank-1):
            tau=Tau.from_zero_weights(S, V)
            print('tau',tau)
            resnd+=[tau,tau.opposite]
            if tau.is_dominant :
                res.append(tau)
            if tau.opposite.is_dominant :
                res.append(tau.opposite)
    print('Non dom',list(set(resnd)),'\n')   
    return list(set(res))

weights=[chi for chi in V.weights_of_S(Partition([3,2]))]
Gred=LinGroup([2])
Vred=Representation(Gred,'fermion',2)
for S in find_hyperplanes_reg(weights, Vred, 30,1, [1,1]):
    print(S,Tau.from_zero_weights(S, Vred))

find_1PS(V)
# Pb 1 : opposite
# Pb 2 : pas invariant par S5
