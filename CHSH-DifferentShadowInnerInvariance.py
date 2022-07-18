from ncpol2sdpa import generate_operators, SdpRelaxation
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig


def plot(n, initial_weight, final_weight):
    
    n_vertices = 8
    
    edges_A_CHSH = {(0,1),(2,3),(4,5),(6,7),(1,5),(0,4),(2,6),(3,7)};
    edges_B_CHSH = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6),(3,7)};

    edges_A1 = {(0,1),(2,3),(4,5),(6,7),(0,4),(2,6),(3,7)};
    edges_B1 = {(1,2),(3,4),(5,6),(0,7),(0,4),(2,6),(3,7)};

    edges_A2 = {(0,1),(2,3),(4,5),(6,7),(2,6),(3,7)};
    edges_B2 = {(1,2),(3,4),(5,6),(0,7),(2,6),(3,7)};

    edges_A3 = {(0,1),(2,3),(4,5),(6,7),(1,5),(3,7)};
    edges_B3 = {(1,2),(3,4),(5,6),(0,7),(1,5),(3,7)};

    edges_A4 = {(0,1),(2,3),(4,5),(6,7),(3,7)};
    edges_B4 = {(1,2),(3,4),(5,6),(0,7),(3,7)};

    edges_A5 = {(0,1),(2,3),(4,5),(6,7)};
    edges_B5 = {(1,2),(3,4),(5,6),(0,7)};
    
    arr1 = np.array(initial_weight);
    arr2 = np.array(final_weight);

    file = open('InnerInvarianceDifferentShadow.txt','w');
    
    epsilon = np.linspace(0,1,n);
    
    i = 0;
    t = 0;
    
    while (i < n):
   
        s = np.add((1-epsilon[i])*arr1, (epsilon[i])*arr2);
        
        theta_CHSH = lovasz(n_vertices, edges_A_CHSH, edges_B_CHSH, s);
        
        t = t + 1
        print('otimização ' + str(t) + '\n')
        
        theta1 = lovasz(n_vertices, edges_A1, edges_B1, s);
        
        t = t + 1
        print('otimização ' + str(t) + '\n')
        
        theta2 = lovasz(n_vertices, edges_A2, edges_B2, s);
        
        t = t + 1
        print('otimização ' + str(t) + '\n')
        
        theta3 = lovasz(n_vertices, edges_A3, edges_B3, s);
        
        t = t + 1
        print('otimização ' + str(t) + '\n')
        
        theta4 = lovasz(n_vertices, edges_A4, edges_B4, s);
        
        t = t + 1
        print('otimização ' + str(t) + '\n')
        
        theta5 = lovasz(n_vertices, edges_A5, edges_B5, s);
        
        t = t + 1
        print('otimização ' + str(t) + '\n')
        
        print('\n' + '\n' + 'epsilon = ' + str(epsilon[i]) + '\n' + 'CHSH = ' + str(theta_CHSH) + '\n' + '34,34 = ' + str(theta1) + '\n' + '114,33 = ' + str(theta2) + '\n' + '33,33 = ' + str(theta3) + '\n' + '113,113 = ' + str(theta4) + '\n' + '1111,1111 = ' + str(theta5) + '\n' + '\n');
        
        file.write(str(epsilon[i]) + "  " + str(theta_CHSH) + "   " + str(theta1) + "   "  + str(theta2) + "   " + str(theta3) + "   " + str(theta4) + "   " + str(theta5) + '\n');
        
        i = i + 1
        
    file.close();
        
    data = np.loadtxt('InnerInvarianceDifferentShadow.txt');

    analytical_epsilon = np.linspace(0,1,100)
    
    alphaCHSH = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A_CHSH, edges_B_CHSH, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A_CHSH, edges_B_CHSH, final_weight)
    alpha1 = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A1, edges_B1, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A1, edges_B1, final_weight)
    alpha2 = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A2, edges_B2, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A2, edges_B2, final_weight)
    alpha3 = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A3, edges_B3, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A3, edges_B3, final_weight)
    alpha4 = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A4, edges_B4, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A4, edges_B4, final_weight)
    alpha5 = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A5, edges_B5, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A5, edges_B5, final_weight)
    
    x = data[:, 0];
    y = data[:, 1];
    #y1 = data[:, 2];
    #y2 = data[:, 3];
    #y3 = data[:, 4];
    #y4 = data[:, 5];
    #y5 = data[:, 6];
    plt.plot(analytical_epsilon,alphaCHSH);
    #plt.plot(analytical_epsilon,alpha1);
    #plt.plot(analytical_epsilon,alpha2);
    #plt.plot(analytical_epsilon,alpha3);
    #plt.plot(analytical_epsilon,alpha4);
    #plt.plot(analytical_epsilon,alpha5);
    #plt.plot(x, y,'.', label=r'$\mathcal{G}_{CHSH}$');
    #plt.plot(x, y1,'x', label=r'$\mathcal{G}_{34,34}$');
    #plt.plot(x, y2,'+', label=r'$\mathcal{G}_{114,33}$');
    #plt.plot(x, y3,'x', label=r'$\mathcal{G}_{33,33}$');
    #plt.plot(x, y4,'+', label=r'$\mathcal{G}_{113,113}$');
    #plt.plot(x, y5,'.', label=r'$\mathcal{G}_{1111,1111}$');
    plt.gca().set_xlabel(r'$\epsilon$')
    plt.gca().set_ylabel(r'$\theta(\mathcal{G} ,\omega^\epsilon)$')
    plt.subplots_adjust(top=1.5)
    plt.legend()
    plt.figure()
    
def lovasz(n_vertices, edges_A, edges_B, b): 

	level = 1
    
	"Adjacency matrices"
	adj_matrix_A = np.zeros((n_vertices,n_vertices))
	adj_matrix_B = np.zeros((n_vertices,n_vertices))
	for i in range(n_vertices):
		for j in range(n_vertices):
			if (i,j) in edges_A:
				adj_matrix_A[i,j] = 1
			if (i,j) in edges_B:
				adj_matrix_B[i,j] = 1

	"Set operators"
	A = generate_operators('A', n_vertices, hermitian=True)
	B = generate_operators('B', n_vertices, hermitian=True)

	"Set objective"
	obj = -sum([b[i]*A[i]*B[i] for i in range(n_vertices)]) #sum of weighted components of the behaviour

	"Substitutions"
	subs = {A[i]**2:A[i] for i in range(n_vertices)}  #conditions of projectors
	subs.update({B[i]**2:B[i] for i in range(n_vertices)})  #conditions of projectors
	(subs.update({B[j]*A[i]:A[i]*B[j] for i in range(n_vertices) for j in
		    range(n_vertices)}))  #symmetry
	(subs.update({A[i]*A[j]:0 for i in range(n_vertices) for j in
		    range(n_vertices) if adj_matrix_A[i,j] == 1}))  #orthogonality relation
	(subs.update({B[i]*B[j]:0 for i in range(n_vertices) for j in
		    range(n_vertices) if adj_matrix_B[i,j] == 1}))  #orthogonality relation

	"Extra monomials"
	extra = ([A[i]*B[j] for i in range(n_vertices) for j in range(n_vertices)])

	"Set problem"
	sdpRelaxation = SdpRelaxation(A+B, verbose=1);
	sdpRelaxation.get_relaxation(level, objective=obj,
		                          substitutions=subs, extramonomials=extra);

	"Solve"
	sdpRelaxation.solve(solver = 'mosek')

	"Final"
	return(-sdpRelaxation.primal)


def pondered_independence_number(n_vertices, edges_A, edges_B, weight):
    
    g = ig.Graph()
    
    g.add_vertices(n_vertices)
    
    edges = list(edges_A) + list(edges_B)
    
    g.add_edges(edges)
    
    arrays_maximais = np.array(g.maximal_independent_vertex_sets(), dtype = object)
    
    i = 0
    j = 0
    
    alpha = np.zeros(shape=(len(arrays_maximais),1))
    
    while(i < len(arrays_maximais)):
        
        j = 0
        
        while(j < len(arrays_maximais[i])):
            
            alpha[i] = alpha[i] + weight[arrays_maximais[i][j]]
            
            j = j + 1
            
        i = i + 1

    return(np.max(alpha))

#Example
#initial_weight = [1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8]
#final_weight = [0,0,0,1/5,1/5,1/5,1/5,1/5]
#n = 15
#plot(n, initial_weight, final_weight)
