from ncpol2sdpa import generate_operators, SdpRelaxation
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig


def plot(n_vertices, edges_A1, edges_B1, edges_A2, edges_B2,  n, initial_weight, final_weight):
    
    arr1 = np.array(initial_weight);
    arr2 = np.array(final_weight);

    file = open('QSet.txt','w')
    
    epsilon = np.linspace(0,1,n)
    
    i = 0;
    
    while (i < n):
   
        s = np.add((1-epsilon[i])*arr1, (epsilon[i])*arr2);
        
        theta1 = lovasz(n_vertices, edges_A1, edges_B1, s);
        theta2 = lovasz(n_vertices, edges_A2, edges_B2, s);
        
        print('epsilon = ' + str(epsilon[i]) + '   ||   ' +'theta1 = ' + str(theta1) + '    //    ' + str(theta2) + '\n');
        
        file.write(str(epsilon[i]) + "  " + str(theta1) + "   " + str(theta2) + '\n');
        
        i = i + 1
        
    file.close();
        
    data = np.loadtxt('QSet.txt');
    
    analytical_epsilon = np.linspace(0,1,100)
    
    alpha1 = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A1, edges_B1, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A1, edges_B1, final_weight)
    alpha2 = (1-analytical_epsilon)*pondered_independence_number(n_vertices, edges_A2, edges_B2, initial_weight) + analytical_epsilon*pondered_independence_number(n_vertices, edges_A2, edges_B2, final_weight)
    
    x = data[:, 0];
    y1 = data[:, 1];
    y2 = data[:, 2];
    plt.plot(analytical_epsilon,alpha1,'orange', label = r'$\alpha(\mathcal{G}_{CHSH}), \alpha(\mathcal{G}_{43,44})$');
    plt.plot(analytical_epsilon,alpha2,'cyan');
    plt.plot(x, y1,'x', color = 'blue', label = r'$\theta(\mathcal{G}_{CHSH})$');
    plt.plot(x, y2,'+', color = 'red', label = r'$\theta(\mathcal{G}_{43,44})$');
    plt.gca().set_xlabel(r'$\epsilon$')
    plt.gca().set_ylabel(r'$\alpha(\mathcal{G} ,\omega^\epsilon), \theta(\mathcal{G} ,\omega^\epsilon)$')
    plt.subplots_adjust(top=1.5)
    plt.legend()
    plt.figure()
    
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


#Example: Evaluates and plot the curves for alpha and theta (the classical and quantum upper bound of a Bell/NC inequality, respectively) for a given convex combination of initial and final weights.
#n_vertices = 8
#edges_A1 = {(0,1),(2,3),(4,5),(6,7),(3,7)};
#edges_B1 = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6),(3,7)};
#edges_A2 = {(0,1),(2,3),(4,5),(6,7),(3,7)};
#edges_B2 = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6)};
#final_weight = [0,0,0,1/5,1/5,1/5,1/5,1/5]
#initial_weight = [0,1/8,1/8,1/8,1/8,1/8,1/8,1/8]
#n = 4
#plot(n_vertices, edges_A1, edges_B1, edges_A2, edges_B2, n, initial_weight, final_weight)
    

        
        
    
    
