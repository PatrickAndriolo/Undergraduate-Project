from ncpol2sdpa import generate_operators, SdpRelaxation
import numpy as np
import matplotlib.pyplot as plt


def plot(n_vertices, edges_A, edges_B, n, initial_weight, final_weight):
    
    arr1 = np.array(initial_weight);
    arr2 = np.array(final_weight);

    file = open('UncoloredLovasz.txt','w')
    
    epsilon = np.linspace(0,0.3,n)
    
    i = 0;
    
    while (i < n):
   
        s = np.add((1-epsilon[i])*arr1, (epsilon[i])*arr2);
        
        colored_theta = colored_lovasz(n_vertices, edges_A, edges_B, s);
        
        uncolored_theta = uncolored_lovasz(n_vertices, edges_A, edges_B, s);
        
        print('\n' + '\n' + 'epsilon = ' + str(epsilon[i]) + '\n' +'colored_theta = ' + str(colored_theta) + '\n' +'uncolored_theta = ' + str(uncolored_theta) + '\n' + '\n');
        
        file.write(str(epsilon[i]) + "  " + str(colored_theta) + "  " + str(uncolored_theta) + '\n');
        
        i = i + 1
        
    file.close();
        
    data = np.loadtxt('UncoloredLovasz.txt');
    
    x = data[:, 0];
    y_colored = data[:, 1];
    y_uncolored = data[:, 2]
    plt.plot(x, y_colored,'x', color = 'red', label=r'$\theta_{CHSH}$');
    plt.plot(x, y_uncolored,'.', color = 'black', label=r'$\vartheta_{CHSH}$');
    plt.gca().set_xlabel(r'$\epsilon$')
    plt.gca().set_ylabel(r'$\vartheta(G ,\omega^\epsilon) \ , \ \theta(\mathcal{G} ,\omega^\epsilon)$')
    plt.subplots_adjust(top=1.5)
    plt.legend()
    plt.figure()

def colored_lovasz(n_vertices, edges_A, edges_B, b): 

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

def uncolored_lovasz(n_vertices, edges_A, edges_B, b): 
    
    edges_A = edges_A.union(edges_B)
    edges_B = {}

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


#Example: Evaluate and plot the Lovasz numbers (''colored'' and ''uncolored'') for CHSH on the convex combination of initial and final weight.
#n_vertices = 8
#edges_A = {(0,1),(2,3),(4,5),(6,7),(1,5),(0,4),(2,6),(3,7)}
#edges_B = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6),(3,7)}
#initial_weight = [1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8]
#final_weight = [0,0,0,1/5,1/5,1/5,1/5,1/5]
#n = 15
#plot(n_vertices, edges_A, edges_B, n, initial_weight, final_weight)
    

        
        
    
    
