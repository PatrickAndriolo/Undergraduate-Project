from ncpol2sdpa import generate_operators, SdpRelaxation
import numpy as np
import matplotlib.pyplot as plt


def plot(n, initial_weight, final_weight):
    
    n_vertices = 8
    
    edges_A_CHSH = {(0,1),(2,3),(4,5),(6,7),(1,5),(0,4),(2,6),(3,7)};
    edges_B_CHSH = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6),(3,7)};

    edges_A1 = {(0,1),(2,3),(4,5),(6,7),(1,5),(2,6),(3,7)};
    edges_B1 = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6),(3,7)};
    
    edges_A2 = {(0,1),(2,3),(4,5),(6,7),(1,5),(2,6)};
    edges_B2 = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6),(3,7)};

    edges_A3 = {(0,1),(2,3),(4,5),(6,7),(2,6)};
    edges_B3 = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(2,6),(3,7)};

    edges_A4 = {(0,1),(2,3),(4,5),(6,7),(2,6)};
    edges_B4 = {(1,2),(3,4),(5,6),(0,7),(1,5),(0,4),(3,7)};

    arr1 = np.array(initial_weight);
    arr2 = np.array(final_weight);

    file = open('InnerInvarianceSameShadow.txt','w');
    
    epsilon = np.linspace(0,1,n);
    
    i = 0;
    t = 0;
    
    while (i < n):
   
        s = np.add((1-epsilon[i])*arr1, (epsilon[i])*arr2);
        
        theta_CHSH = lovasz(n_vertices, edges_A_CHSH, edges_B_CHSH, s);
        
        t = t + 1
        print('otimização ' + str(t) + '/' + str(5*n) + '\n')
        
        theta1 = lovasz(n_vertices, edges_A1, edges_B1, s);
        
        t = t + 1
        print('otimização ' + str(t) + '/' + str(5*n) + '\n')
        
        theta2 = lovasz(n_vertices, edges_A2, edges_B2, s);
        
        t = t + 1
        print('otimização ' + str(t) + '/' + str(5*n) + '\n')
        
        theta3 = lovasz(n_vertices, edges_A3, edges_B3, s);
        
        t = t + 1
        print('otimização ' + str(t) + '/' + str(5*n) + '\n')
        
        theta4 = lovasz(n_vertices, edges_A4, edges_B4, s);
        
        t = t + 1
        print('otimização ' + str(t) + '/' + str(5*n) + '\n')
        
        print('\n' + '\n' + 'epsilon = ' + str(epsilon[i]) + '\n' + 'CHSH = ' + str(theta_CHSH) + '\n' + '34,44 = ' + str(theta1) + '\n' + '114,44 = ' + str(theta2) + '\n' + '33,44 = ' + str(theta3) + '\n' + '113,44 = ' + str(theta4) + '\n' + '1111,44 = ' + '\n');
        
        file.write(str(epsilon[i]) + "  " + str(theta_CHSH) + "   " + str(theta1) + "   "  + str(theta2) + "   " + str(theta3) + "   " + str(theta4) + "   " + '\n');
        
        #print('\n' + '\n' + 'epsilon = ' + str(epsilon[i]) + '\n' + '43,44 = ' + str(theta1) + '\n' + '44,33 = ' + str(theta2) + '\n')
        
        #file.write(str(epsilon[i]) + "  " + str(theta1) + "   "  + str(theta2) + '\n')
        
        i = i + 1
        
    file.close();
    
        
    data = np.loadtxt('InnerInvarianceSameShadow.txt');

    x = data[:, 0];
    y = data[:, 1];
    y1 = data[:, 2];
    y2 = data[:, 3];
    y3 = data[:, 4];
    y4 = data[:, 5];
    plt.plot(x, y,'x', label=r'$\mathcal{G}_{CHSH}$');
    plt.plot(x, y1,'x', label=r'$\mathcal{G}_{43,44}$');
    plt.plot(x, y2,'x', label=r'$\mathcal{G}^{2}_{33,44}$');
    plt.plot(x, y3,'x', label=r'$\mathcal{G}_{311,44}$');
    plt.plot(x, y4,'x', label=r'$\mathcal{G}_{311,43}$');
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


#Exemplo:
    
#initial_weight = [1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8]
#final_weight = [0,0,0,1/5,1/5,1/5,1/5,1/5]
#n = 15
#plot(n, initial_weight, final_weight)
