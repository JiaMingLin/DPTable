import math
import random
import logging
import operator
import time
import numpy as np
import cvxpy as cvx
from pathos.multiprocessing import ProcessingPool as Pool


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
"""
Util functions defines
"""
class MarginalOptimization:

    def __init__(self, input_file, clusters_num, _lambda, output_file=None):
        

        self.node_card, self.cliques = self.get_jt(input_file)

        self._lambda = float(_lambda)
        self.output_file = output_file
        self.max_iter = 20

        self.nodes_num = len(self.node_card)
        self.cliques_num = len(self.cliques)
        self.clusters_num = clusters_num
        
        logging.info("Node numbers: %d. Cliques number: %d. Clusters number: %d" % (self.nodes_num, self.cliques_num, self.clusters_num))

    """
    The construtor of difference operator.
    :param m: 
        Number of clusters
    """ 
    def construct_difference(self):
        m = self.clusters_num
        M = np.zeros((m,m*(m-1)/2), dtype=int)
        count = 0
        for i in range(m):
            M[i,count:(count+m-i-1)] = 1
            M[i+1:m, count:(count+m-i-1)] = -np.eye(m-i-1)
            count=count+m-i-1;
        return M

    """
    The matrix representation of junction tree
    :param d: 
          Number of attributes.
    :param cliques:
        cliques.
    """ 
    def jt_rep(self):
        cliques = self.cliques
        n = self.cliques_num
        d = self.nodes_num

        O = np.zeros((d,n),dtype=int)
        for i in range(n):
            O[cliques[i], i] = 1
        
        return O


    def get_jt(self, jt_filepath):
        node_card = []; cliques = []
        with open(jt_filepath, 'r') as fin:
            data = fin.read().splitlines(True)
            node_card_option = 1
            for line in data:
                if len(line) ==0:
                    continue
                if '-' in line:
                    node_card_option = 0
                elif node_card_option == 1:                    
                    node_card.append(int(line.split(' ')[1]))
                else:                                   
                    cliques.append([int(col)-1 for col in line.split(' ') if len(col)>0 and col != '\n'])

        return node_card, cliques

    """
    Product of the attributes' domain size in the each clique
        : param node_card: 
            the attributes' domain size in array.
        : param cliques:
            the cliques and attributes bucket.
    """
    def log_p_func(self):
        node_card = self.node_card
        cliques = self.cliques
        log_p_arr = []
        for i in range(len(cliques)):
            value = sum(np.log(list(node_card[k] for k in cliques[i])))
            log_p_arr.append(value)
        return np.asarray(log_p_arr).T

    """
    Do a marginal optimization
        :param d: 
            Number of attributes
        :param m:
            Number of cluster
        :param _lambda:
            Constant for balance
        :param max_iter:
            Maximun iteration
    """
    def marginal_optimization(self, seed = None):
        logging.debug("Starting to merge marginals")
        # get number of cliques: n
        node_card = self.node_card; cliques = self.cliques
        d = self.nodes_num; n = self.cliques_num; m = self.clusters_num
    
        # get the junction tree matrix representation: O
        O = self.jt_rep()
        
        # get log_p is the array of numbers of sum(log(attribute's domain))
        log_p = self.log_p_func()
    
        # get log_node_card: log(C1), log(C2), ..., log(Cd)
        log_node_card = np.log(node_card)
    
        # get value of sum_log_node_card: log(C1 * C2 *...* Cd)
        sum_log_node_card = sum(log_node_card)
    
        # get the difference operator M on cluster number: m
        M = self.construct_difference()
        # initial a seed Z
        prev_Z = seed
        if prev_Z is None:
            prev_Z = np.random.rand(n,m)        
    
        # run the convex optimization for max_iter times
        logging.debug("Optimization starting...")
        for i in range(self.max_iter):
            logging.debug("The optimization iteration: "+str(i+1))
            # sum of row of prev_Z
            tmp1 = cvx.sum_entries(prev_Z, axis=0).value
        
            # tmp2 = math.log(tmp1)-1+sum_log_node_card
            tmp2 = np.log(tmp1)-1+sum_log_node_card

            # tmp3: difference of pairwise columns = prev_Z * M
            tmp3 = np.dot(prev_Z,M)
            # convex optimization
            Z = cvx.Variable(n,m)
            t = cvx.Variable(1,m)
            r = cvx.Variable()
        
            objective = cvx.Minimize(cvx.log_sum_exp(t)-self._lambda*r)
            constraints = [
                Z >= 0,
                Z*np.ones((m,1),dtype=int) == np.ones((n,1), dtype=int),
                r*np.ones((1,m*(m-1)/2), dtype=int) - 2*np.ones((1,n), dtype=int)*(cvx.mul_elemwise(tmp3, (Z*M))) + cvx.sum_entries(tmp3 * tmp3, axis=0) <= 0,
                np.ones((1,n),dtype=int)*Z >= 1,
                log_p*Z-t-np.dot(log_node_card,O)*Z+tmp2+cvx.mul_elemwise(np.power(tmp1,-1), np.ones((1,n), dtype = int)*Z) == 0
            ]
            prob = cvx.Problem(objective, constraints)
            result = prob.solve(solver='SCS',verbose=False)
            prev_Z[0:n,0:m] = Z.value

        return prev_Z, O

    def total_variance(self, new_clusters):
        node_card = self.node_card
        size_prod = self.clique_domain_size(new_clusters)
        clusters_num = len(new_clusters[0])
        return 2 * (clusters_num)**2 * sum(size_prod.values())


    def clique_domain_size(self, given_cliques):
        node_card = self.node_card
        size_prod = dict()
        for k in range(len(given_cliques)):
            c = given_cliques[k]
            attr_domains = [node_card[i] for i in c]
            size_prod[k] = np.prod(attr_domains)
        return size_prod

    """
    Obtain the final result according to the representation matrix, Z, O.
        :param Z:
            the matrix representation for cliques and clusters.
        :param O:
            the matrix representation for attributes and cliques.
    """
    def get_final_result(self, Z, O):
        logging.debug("Generating results...")
        # find the index of max element in each row of Z
        index = map(lambda row: row.tolist().index(max(row)), Z) 
        # get the row length of Z
        m = len(Z[0])
    
        cluster_ls = []
    
        for icols in range(m):
            tmp_index = [i for i, j in enumerate(index) if j == icols]
            if tmp_index is not None and len(tmp_index) != 0:
                tmp = np.sum(O[:,tmp_index],axis=1)
                cluster_ls.append([i for i, j in enumerate(tmp) if j > 0])
            
        return cluster_ls

    def print_result(self, cluster_ls):
        result = [' '.join(str(col+1) for col in cluster) for cluster in cluster_ls]
        contents = '\n'.join(result)
        if self.output_file != None:
            with open(self.output_file, 'w+') as fout:
                fout.writelines(contents)
        else:
            print contents


def usage():
    print 'python marginal-optimization.py <input_file> <cluster_number> <lambda> <output_file> <random_seeds_num>'

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 6:
        usage()
        sys.exit(1)

    input_file = sys.argv[1]
    clusters_num = int(sys.argv[2])
    _lambda = sys.argv[3]
    output_file = sys.argv[4]
    seed_num = int(sys.argv[5])

    p = Pool(20)
    marginal_opt = MarginalOptimization(input_file, clusters_num, _lambda, output_file)
    cliques_num = marginal_opt.cliques_num

    random_seeds = [np.random.rand(cliques_num, clusters_num) for i in range(seed_num)]
    logging.info("Random seed number: %d" % (seed_num))
    logging.info("Starting to merge cliques...")

    start = time.time()
    z_and_o = p.map(marginal_opt.marginal_optimization, random_seeds)

    merged_clusters = map(lambda zo: marginal_opt.get_final_result(zo[0], zo[1]), z_and_o)
    variances_ls = map(lambda margin: marginal_opt.total_variance(margin), merged_clusters)
    margin_variance = zip(variances_ls, merged_clusters)
    sorted_var = sorted(margin_variance, key = operator.itemgetter(0))
    opt_cluster = sorted_var[0][1]

    end = time.time()
    logging.info("Cliques merge job done in %d seconds." % (end-start))

    marginal_opt.print_result(opt_cluster)
    logging.info("Cliques Merged result printed!")

