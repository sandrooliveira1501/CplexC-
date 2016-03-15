#!/usr/bin/python2.6
# -*- coding: iso-8859-1 -*-

 #   This software was developed by Zanoni Dias, Ulisses Dias
 #
 #   It should not be redistributed or used for any commercial purpose
 #   without written permission from authors
 #
 #   release date: nov 15, 2011
 #
 # This software is experimental in nature and is
 # supplied "AS IS", without obligation by the authors to provide
 # accompanying services or support.  The entire risk as to the quality
 # and performance of the Software is with you. The authors
 # EXPRESSLY DISCLAIM ANY AND ALL WARRANTIES REGARDING THE SOFTWARE,
 # WHETHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES
 # PERTAINING TO MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
 # 
 #

 # If you use this softwore anytime in your work, please cite the
 # following paper:
 # 
 # DIAS, U. ; DIAS, Z. . An Improved 1.375-Approximation Algorithm for
 # the Transposition Distance Problem. In: International Conference on
 # Bioinformatics and Computational Biology (ACM-BCB'2010), 2010,
 # Niagara Falls, NY, USA. Proceeding of the 1st ACM International
 # Conference on Bioinformatics and Computational Biology
 # (ACM-BCB'2010). New York, NY, USA : ACM, 2010. p. 334-337.

import sys
import copy
import math
import elias2005

sys.setrecursionlimit(2000)
class dias2010 :
    def __init__(self, param_permutation) :
        self.input = param_permutation
        aux =  param_permutation.split(",")
        permutation = []
        for item in aux :
            permutation.append(int(item))
        self.graph = cycle_graph(permutation)        


    def analyze_transposition(self, cycle_graph, i,j,k) :
        inner_graph = copy.deepcopy(cycle_graph)
        inner_graph.transposition(i, j, k)
        (next_step, next_trans) = self.sort_steps(inner_graph, want_trace = True)
        return (next_step, inner_graph, next_trans)

    def sort(self) :
        distance = 0
        sequence = []
        self.trace = ""
        graph = self.graph
        while not graph.is_ordered() :
            try :
                transpositions = self.sort_steps(graph)
                for transposition in transpositions :
                    distance = distance + 1
                    (i,j,k) = transposition
                    sequence.append(transposition)
                    graph.transposition(i,j,k)
            except :
                print "Erro %s" % graph.__str__()
                print "dias2010_1.375 - %s - %i - %s" % (self.input, 
                                                          distance, 
                                                          sequence)
                sys.exit() 
        print "dias2010_1.375 - %s - %i - %s" % (self.input, 
                                                  distance, 
                                                  sequence)

    def sort_steps(self, cycle_graph, want_trace = False) :
        # Choose the triple that can modify a bad oriented
        # cycle into an oriented cycle allowing a valid 2-move.
        def choose_bad_oriented_modifying_triple(triples, bad_oriented_cycles) :
            modifying_triples = []
            for triple in triples :
                [i,j,k] = triple                
                for bad_cycle in bad_oriented_cycles :
                    [x,y,a,z,b] = bad_cycle.get_cycle()
                    test1 = False
                    test2 = False
                    for el in [x,y,a,z,b] :
                        if i.index < el.index < j.index :
                            test1 = True
                        if j.index < el.index < k.index :
                            test2 = True
                    if (
                        (not (i.index < a.index and
                              k.index > x.index)) and
                         test1 and test2) :
                        modifying_triples.append([i,j,k])
            return modifying_triples

        # Choose the best triple based on the lookahead
        # Sometimes two lookaheds are necessary (in case of draw)

        def choose_best_triple(triples, edges, graph = cycle_graph, depth = 2) :
            return choose_best_triple_1(triples, edges)
            local_depth = depth 
            local_triples = []
            for triple in triples :
                [a,b,c] = [triple[0].index, triple[1].index, triple[2].index]
                local_triples.append([graph, [a,b,c], [a,b,c]])
            best_move = []

            while local_depth > 0 :
                local_depth = local_depth - 1
                scores = []
                for triple in local_triples :
                    [i,j,k] = [triple[2][0],triple[2][1],triple[2][2]]
                    scores.append(self.analyze_transposition(triple[0], i,j,k))

                best_triple = 0
                    
                for i in range(len(scores)) :
                    [a,b,c] = local_triples[i][1]

                    if scores[i][0] < scores[best_triple][0] :
                        best_triple = i
                            
                best_move = local_triples[best_triple][1]
                if scores[best_triple][0] == 8 :
                    break

                next_triples = []
                for i in range(len(scores)) :
                    if scores[i][0] == scores[best_triple][0] :

                        next_triples.append([scores[i][1],local_triples[i][1],scores[i][2]])

                if len(next_triples) < 2 :
                    break
                local_triples = next_triples

            return [[best_move[0], 
                     best_move[1], 
                     best_move[2]]]

        # Choose the best triple based on the number of 
        # of edges that are turned right
        def choose_best_triple_1(triples, edges) :
            # Classifying triples
            scores = []
            for triple in triples :
                [i,j,k] = [triple[0].index,triple[1].index,triple[2].index]
                score = 0
                for edge in edges :
                    if edge[0] < edge[1] :
                        if i < edge[0] <= j < edge[1] <= k :
                            score = score - 1
                    else :
                        if i < edge[1] <= j < edge[0] <= k :
                            score = score + 1
                scores.append(score)

            # Choosing the best triple
            best_triple = 0            
            for i in range(len(scores)) :
                [a,b,c] = triples[i]

                if scores[i] > scores[best_triple] :
                    best_triple = i

            transposition = triples[best_triple]
            return [[transposition[0].index, 
                     transposition[1].index, 
                     transposition[2].index]]


        valid_2_moves = []
        good_0_moves  = []
        oriented_cycles = cycle_graph.get_oriented_cycles()
        bad_oriented_cycles = cycle_graph.get_bad_oriented_cycles()


        ## STEP 1: Verify if there is one triple shuffling a nonoriented cycle
        for cycle in oriented_cycles :
            triples =  cycle_graph.find_oriented_triple(cycle, all = True)

            for triple in triples :
                [i,j,k] = triple
                # If the valid 2-move is interleaved with a nonoriented long cycle
                if cycle_graph.is_shuffling_transposition(triple[0], triple[1], triple[2]) :
                    if want_trace :
                        return [1,[i.index, j.index, k.index]]
                    self.trace = self.trace + "1,"
                    return [[i.index, j.index, k.index]]
                valid_2_moves.append(triple)

        ## STEP 2: Verify if there is one triple transforming a bad
        ## oriented cycle into a good oriented cycle
        transforming_triples = choose_bad_oriented_modifying_triple(valid_2_moves,
                                                                   bad_oriented_cycles)
        if len(transforming_triples) > 0 :
            [i,j,k] = transforming_triples[0]
            if want_trace :
                return [2,[i.index,j.index,k.index]]
            self.trace = self.trace + "2,"
            return [[i.index, j.index, k.index]]
                



        edges = cycle_graph.get_edges()

        ## STEP 3: We cannot guarantee the next valid 2-move, so try
        ## the best triple. The best triple is the one which generates
        ## the most number of right edges.

        # There is a valid 2-move, so, lets choose one
        if len(valid_2_moves) > 0 :
            if want_trace :
                triple = choose_best_triple_1(valid_2_moves,edges)
                return [3, triple[0]]
            self.trace = self.trace + "3,"
            return choose_best_triple(valid_2_moves,edges)

        ## Next: There is no valid 2-move, lets find a good 0-move
        ## STEP 4: Choose the triple that interleaves with a
        ## nonoriented cycle and generates the most number of right
        ## edges
        (shuf_good_0_moves, non_shuf_good_0_moves) = cycle_graph.find_good_0_moves()

        if len(shuf_good_0_moves) > 0 :
            if want_trace :
                triple = choose_best_triple_1(shuf_good_0_moves, edges)
                return [4, triple[0]]
            self.trace = self.trace + "4,"
            return choose_best_triple(shuf_good_0_moves, edges)

        ## STEP 5: Verify if there is one non shuffling good 0-move
        ## transforming a bad oriented cycle into a good oriented
        ## cycle
        transforming_triples = choose_bad_oriented_modifying_triple(non_shuf_good_0_moves,
                                                                    bad_oriented_cycles) 
        if len(transforming_triples) > 0 :
            if want_trace :
                triple = choose_best_triple_1(transforming_triples, edges)
                return [5,triple[0]]
            self.trace = self.trace + "5,"
            return choose_best_triple(transforming_triples, edges)



        ## STEP 6 : We cannot guarantee the next valid 2-move, so try
        ## the best triple. The best triple is the one which generates
        ## the most number of right edges.
        if len(non_shuf_good_0_moves) > 0 :
            if want_trace :
                triple = choose_best_triple_1(non_shuf_good_0_moves, edges)
                return [6,triple[0]]
            self.trace = self.trace + "6,"
            return choose_best_triple(non_shuf_good_0_moves, edges)

        ## STEP 7: There are at least two oriented cycles, so we can
        ## apply a (4,3)-sequence.
        if len(bad_oriented_cycles) > 1 :
            self.trace = self.trace + "7,"

            [cx,cy,ca,cz,cb] = bad_oriented_cycles[0].get_cycle()
            [dx,dy,da,dz,db] = bad_oriented_cycles[1].get_cycle()
                
            for triple in [[cy,cz,cx],[ca,cb,cy],[cb,cy,cx]] :
                [i,j,k] = [triple[0].index, triple[1].index, triple[2].index]
            
                if (i < db.index < j < dz.index < k < dx.index or
                    db.index < i < dz.index < j < dx.index < k or 
                    i < da.index < j < db.index < k < dx.index) : 
                    if want_trace : 
                        return [7,[i,j,k]]
                    return [[i,j,k]]
                
            for triple in [[dy,dz,dx], [da,db,dy],[db,dy,dx]] :
                [i,j,k] = [triple[0].index, triple[1].index, triple[2].index]
                if (i < cb.index < j < cz.index < k < cx.index or
                    cb.index < i < cz.index < j < cx.index < k or
                    i < ca.index < j < cb.index < k < cx.index) :
                    if want_trace :
                        return [7,[i,j,k]]
                    return [[i,j,k]]

            if (ca.index < da.index) :
                if want_trace :
                    return [7,[ca.index,da.index,dz.index]]
                return [[ca.index, da.index, dz.index]]
            else :
                if want_trace :
                    return [7,[da.index,ca.index,cz.index]]
                return [[da.index, ca.index, cz.index]]
        
        if want_trace : 
            return [8,[0,0,0]]
            
        ## STEP 8: There is no valid 2-moves nor good 0-moves
        ## available. So, it is better to use the Elias and Hartman
        ## algorithm
        self.trace = self.trace + "8,"
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        cycle_graph.transform_into_simple_permutation()
        transpositions = cycle_graph.get_database_move()
        if transpositions :
            return transpositions

class cycle_graph(elias2005.cycle_graph) :
    def __init__(self, permutation) :
        n = len(permutation)
        self.n = n
        self.permutation = permutation
        node_list_value = [-1 for x in range(n+2)]
        node_list_index = [-1 for x in range(n+2)]

        # Lemma 9: creating ap and ab
        previous_node = cycle_graph_node(0,0)
        node_list_value[0] = previous_node
        node_list_index[0] = previous_node

        for i in range(n) :
            local_node = cycle_graph_node(i+1, int(permutation[i]))
            node_list_value[int(permutation[i])] = local_node
            node_list_index[i+1] = local_node
            previous_node.set_pointers(ab = local_node)
            local_node.set_pointers(ap = previous_node)
            previous_node = local_node

        local_node = cycle_graph_node(n+1, n+1)
        node_list_value[n+1] = local_node
        node_list_index[n+1] = local_node
        previous_node.set_pointers(ab = local_node)
        local_node.set_pointers(ap = previous_node)

        # Lemma 10: creating ac
        for i in range(n+1) :
            node = node_list_value[i]
            node.ac = node_list_value[node.value + 1]
            node.ac.rac = node
        
        self.begin_index = node_list_index[0]
        self.end_index = node_list_index[-1]
        self.begin_value = node_list_value[0]
        self.end_value = node_list_value[-1]
            
        self.decompose_cycles(node_list_index)

    def decompose_cycles(self, node_list_index) :
        def is_oriented(cycle) :
            previous = cycle[0]
            for count in range(1,len(cycle)) :
                if previous < cycle[count] :
                    return 1
                previous = cycle[count]
            return 0

        n = self.n
        self.clean_visit()
        # Theorem 5: decomposing cycles
        cycles = self.get_cycles()
                
        for i in range(len(cycles)) :
            cycle = cycles[i]
            ccycle = len(cycle)
            tcycle = is_oriented(cycle)
            ncycle = i+1
            cbegin = node_list_index[cycle[-1]]
            cend = node_list_index[cycle[0]]
            for element in cycle :
                node_list_index[element].ccycle = ccycle
                node_list_index[element].tcycle = tcycle
                node_list_index[element].ncycle = ncycle
                node_list_index[element].cbegin = cbegin
                node_list_index[element].cend = cend

    #######################################################################
    ######################## Metodos necessarios ##########################
    #######################################################################
    #ap     : points to the record that stores pi_{i-1}, 1 <= i <= n+1
    #ab     : points to the record that stores pi_{i+1}, 0 <= i <= n
    def transposition(self, i, j, k) :
        node = self.begin_index
        size = self.n + 2
        count = 0
        node_i = 0
        node_j = 0
        node_k = 0
        while node :
            if count == i :
                node_i = node
            if count == j :
                node_j = node
            if count == k :
                node_k = node
            node = node.ab
            count = count + 1

        node_j.ap.ab = node_k
        node_k.ap.ab = node_i
        node_i.ap.ab = node_j

        aux = node_j.ap
        node_j.ap = node_i.ap
        node_i.ap = node_k.ap
        node_k.ap = aux


        node_list_index = []
        node = self.begin_index
        count = 0
        while node:
            node.index = count
            node_list_index.append(node)
            node = node.ab
            count = count + 1

        self.decompose_cycles(node_list_index)

    #######################################################################
    ######################### Metodos auxiliares## ########################
    #######################################################################
    def clean_visit(self) :
        node = self.begin_index
        while node :
            node.visit = 0
            node = node.ab

    def is_ordered(self) :
        node = self.begin_index
        count = 0
        while node :
            if node.value != count :
                return False
            node = node.ab
            count = count + 1
        return True
            
    def print_nodes(self) :
        if self.begin_index:
            print("begin_index = %i" % self.begin_index.index)
        if self.end_index :
            print("end_index   = %i" % self.end_index.index)
            
        node = self.begin_index 
        while node :
            print(node)
            node = node.ab

    def get_cycles(self) :
        cycles = []
        self.clean_visit()
        n = self.n
        node = self.end_index
        while node.ap :
            cycle = []
            local_node = node
            while not local_node.visit :
                cycle.append(local_node.index)
                local_node.visit = 1
                local_node =  local_node.ap.ac
            if len(cycle) :
                cycles.append(cycle)
            node = node.ap
        return cycles

    def __str__(self) :
        str = ""
        node = self.begin_index
        while node :
            str = str + "%s," % node.value
            node = node.ab
        return str

    def get_long_non_oriented_cycle(self) :
        node = self.begin_index
        while node :
            if node.ccycle > 2 and not node.tcycle :
                return node.cend
            node = node.ab
        return 0

    def get_oriented_cycles(self) :
        cycles = []
        self.clean_visit()
        node = self.end_index
        while node :
            if node.ccycle > 2 and node.tcycle and not node.cend.visit  :
                node.visit = 1
                cycles.append(node)
            node = node.ap
        return cycles

    def get_bad_oriented_cycles(self) :
        cycles = []
        self.clean_visit()
        node = self.end_index 
        while node :
            if node.ccycle == 5 and node.tcycle and not node.cend.visit :
                node.visit = 1
                [x,y,a,z,b] = node.get_cycle()
                if a.index < b.index < y.index < z.index < x.index :
                    cycles.append(node)
            node = node.ap
        return cycles
        
        

    # Number of gray edges that must be traversed to go from r to t
    def distance(self, r, t) :
        count = 0
        node = r
        while (node.index != t.index) :
            count = count + 1
            node = node.ap.ac
        return count

    # Split the cycle in 3 based on the transposition
    # i, j and k must be a triple
    def split_cycle(self, i, j, k) :
        if i.ncycle == j.ncycle == k.ncycle :
            splits = [k,i,j]
            cycles = [[],[],[]]
            for count in range(0,3) :
                node = splits[count]
                while node.index != splits[count+1].index :
                    cycles[count].append(node)
        else :
            return [[],[],[]]
        

    #######################################################################
    #################### Metodos retirados de lemas #######################
    #######################################################################


    def find_shuffling_edges(self, i) :
        array_cycle = {}
        j = i.ap
        last_visited = i
        while j.index != 0 :
            if j.ncycle == i.ncycle :
                last_visited = j
            elif j.ccycle > 2 and j.tcycle == 0 : # long non-oriented cycle {C' candidate}
                if not array_cycle.has_key(j.cend) :
                    array_cycle[j.cend] = [last_visited]
                elif array_cycle[j.cend][-1].index != last_visited.index :
                    array_cycle[j.cend].append(last_visited)
            j = j.ap
        return array_cycle

    # Verify if three indices (i,j,k) are a shuffling_transposition 
    def is_shuffling_transposition(self, i, j, k) :
        array_cycle = {}
        last_visited = cycle_graph_node(3000,3000)
        node = self.end_index
        while node.index != 0 :
            if (node.index == i.index or
                node.index == k.index or 
                node.index == j.index) :
                last_visited = node
            elif node.ccycle > 2 and node.tcycle == 0 : # long non-oriented cycle
                if not array_cycle.has_key(node.cend) :
                    array_cycle[node.cend] = [last_visited]
                elif array_cycle[node.cend][-1].index != last_visited.index :
                    array_cycle[node.cend].append(last_visited)
            node = node.ap

        for cend in array_cycle.keys() :
            if len(array_cycle[cend]) > 3 :
                return True
            if  (array_cycle[cend][0].index != 3000 and 
                 len(array_cycle[cend]) == 3) :
                return True
            elif (array_cycle[cend][0].index == 3000 and 
                  len(array_cycle[cend]) == 3  and 
                  i.index <= cend.rac.ab.index) :
                return True
        return False

    # Lemma 13: (zero_move = False) Finds oriented triple with two odd distances
    def find_oriented_triple(self, oriented_cycle, all = False, zero_move = False) :
        all_triples = []
        k = oriented_cycle.cend
        starting = True
        size = oriented_cycle.ccycle
        while k.index != oriented_cycle.cend.index or starting :
            starting = False
            i = k.ap.ac
            dist_k_i = 1
            while i.ap.ac.index != k.index and i.index != k.index :
                j = i.ap.ac
                dist_i_j = 1
                while j.index != k.index :
                    num_odd = ( (dist_i_j % 2) +
                                (dist_k_i % 2) +
                                ((oriented_cycle.ccycle-(dist_k_i + dist_i_j)) % 2) )
                    if (i.index < j.index < k.index and
                        (zero_move or num_odd >= 2)) :
                        if all :
                            all_triples.append([i,j,k])
                        else :
                            return (i, j, k)
                    j = j.ap.ac
                    dist_i_j = dist_i_j + 1

                i = i.ap.ac
                dist_k_i = dist_k_i + 1
            k = k.ap.ac
        if all :
            return all_triples
        else :
            return 0,0,0

    # all cycles on graph must be non-oriented 
    # begin.index < end.index
    # return (x,y) : (D = (x ... c ... d ... y) such (i1, ik) and (c, d) interleaves 
    def find_interleaving_pair(self, begin, end) :
        node = end.ap
        while node.index != begin.index :
            if not node.tcycle :
                if node.index == node.cend.index and node.cbegin.index < begin.index :
                    return (node.cbegin, node.cend)
                if node.index == node.cbegin.index and node.cend.index > end.index :
                    return (node.cbegin, node.cend)
            node = node.ap
        return (0, 0)


    def oriented_even_transposition(self, cycle, triples, transposition) :
        def cmp(triple1, triple2) :
            return triple2[-1] - triple1[-1]

        if cycle.ccycle % 2 == 0 :
            if not transposition :
                return True

            for triple in triples :
                sum_even = 0
                dist_1 = self.distance(triple[0], triple[1])
                dist_2 = self.distance(triple[1], triple[2])
                dist_3 = self.distance(triple[2], triple[0])
                if dist_1 % 2 == 0 :
                    sum_even = sum_even + dist_1
                if dist_2 % 2 == 0 :
                    sum_even = sum_even + dist_2
                if dist_3 % 2 == 0 :
                    sum_even = sum_even + dist_3
                triple.append(sum_even)
            triples.sort(cmp)
            (i,j,k,w) = triples[0][0:4]
            if w > 2 :
                return i, j, k

            # Encontrando ciclos pares internos
            even_cycles = []
            node = cycle.cend
            while node.index != 0 :
                if node.ncycle != cycle.ncycle and node.ccycle % 2 ==0 :
                    even_cycles.append([node.index])
                    inner_node = node.ap.ac
                    while inner_node.index != node.index :
                        even_cycles[-1].append(inner_node.index)
                        inner_node = inner_node.ap.ac
                node = node .ap
            
            # Verificando quem sera o mais intercalado
            for triple in triples :
                Cpares = list(even_cycles)
                (i,j,k) = triple[0:3]
                splits = [k,i,j]
                cycles = [[],[],[]]
                for count in range(0,3) :
                    node = splits[count]
                    while node.index != splits[(count+1)%3].index :
                        cycles[count].append(node.index)
                        node = node.ap.ac
                    if len(cycles[count]) % 2 == 0 :
                        Cpares.append(cycles[count])
                
                weight = 0
                for cycle in cycles :
                    if len(cycle) % 2 != 0 and len(cycle) > 2 :
                        for Cpar in Cpares :
                            Cimpar = list(cycle)
                            Cimpar.sort()
                            count_P = 0 #contador par
                            size_P = len(Cpar)
                            count_I = 0 #contador impar
                            size_I = len(Cimpar)
                            while count_P < size_P and count_I < size_I :
                                if Cpar[count_P] < Cimpar[count_I] :
                                    if Cimpar[count_I] == Cpar[count_P] + 1 :
                                        weight = weight + 1
                                    count_P = count_P + 1
                                else :
                                    if Cpar[count_P] == Cimpar[count_I] + 1 :
                                        weight = weight + 1
                                    count_I = count_I + 1
                triple.append(weight)                       
            triples.sort(cmp)
                            
                    

            (i,j,k) = triples[0][0:3]
            return i, j, k
        else :
            if transposition :
                return (0,0,0)
            else : 
                return False

    def oriented_inner_non_oriented_transposition(self, cycle, triples, transposition) :
        def cmp(triple1, triple2) :
            return triple2[-1] - triple1[-1]
            
        if len(triples) > 0 :
            if not transposition :
                return True
            #for triple in triples :
            #    triple.append(triple[2].index-triple[0].index)
            #triples.sort(cmp)
            right = [[cycle.cbegin, cycle.cend]]
            left = []
            node = cycle.cend.ap.ac
            while node.index != cycle.cend.index :
                if node.index < node.rac.ab.index :
                    left.append([node.rac.ab, node])
                else :
                    right.append([node.rac.ab, node])
                node = node.ap.ac
            
            for triple in triples :
                weight = 0
                (i,j,k) = triple[0:3]
                for edge in left :
                    if i.index <= edge[1].index < j.index < edge[0].index <= k.index :
                        weight = weight + 1

                for edge in right :
                    if i.index < edge[0].index <= j.index <= edge[1].index < k.index :
                        weight = weight - 1
                triple.append(weight)
            triples.sort(cmp)

            (i,j,k) = triples[0][0:3]
            return (i,j,k)
            
        else :
            if transposition :
                return (0,0,0)
            else :
                return False

    def oriented_generic_transposition(self, cycle, triples) :
        if len(triples) > 0 :
            (i,j,k) = triples[0][0:3]
            return (i, j, k)
        else :
            return (0,0,0)


    def oriented_zero_move_crossed_even_transposition(self, cycle, triples, transposition) :
        for triple in triples :
            (y,z,x) = triple[0:3]
            a = z.rac.ab
            b = z.ap.ac
            even_cycles_transpositions = []
            if self.distance(b,x) % 2 == 0 :
                even_cycles_transpositions.append([a,1,y.index,self.end_index.index])
                even_cycles_transpositions.append([b, 1, y.index, x.index])
            else :
                even_cycles_transpositions.append([a, 1, y.index, x.index])
                even_cycles_transpositions.append([b, 1, y.index,  self.end_index.index])
            even_cycles_transpositions.append([z, 1, x.index, self.end_index.index])
            even_cycles_transpositions.append([y, 1, z.index, self.end_index.index])
            even_cycles_transpositions.append([a, 1, b.index, self.end_index.index])

            for ect in even_cycles_transpositions :
                (i,j,k) = self.moving_even_cycles_transposition([ect[0]],ect[1], ect[2],ect[3])
                if i and j and k :
                    if transposition :
                        return (i, j, k)
                    else :
                        return True
        if transposition :
            return (0,0,0)
        else :
            return False

    
    def oriented_zero_move_crossed_oriented_transposition1(self, triples_1, triples_2, transposition) :
        (y, z, x) = triples_1[0][0:3]
        a = z.rac.ab
        b = z.ap.ac
        other_cycle_transpositions = []
        if self.distance(b,x) % 2 == 0 :
            other_cycle_transpositions.append([a,1,y.index,self.end_index.index])
            other_cycle_transpositions.append([b, 1, y.index, x.index])
        else :
            other_cycle_transpositions.append([a, 1, y.index, x.index])
            other_cycle_transpositions.append([b, a.index, y.index,  self.end_index.index])
        other_cycle_transpositions.append([z, b.index, x.index, self.end_index.index])
        other_cycle_transpositions.append([y, a.index, z.index, self.end_index.index])
        other_cycle_transpositions.append([a, 1, b.index, x.index])

        for triple in triples_2 :
            (i,j,k) = triple[0:3]
            for oct in other_cycle_transpositions :
                (move,start,through,until) = oct[0:4]
                if start <= i.index < move.index < j.index < through < k.index <= until :
                    if transposition :
                        return (i, j, k)
                    else :
                        return True
        if transposition :
            return (0,0,0)
        else :
            return False

        
    def oriented_zero_move_crossed_oriented_transposition2(self, triples_1, triples_2, transposition) :
        (y1,z1,x1) = triples_1[0][0:3]
        (y2,z2,x2) = triples_2[0][0:3]

        (a1,b1) = (z1.rac.ab, z1.ap.ac)
        (a2,b2) = (z2.rac.ab, z2.ap.ac)
        transp = []
        if (a1.index < a2.index) :
            transp = [a1,a2,z2]
        elif (a2.index < a1.index) :
            transp = [a2,a1,z1]

        if len(transp) == 3 :
            if transposition :
                return (transp[0], transp[1], transp[2])
            else :
                return True
        else :
            if transposition :
                return (0,0,0)
            else :
                return False

    def oriented_zero_move_generic_transposition(self,cycle, triples, transposition) :
        if not transposition :
            return True
        for triple in triples :
            (y,z,x) = triple[0:3]

            a = z.rac.ab
            b = z.ap.ac
            if b.index < a.index :
                return (b, a, z)
            else :
                return (b, z, x)
        return (0,0,0)


    def find_good_0_moves(self) :
        shuffling_triples = []
        non_shuffling_triples = []
        k = self.end_index
        while k.ap.index != 0 :
            if k.ccycle % 2 == 0 : 
                j = k.ap 
                while j.index != 0 :
                    if j.ccycle % 2 == 0 :
                        i = j.ap
                        while i.index != 0 :
                            if not i.tcycle and i.ccycle % 2 == 0 :
                                if not (k.ncycle == j.ncycle == i.ncycle) :
                                    if j.ncycle == k.ncycle  and self.distance(k,j) % 2 == 1 :
                                        if self.is_shuffling_transposition(i,j,k) :
                                            shuffling_triples.append([i,j,k])
                                        non_shuffling_triples.append([i,j,k])
                                    elif i.ncycle == j.ncycle and self.distance(j,i) % 2 == 1 :
                                        if self.is_shuffling_transposition(i,j,k) :
                                            shuffling_triples.append([i,j,k])
                                        non_shuffling_triples.append([i,j,k])
                                    elif i.ncycle == k.ncycle and self.distance(k,i) % 2 == 1 :
                                        if self.is_shuffling_transposition(i,j,k) :
                                            shuffling_triples.append([i,j,k])
                                        non_shuffling_triples.append([i,j,k])
                            i = i.ap
                    j = j.ap
            k = k.ap
        return (shuffling_triples, non_shuffling_triples)

    def classify_non_oriented_shuffling_triples(self, triples) :
        for triple in triples :
            triple.append([])
            
        node = self.end_index
        self.clean_visit
        while node.index != 0 :
            if not node.tcycle and node.ccycle > 2 and not node.cend.visit :
                node.visit = 1
                for triple in triples :
                    if not (
                        (triple[2].index < node.cbegin.index) or
                        (triple[0].index > node.index) or
                        (triple[2].index > node.index and triple[0].index < node.cbegin.index)) :
                        inner = triple[2].ap
                        found = False
                        while not found and inner.index != triple[1].index :
                            if inner.ncycle == node.ncycle :
                                found = True
                            inner = inner.ap

                        if found :
                            inner = triple[1].ap
                            found = False
                            while not found and inner.index != triple[0].index :
                                if inner.ncycle == node.ncycle :
                                    found = True
                                inner = inner.ap
                            if found :
                                triple[3].append(node)
            node = node.ap
        return triples 

    def classify_non_oriented_triples(self, triples) :
        def cmp(triple1, triple2) :
            if len(triple1[-2]) > 0 and  len(triple2[-2]) > 0 :
                return triple2[-1] - triple1[-1]
            if len(triple1[-2]) == 0 and len(triple2[-2]) == 0 :
                return triple2[-1] - triple1[-1]
            return len(triple2[-2]) - len(triple1[-2])

        node = self.end_index
        right = []
        left = []
        while node.index != 0 :
            if node.index < node.rac.ab.index :
                left.append([node.rac.ab, node])
            else :
                right.append([node.rac.ab, node])
            node = node.ap

        for triple in triples :
            weight = 0
            (i,j,k) = triple[0:3]
            for edge in left :
                if i.index < edge[1].index < j.index < edge[0].index < k.index :
                    weight = weight + 2

            for edge in right :
                if i.index < edge[0].index < j.index < edge[1].index < k.index :
                    weight = weight - 1
            triple.append(weight)
        triples.sort(cmp)
        return triples
            

 
    def classify_non_oriented_crossing_triples(self, triples) :
        def cmp(triple1, triple2) :
            return triple2[-1] - triple1[-1]

        def classify_interleaved_edges(interleaved, shuffled) :
            count_edges  = [] 
            node = interleaved.cend.ap
            open = interleaved.cend
            edges = []
            while node.index >= interleaved.cbegin.index :
                if node.ncycle == shuffled.ncycle :
                    edges.append(node)
                elif node.ncycle == interleaved.ncycle :
                    count_edges.append([open,node,edges])
                    open = node
                    edges = []
                node = node.ap
            return count_edges

        crossed =  {}
        for triple in triples :
            if len(triple[3]) > 0 :
                weight = 0
                half = 0
                for shuffled in triple[3] :
                    if not crossed.has_key(shuffled.ncycle) :
                        node = shuffled.cend
                        inner_crossed = {}
                        while node.index != shuffled.cbegin.index :
                            if node.ncycle != shuffled.ncycle and not node.tcycle and node.ccycle > 2 :
                                if not inner_crossed.has_key(node.ncycle) :
                                    inner_crossed[node.ncycle] = classify_interleaved_edges(node, shuffled)
                            node = node.ap
                        crossed[shuffled.ncycle] = inner_crossed

                    (i,j,k) = triple[0:3]
                    for cross in crossed[shuffled.ncycle] :
                        for item in crossed[shuffled.ncycle][cross] :
                            (b,a,c) = item
                            if len(c) > 0 :
                                # retirando aresta
                                if (a.index < i.index < c[-1].index <= c[0].index < j.index < b.index) :
                                    weight = weight - 1
                                
                                # entrelacada com ciclo criado
                                if  j.index < b.index < k.index :
                                    half = 1

                                if i.index < a.index < j.index :
                                    half = 2

                                if half == 1  and i.ncycle == j.ncycle == b.ncycle :
                                    weight = 3

                                if half == 2 and j.ncycle == k.ncycle == b.ncycle :
                                    weight = 3
                                    
                            else :
                                #entrelacada com ciclo existente
                                if a.index < k.index < b.index :
                                    weight = weight + 3
                                elif i.index < a.index < j.index < b.index:
                                    weight  = weight + 1
                                elif a.index < i.index < b.index < j.index :
                                    weight = weight + 1
                triple.append(weight)
            else :
                triple.append(-10)

        triples.sort(cmp)
        return triples

    # Mover move, passar por through e pesquisar ateh until
    def moving_even_cycles_transposition(self, move, start, through, until) :
        node = self.end_index
        while node.index > until :
            node = node.ap

        after_edges = []
        while node.index > through :
            if node.ccycle % 2 == 0  and not node.tcycle :
                after_edges.append(node)
            node = node.ap

        if len(after_edges) == 0 :
            return (0, 0, 0)

        node = node.ap
        inner_edges = []
        while node.index > move[-1].index :
            if node.ccycle % 2 == 0 and not node.tcycle :
                inner_edges.append(node)
            node = node.ap

        if len(inner_edges) == 0 :
            return (0, 0, 0)
        
        node = move[0].ap 
        while node.index >= start :
            if node.ccycle % 2 == 0 and not node.tcycle :
                i = node
                for j in inner_edges :
                    for k in after_edges :
                        if not (k.ncycle == j.ncycle == i.ncycle) :
                            if j.ncycle == k.ncycle  and self.distance(k,j) % 2 == 1 :
                                return (i, j, k)
                            elif i.ncycle == j.ncycle and self.distance(j,i) % 2 == 1 :
                                return (i, j, k)
                            elif i.ncycle == k.ncycle and self.distance(k,i) % 2 == 1 :
                                return (i, j, k)
            node = node.ap
        return (0, 0, 0)

    def strongly_even_cycles_transposition(self) :
        # direction
        # 0 : backward
        # 1 : foreward
        def find_substitute(through, until, not_in_cycle, direction) :
            def next_node(node, direction) :
                if direction :
                    return node.ab
                else :
                    return node.ap
            node = next_node(through, direction)

            while node and node.index != until :
                if  node.index != 0 and node.ncycle != not_in_cycle and node.ccycle % 2 == 0 :
                    return node
                else :
                    node = next_node(node, direction)
            return 0


        i = self.end_index
        while i.index != 0 :
            if i.ccycle > 2 and i.tcycle == 0 and i.ccycle % 2 == 0: 
                array_cycle = self.find_shuffling_edges(i)
                for cend in array_cycle.keys() :
                    size = len(array_cycle[cend])
                    outer_count = size - 1
                    (x,y) = (0,0)
                    while outer_count >= 0 :
                        (x,y) = (0,0)
                        if outer_count == size - 1 :
                            y = array_cycle[cend][outer_count]
                            if y.index != y.cbegin.index :
                                x = y.ap.ac
                            else :
                                outer_count = outer_count - 1
                            
                        if outer_count != size - 1 :
                            x = array_cycle[cend][outer_count+1]
                            y = array_cycle[cend][outer_count+0]

                        if x.index != y.ap.ac.index :
                            x = y.ap.ac

                        node = y.ap 
                        inner = 0
                        while not inner and node.index != x.index :
                            if node.ncycle == cend.ncycle :
                                inner = node
                            node = node.ap
                            
                        (through, until,direction) = (0,0,0)
                        # Setting the through edge from right and left
                        (through_right, through_left) = (0,0)
                        (until_left, until_right) = (0, self.end_index.index+1)

                        if inner.index == inner.cend.index :
                            through_right = self.end_index
                            until_left = inner.cbegin.index
                            if y.index == y.cend.index :
                                until_left = max(until_left, y.cbegin.index)
                        elif y.index == y.cend.index :
                            through_right = self.end_index
                            until_left = y.cbegin.index
                        elif y.rac.ab.index > inner.rac.ab.index :
                            through_right = y.rac.ab
                        else :
                            through_right = inner.rac.ab
                            
                        if inner.index == inner.cbegin.index :
                            through_left = self.begin_index
                            until_right = inner.cend.index
                            if x.index == x.cbegin.index :
                                until_right = min(until_right, x.cend.index)
                        elif x.index == x.cbegin.index :
                            through_left = self.begin_index
                            until_right = x.cend.index
                        elif x.ap.ac.index < inner.ap.ac.index :
                            through_left = x.ap.ac
                        else :
                            through_left = inner.ap.ac

                        z = find_substitute(through_left, until_left, x.ncycle, 0)
                        if z :
                            return (z,x,y)
                        else :
                            z = find_substitute(through_right, until_right,
                                                    x.ncycle, 1)
                            if z :
                                return (z,x,y)
                        outer_count = outer_count - 1

            i = i.ap
        return 0,0,0

    def even_cycles_transposition(self) :
        cycles = []
        self.clean_visit()

        node = self.end_index
        while node :
            if not node.tcycle and node.ccycle > 2 and not node.cend.visit :
                node.visit = 1
                size = len(cycles)
                if size == 0 :
                    cycles.append(node)
                elif node.cbegin.index < cycles[-1].cbegin.index :
                    if node.ccycle % 2 == 1 or cycles[-1].ccycle % 2 == 0 :
                        cycles.append(node)
                    else :
                        count = size - 1
                        while count >= 0 and cycles[count].ccycle % 2 == 1 :
                            count = count - 1
                        cycles.insert(count, node)
                    cycles.append(node)
                else :
                    count = size - 1
                    while count >= 0 and node.cbegin.index > cycles[count].cbegin.index :
                        count = count - 1
                    cycles.insert(count, node)
            node = node.ap

        for cycle in cycles :
            node = cycle.cbegin
            while node.index != cycle.cend.index :
                #(start, until) = 1, self.end_index.index
                (start,until) = (0,0)
                if node.index == node.cbegin.index :
                    start = 1
                    until = node.cend.index - 1
                else :
                    start = cycle.cbegin.index
                    until = self.end_index.index


                through = node.rac.ab
                while through.index != cycle.cbegin.index :
                    if through.index < until :
                        (i,j,k) = self.moving_even_cycles_transposition([node], 
                                                                        start, through.index, until)
                        if i and j and k :
                            return (i,j,k)
                    through = through.rac.ab
                node = node.rac.ab
        return (0,0,0)
                    
#     # Cycle must be long
    def no_shuffling_transposition(self, cycle) :
        (y, x) = self.find_interleaving_pair(cycle.cbegin, cycle.cend)
        if not (y and x) :
            return (0,0,0)
        s = cycle.cend.ap.ac

        if y.index > s.index :
            (u,v) = self.find_interleaving_pair(cycle.cbegin, s)
            if not (u and v) :
                return (0,0,0)
            if v.index < s.index :
                return (y, x, v)
            elif v.index < y.index :
                return (y, x, u) 
            elif v.index < cycle.cend.index :
                return (u, v, x)
            else :
                return (y, x, u)
        elif x.index < s.index :
            (u,v) = self.find_interleaving_pair(s, cycle.cend)
            if not (u and v) :
                return (0, 0, 0)
            if u.index > x.index :
                return (y, x, v)
            else :
                return (u, v, y)
        elif x.index < cycle.cend.index :
            (u,v) = self.find_interleaving_pair(cycle.cbegin,s)
            if not (u and v) :
                return (0, 0, 0)
            if v.index < s.index :
                return (u, v, x)
            else :
                return (u, v, y)
        else :
            (u,v) = self.find_interleaving_pair(s, cycle.cend)
            if not (u and v) :
                return (0, 0, 0)
            return (u, v, x)
        return (0,0,0)

    def two_cycles_transposition(self) :
        node = self.end_index
        C = 0
        while node :
            if node.ccycle == 2 :
                if not C :
                    C = node
                elif node.ncycle == C.ncycle :
                    C = 0
                else :
                    return (node.cbegin, node.cend, C.cend)
            node = node.ap
        return (0, 0, 0)


class cycle_configuration_graph(cycle_graph) :
    def __init__(self, cycles, shift = 0, mirror = 0) :
        def is_oriented(cycle) :
            previous = cycle[0]
            for count in range(1,len(cycle)) :
                if previous < cycle[count] :
                    return 1
                previous = cycle[count]
            return 0

        self.n = 0
        for cycle in cycles :
            self.n = self.n + len(cycle)

        n = self.n
        self.num_cycles = len(cycles)
          
        # Creating ap and ab
        node_list_index = []
        node_list_index = [cycle_graph_node(i,-1) for i in range(n+1)]
        self.begin_index = node_list_index[0]
        self.end_index = node_list_index[-1]

        for i in range(n) :
            node_list_index[i].ab       = node_list_index[(i+1)]
            node_list_index[(i+1)].ap = node_list_index[i]

        # Creating ac and rac
        for i in range(len(cycles)) :
            cycle = cycles[i]

            ccycle = len(cycle)
            tcycle = is_oriented(cycle)
            ncycle = i+1
            cbegin = node_list_index[cycle[-1]]
            cend = node_list_index[cycle[0]]

            for j in range(ccycle) :
                node_list_index[cycle[j]].ccycle = ccycle
                node_list_index[cycle[j]].tcycle = tcycle
                node_list_index[cycle[j]].ncycle = ncycle
                node_list_index[cycle[j]].cbegin = cbegin
                node_list_index[cycle[j]].cend = cend

                node_list_index[cycle[j]].rac                    = node_list_index[(cycle[(j-1)%ccycle]-1)]
                node_list_index[(cycle[(j-1)%ccycle]-1)].ac = node_list_index[cycle[j]]

    def get_cycles1(self) :
        self.clean_visit()
        node = self.begin_index.ab 
        cycles = []
        while node :
            cycle = []
            local_node = node
            while not local_node.visit :
                cycle.append(local_node.index)
                local_node.visit = 1
                local_node =  local_node.ap.ac
            if len(cycle) :
                cycles.append(cycle)
            node = node.ab
        return cycles


    def is_ordered(self) :
        node = self.begin_index.ab 
        while node :
            if node.ap.ac.index != node.index :
                return False
            node = node.ab
        return True

class cycle_graph_node(elias2005.cycle_graph_node) :

    #value  : stores pi_i
    #index  : stores the black edge i, 0 <= i <= n+1
    #ap     : points to the record that stores pi_{i-1}, 1 <= i <= n+1
    #ab     : points to the record that stores pi_{i+1}, 0 <= i <= n
    #ac     : points to the record that stores i + 1,    0 <= i <= n
    #visit  : indicates if the record has been visited
    #ccycle : stores the length of the cycle to which pi_i belongs
    #tcycle : stores 1 if the cycle is oriented or 0 if it is non-oriented
    #ncycle : is a unique identifier for the cycle
    #pcycle : stores the cycle position of the element in the canonical representation
    #cbegin : stores the black edge i_1 in the canonical representation
    #cend   : stores the black edge i_k in the canonical representation
    def __init__(self, index, value) :
        self.index, self.value = index, value
        self.ap, self.ab, self.ac, self.rac = 0,0,0,0
        self.visit = 0
        self.ccycle, self.tcycle, self.ncycle = 0,0,0
        self.cbegin, self.cend = 0,0

    def set_pointers(self, ap = 0, ab = 0, ac = 0) :
        if ap :
            self.ap = ap
        if ab :
            self.ab = ab
        if ac :
            self.ac = ac
    
    def get_cycle(self) :
        last_node = self.cend
        cycle = [last_node]
        last_node = last_node.ap.ac
        while last_node.index != self.cend.index :
            cycle.append(last_node)
            last_node = last_node.ap.ac
        return cycle

    def __str__(self) :
        str =  "pi_{%i} = %i|{ccycle=%i, tcycle=%i, ncycle=%i}" % (
            self.index, self.value, self.ccycle, self.tcycle, self.ncycle)

        if self.cbegin :
            str = str + ", cbegin=%i" % self.cbegin.index

        if self.cend :
            str = str + ", cend=%i" % self.cend.index

        if self.ab :
            str = str + ", ab=%i" % self.ab.index

        if self.ap :
            str = str + ", ap=%i" % self.ap.index

        if self.ac :
            str = str + ", ac=%i" % self.ac.index

        if self.rac :
            str = str + ", rac=%i" % self.rac.index

        str = str + "}"

        return str

if __name__ == "__main__" :
    dias = dias2010(sys.argv[1])
    dias.sort()

