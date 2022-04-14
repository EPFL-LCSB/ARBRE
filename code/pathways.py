__author__ = 'anastasia'
import time
import multiprocessing
from itertools import islice
import networkx as nx
import numpy as np
import pandas as pd

class Network():
    def __init__(self, data):
        self.data = data
        self.createGraphFile()

    def createGraphFile(self):
        self.G = nx.Graph()

        print("Load network...")

        df_network = self.data.df_network[self.data.df_network['score'] >= self.data.CAR_threshold]
        df_network = df_network[df_network['max_bridgit'] >= self.data.bridgit_threshold]

        for index, row in df_network.iterrows():
            self.G.add_edges_from([(int(row['source']), int(row['target']),
                                        {'id': row['UID of pair'], 'car': row['score'],
                                         'dist': row['dist'], 'dist_known':row['dist_known'],
                                         'dist_exp': row['dist_exp'], 'dist_exp_known':row['dist_exp_known']})])

        compoundsPassFilter = self.getFittingCompounds()

        self.G = self.getSubgraphCompounds(compoundsPassFilter)

    def getFittingCompounds(self):
        """
        This function truncates the list of compounds to the compounds with the properties defined in the compound_parameters.csv file
        :return: list of LCSB IDs of the compounds that fit the parameter properties
        """
        df_compound_limits = self.data.df_compound_properties_limits
        query_list = []
        for index, row in df_compound_limits.iterrows():
            # e.g.'(Num_atoms>1&Num_atoms<1000) and (MW>0&MW<1000)'
            query_list.append(row['Code']+'>='+str(row['Min'])+'&'+row['Code']+'<='+str(row['Max']))

        query = '&'.join(query_list)

        df_filtered_compounds = self.data.df_compounds.query(query, engine="python")

        cUIDs = df_filtered_compounds['cUID'].to_list()
        cUIDs = list(set(cUIDs)-{1467865652}-{1467865841}-set(self.data.excludecompounds)) # exclude CoA and acetyl CoA and compounds in excludelists/compounds.txt

        return cUIDs

    def getSubgraphCompounds(self, cUIDs):
        H = self.G.subgraph(cUIDs)
        return H


class Pathways():
    def __init__(self, dat, G):
        self.G = G
        self.data = dat

    def findPathwaysSet(self):
        print('Searching for pathways')
        t2 = time.time()
        if self.data.use_source_compounds_list == 'no':
            initialPathways = self.findInitialPathways()
        elif self.data.use_source_compounds_list == 'yes':
            self.data.getSourceCompounds()
            self.findSourceMetabolitesInGraph()
            initialPathways = self.findPathwaysFromSeveralCompounds()
        t3 = time.time()

        print('time to find raw pathways:', t3 - t2)

        # write down the initial pathways set to check later if it is recovered in the final network
        init_pathways_file = open(self.data.pathway_file, 'w')
        init_pathways_file.write('length,cars,avr_car,known_steps,total_known_steps,percent_known,pathway_intermediates\n')
        for path in initialPathways:
            CARs = [self.G[path[i]][path[i+1]]['car'] for i in range(len(path)-1)]
            avg_CAR = np.mean(CARs)
            known_steps = [1 if self.G[path[i]][path[i+1]]['id'] in self.data.known_steps_all else 0 for i in range(len(path)-1)]
            total_known_steps = sum(known_steps)
            percent_known_steps = total_known_steps/(len(path)-1)
            init_pathways_file.write('{},{},{},{},{},{},{}\n'.format(
                len(path)-1,
                '|'.join([str(i) for i in CARs]),
                '%.2f'%avg_CAR,
                '|'.join([str(i) for i in known_steps]),
                total_known_steps,
                '%.2f'%percent_known_steps,
                '->'.join([str(i) for i in path])))
        init_pathways_file.close()

        return initialPathways

    # ===============================================================================
    #
    # Graph search towards the list of selected compounds ---------->

    def findPathwaysFromSeveralCompounds(self):

        ts0 = time.time()
        PathwaysList = []

        # pathway search is parallelized
        pool = multiprocessing.Pool()
        PathwaysList.extend(pool.map(self.findShortestPathwaysFromFromSourceMetabolite, self.sourceMetabolites))
        # pool has to be closed to avoid OSError: [Errno 24] Too many open files
        pool.close()

        #for boundary in total_boundary:
        #    print('Boundary', boundary)
        #    PathwaysList.extend(self.graph.findShortestPathwaysToModel(boundary))
        ts1 = time.time()
        print('Time to find all pathways: ', ts1 - ts0)
        totalPathwaysList = []
        for path in PathwaysList:
            totalPathwaysList.extend(path)

        #self.stats_file.write('{}\t'.format(len(totalPathwaysList))) # number of pathways found for stats
        #self.stats_file.write('{}\t'.format(ts1 - ts0)) #  time to find all pathways for stats

        # write down the initial pathways set to check later if it is recovered in the final network
        # if the file does not exist means here we found an initial pathways set without defined precursor
        init_pathways_file = open(self.data.pathway_file, 'w')
        init_pathways_file.write('length,cars,avr_car,known_steps,total_known_steps,percent_known,pathway_intermediates\n')
        for path in totalPathwaysList:
            CARs = [self.G[path[i]][path[i+1]]['car'] for i in range(len(path)-1)]
            avg_CAR = np.mean(CARs)
            known_steps = [1 for i in range(len(path)-1) if self.G[path[i]][path[i+1]]['id'] in self.data.known_steps_all]
            total_known_steps = sum(known_steps)
            percent_known_steps = total_known_steps/(len(path)-1)
            init_pathways_file.write('{},{},{},{},{},{},{}\n'.format(
                len(path)-1,
                '|'.join([str(i) for i in CARs]),
                '%.2f'%avg_CAR,
                '|'.join([str(i) for i in known_steps]),
                total_known_steps,
                '%.2f'%percent_known_steps,
                '->'.join([str(i) for i in path])))
        init_pathways_file.close()

        return totalPathwaysList

    #############            Functions for pathway search within Graph             ###############

    def setNumberOfPathways(self, num_pathways):
        """
        This function sets number of shortest pathways that will be found within the initial graph
        :param num_pathways:
        """
        self.number_of_pathways = num_pathways

    def findSourceMetabolitesInGraph(self):
        self.sourceMetabolites = set(self.data.sourceCompounds).intersection(set(self.G.nodes()))
        print("Number of defined source metabolites in the network:", len(self.sourceMetabolites))

    def findInitialPathways(self):
        allPathways = self.k_shortest_paths(self.data.main_precursor, self.data.main_target, self.data.num_shortest_pathways, self.data.distance_type)
        return allPathways

    def findShortestPathwaysFromFromSourceMetabolite(self, sourceMetabolite):
        #collecting all shortest pathways to a target from each precursor in the model
        allPathways = []
        allPathways.extend(self.k_shortest_paths(sourceMetabolite, self.data.main_target, self.data.num_shortest_pathways, self.data.distance_type))
        return allPathways

    def k_shortest_paths(self, source, target, k, weight):
        return list(islice(nx.shortest_simple_paths(self.G, source, target, weight=weight), k))
