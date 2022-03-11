__author__ = 'anastasia'

import pandas as pd
import networkx as nx
import os
import matplotlib.pyplot as plt

from data import Data
from pathways import Network

class Analyses():

    def __init__(self):
        self.data = Data()
        self.data.readParametersFile()
        self.createGraphFile()
        self.df_network = self.data.df_network[self.data.df_network['score'] >= self.data.CAR_threshold]
        self.df_network = self.df_network[self.df_network['max_bridgit'] >= self.data.bridgit_threshold]

    def drawGraphCheck(self):
        with open(self.data.path_compounds_file) as f:
            # list of compounds around which we expand the network
            nodes_init = [int(i.strip()) for i in f.readlines()]

        subG = self.G.subgraph(nodes_init)
        pos=nx.spring_layout(subG)
        nx.draw(subG, pos)
        nx.draw_networkx_labels(subG, pos)
        plt.show()

    def analyseRemovingCompounds(self):
        """
        This function analyses how many compounds are how many steps away after removing each of the shikimate pathway compounds
        :return:
        """

        dict_comp_gen = dict()

        with open('../projects/'+self.data.projectname+'/excludelists/compounds_one_by_one.txt') as f:
            # list of LCSB compound IDs of compounds that have to be removed
            self.compoundsExclude = [int(i.strip()) for i in f.readlines()]

        with open(self.data.path_compounds_file) as f:
            # list of compounds around which we expand the network
            nodes_init = [int(i.strip()) for i in f.readlines()]

        self.nodes_current = list(set(nodes_init).intersection(self.G.nodes())) # filtering only the compounds that are in the graph
        nodes_init_in_graph = self.nodes_current

        self.compoundsExclude = list(set(self.compoundsExclude).intersection(self.G.nodes()))

        if len(self.nodes_current)==0:
            print('the provided compounds are not in graph')
            exit()

        for cmp_init in nodes_init_in_graph:
            print(cmp_init)
            compoundsExcludeTemp = set(self.compoundsExclude) - {cmp_init} # do not exclude the compound around which we expand! => diagonal empty
            for cmp_exc in compoundsExcludeTemp:
                self.nodes_current = [cmp_init]
                G_temp = self.G.copy()
                G_temp.remove_node(cmp_exc) # remove the compound

                self.nodes_total = []
                self.nodes_total.extend(self.nodes_current)

                # expansion ->
                for i in range(self.data.generations_expansion+1):
                    # for how many generations is expanded defined in the parameter file
                    self.getNeighborsGraph(G_temp)
                    dict_comp_gen[(cmp_init,cmp_exc,i)]=len(self.nodes_current)

        # writing the output table
        dict_out = dict()
        for cmp_init in nodes_init_in_graph:
            dict_out[cmp_init]=dict()
            compoundsExcludeTemp = set(self.compoundsExclude) - {cmp_init}
            for cmp_exc in compoundsExcludeTemp:
                str_gen = []
                for i in range(self.data.generations_expansion+1):
                    str_gen.append(str(dict_comp_gen[cmp_init,cmp_exc,i]))
                cmpgens = '->'.join(str_gen)
                dict_out[cmp_init][cmp_exc] = cmpgens

        df = pd.DataFrame.from_dict(dict_out)
        df.to_csv('../output/'+self.data.projectname+'/compound_sensitivity.csv')

    def analyseRemovingReactions(self):

        dict_comp_gen = dict()

        with open('../projects/'+self.data.projectname+'/excludelists/reactions_one_by_one.txt') as f:
            # list of LCSB compound IDs of compounds that have to be removed
            reactionsExclude = [int(i.strip()) for i in f.readlines()]

        with open(self.data.path_compounds_file) as f:
            # list of compounds around which we expand the network
            nodes_init = [int(i.strip()) for i in f.readlines()]

        self.nodes_current = list(set(nodes_init).intersection(self.G.nodes())) # filtering only the compounds that are in the graph
        nodes_init_in_graph = self.nodes_current

        if len(self.nodes_current)==0:
            print('the provided compounds are not in graph')
            exit()

        for cmp_init in nodes_init_in_graph:
            print(cmp_init)
            for rxn_excl in reactionsExclude:
                self.nodes_current = [cmp_init]

                df_reactions_temp = self.data.df_reactions[self.data.df_reactions["rxnUID"]!=rxn_excl]
                edges_uids_keep = set(df_reactions_temp['UID of pair'].to_list())
                df_temp_network = self.df_network[self.df_network['UID of pair'].isin(edges_uids_keep)]
                G_temp = self.createTempGraphFile(df_temp_network)

                self.nodes_total = []
                self.nodes_total.extend(self.nodes_current)

                # expansion ->
                for i in range(self.data.generations_expansion+1):
                    # for how many generations is expanded defined in the parameter file
                    self.getNeighborsGraph(G_temp)
                    dict_comp_gen[(cmp_init,rxn_excl,i)]=len(self.nodes_current)

        # writing the output table
        dict_out = dict()
        for cmp_init in nodes_init_in_graph:
            dict_out[cmp_init]=dict()
            for rxn_excl in reactionsExclude:
                str_gen = []
                for i in range(self.data.generations_expansion+1):
                    str_gen.append(str(dict_comp_gen[cmp_init,rxn_excl,i]))
                cmpgens = '->'.join(str_gen)
                dict_out[cmp_init][rxn_excl] = cmpgens

        df = pd.DataFrame.from_dict(dict_out)
        df.to_csv('../output/'+self.data.projectname+'/reaction_sensitivity.csv')

    def getNeighborsGraph(self, G):

        nodes = []
        # extracting network only around neighbors of the main pathway in the main component
        for id, node in enumerate(self.nodes_current):
            nodes.extend(list(G.neighbors(node)))
        self.nodes_current = list(set(nodes) - set(self.nodes_total))
        self.nodes_total.extend(self.nodes_current)

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

    def createTempGraphFile(self, df_network_temp):

        tempG = nx.Graph()

        print("Load network...")

        for index, row in df_network_temp.iterrows():
            tempG.add_edges_from([(int(row['source']), int(row['target']),
                                        {'id': row['UID of pair'], 'car': row['score'],
                                         'dist': row['dist'], 'dist_known':row['dist_known'],
                                         'dist_exp': row['dist_exp'], 'dist_exp_known':row['dist_exp_known']})])

        return tempG

a = Analyses()
a.analyseRemovingCompounds()
a.analyseRemovingReactions()
#a.drawGraphCheck()
