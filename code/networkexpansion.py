__author__ = 'anastasia'
import os
import networkx as nx
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

class NetworkExpansion():

    def __init__(self, dat, G):
        self.G = G
        self.data = dat

    def runExpansion(self):
        dict_edges_gen = dict()
        dict_comp_gen = dict()
        self.data.getPathCompounds()
        self.nodes_init = self.data.pathCompounds
        G_sub = self.G.subgraph(self.nodes_init)
        edges, edges_id = self.getEdges(G_sub)
        edges_total = []
        edges_total.extend(edges_id)

        for edge in edges_id:
            dict_edges_gen[edge]=0
        for node in self.nodes_init:
            dict_comp_gen[node]=0

        G = self.truncateGraphBasedOnCar()

        self.nodes_current = list(set(self.nodes_init).intersection(G.nodes()))
        if len(self.nodes_current)==0:
            print('the provided compounds are not in graph')
            exit()
        self.nodes_total = []
        self.nodes_total.extend(self.nodes_current)

        for i in range(self.data.generations_expansion):
            G_exp = self.getNeighborsGraph(G)
            edges, edges_id = self.getEdges(G_exp)
            edges_new = list(set(edges_id) - set(edges_total))
            edges_total.extend(edges_new)
            for edge in edges_new:
                dict_edges_gen[edge]=i+1
            for node in self.nodes_current:
                dict_comp_gen[node]=i+1

        G_sub_out = G.subgraph(list(set(self.nodes_total)-{1467865652})) # removing CoA

        ### ANALYSIS OF FILTERING OUT PARTS OF THE NETWORK AROUND SHIKIMATE PATHWAY ###
        #df_comp_shiki_synthesis = pd.read_csv('../output/tyrphetrp_case2a/compound_derivatives.csv')
        #comp_path = df_comp_shiki_synthesis[df_comp_shiki_synthesis['rxn_steps_away']==0]['LCSBID'].to_list()
        #comp_path_plus_one = df_comp_shiki_synthesis[df_comp_shiki_synthesis['rxn_steps_away']==1]['LCSBID'].to_list()
        #G_sub_out = G.subgraph(list(set(self.nodes_total)-{1467865652}-set(comp_path))) # removing CoA and shikimate PW
        #listAA = [1467866080,1467866345, 1467867581,1467866192, 1467866836, 1467865856,1467866374,1467866026, 1467867359, 1467870608, 1467867191, 1467866162,1467866492, 1467866580, 1467867529, 1467866388, 1467868061, 1467866565, 1467866617, 1467867986]
        #G_sub_out = G.subgraph(list(set(self.nodes_total)-{1467865652}-set(comp_path)-set(comp_path_plus_one)-set(listAA))) # removing CoA and shikimate PW plus one

        self.writeOutputAndGephiFile(G_sub_out, dict_comp_gen, dict_edges_gen)

    def truncateGraphBasedOnCar(self):
        # Creating subgraph with edges on the cutoff level
        edges = list(self.G.edges())
        edges_id = [self.G.get_edge_data(edge[0], edge[1])['id'] for edge in edges]
        edges_weights = [self.G.get_edge_data(edge[0], edge[1])['car'] for edge in edges]

        G = nx.Graph()
        for id, edge in enumerate(edges):
            if edges_weights[id] > self.data.CAR_threshold:
                G.add_edge(edge[0], edge[1], id = int(edges_id[id]), car=edges_weights[id])
        return G

    def getEdges(self, Graph):
        edges = list(Graph.edges())
        edges_id = [Graph.get_edge_data(edge[0], edge[1])['id'] for edge in edges]
        return edges, edges_id

    def getNeighborsGraph(self, G):
        nodes = []
        # extracting network only around neighbors of the ain pathway in the main component
        for id, node in enumerate(self.nodes_current):
            nodes.extend(list(G.neighbors(node)))
        self.nodes_current = list(set(nodes) - set(self.nodes_total))
        self.nodes_total.extend(self.nodes_current)
        G_sub_neighbors = G.subgraph(self.nodes_total)
        return G_sub_neighbors

    def writeOutputAndGephiFile(self, G_out, dict_comp_gen, dict_edges_gen):
        #if not os.path.exists('../output/'+self.data.projectname+'/visualization'):
        #    os.mkdir('../output/'+self.data.projectname+'/visualization')
        if not os.path.exists(self.data.gephifilesfolder):
            os.mkdir(self.data.gephifilesfolder)

        outfile = open(self.data.gephifilesfolder+'{}_generations.gdf'.format(self.data.generations_expansion), 'w')

        #processing nodes of the subgraph
        outfile.write("nodedef>name VARCHAR,label VARCHAR,labelvisible BOOLEAN,width DOUBLE,height DOUBLE,color VARCHAR\n")
        listDF = list()

        df_nodes = pd.DataFrame(list(G_out.nodes()), columns=['cUID'])
        cmps_df = self.data.df_compounds.merge(df_nodes, on = 'cUID', how='inner')
        patents_all = cmps_df['NUM_PATENTIDS'].to_list()
        patents_all.sort()
        if len(patents_all)>self.data.num_top_patented_compounds_select:
            top_patent_threshold = patents_all[-(self.data.num_top_patented_compounds_select+1)]
        else: top_patent_threshold = 0

        for node in G_out.nodes():
            dict_comp = dict()
            try:
                cmp_df = cmps_df[cmps_df['cUID']==node].iloc[0]
                if node in self.nodes_init or cmp_df['NUM_PATENTIDS'] > top_patent_threshold:
                    labvisib = 'true'
                else:
                    labvisib = 'false'
                name=cmp_df['COMMON_NAME']
                num_patents = cmp_df['NUM_PATENTIDS']
                if cmp_df['NUM_PATENTIDS'] > 0:
                    patents=np.sqrt(cmp_df['NUM_PATENTIDS'])
                else: patents = 1
            except:
                num_patents = 0
                labvisib = 'false'
                patents = 1

            dict_comp['LCSBID'] = node
            dict_comp['name'] = name
            dict_comp['number_of_patents'] = num_patents
            dict_comp['rxn_steps_away'] = dict_comp_gen[node]
            listDF.append(dict_comp)

            outfile.write("{},'{}','{}',{},{},'{}'\n".format(node, name, labvisib, patents, patents, self.data.colorpalette[dict_comp_gen[node]]))

        df = pd.DataFrame(listDF)
        df.to_csv('../output/'+self.data.projectname+'/compound_derivatives.csv', index = False)
        #processing edges of the subgraph
        outfile.write("edgedef>node1 VARCHAR,node2 VARCHAR,directed BOOLEAN, color VARCHAR, weight DOUBLE\n")

        for edge in G_out.edges():
            edge_id = G_out[edge[0]][edge[1]]['id']
            outfile.write("{},{},{},'{}',{}\n".format(edge[0], edge[1], 'false', dict_edges_gen[edge_id], G_out[edge[0]][edge[1]]['car']))

        outfile.close()