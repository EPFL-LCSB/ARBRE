__author__ = 'anastasia'

import pandas as pd
import networkx as nx
import os
import sys

class Data():
    def __init__(self):

        if len(sys.argv) == 1:
            print('Please, set project name as the first argument of the script "ex. python3 main.py scoulerine_case1a"')
            print('Alternatively, you can introduce the project name now:')
            self.projectname = input()

        elif len(sys.argv) == 2:
            self.projectname = sys.argv[1]

        else:
            print('Please, set project name as the first argument of the script "ex. python3 main.py scoulerine_case1a"')
            print('Finishing the execution')
            exit()

        #self.projectname = "norcoclaurine_precursorsList"
        #self.projectname = "scoulerine_case1a"
        #self.projectname = "benzenecarboperoxic_acid_case1b"
        #self.projectname = "tyrphetrp_case2a"
        #self.projectname = "morphine_case2b"

        self.network_file = '../data/network.csv'
        self.compounds_file = '../data/compounds.csv'
        self.reactions_file = '../data/reactions_pairs.csv'
        self.bridgit_predictions_file = '../data/bridgit_predictions.csv'
        self.source_compounds_file = '../projects/'+self.projectname+'/source_compounds.txt'
        self.path_compounds_file = '../projects/'+self.projectname+'/path_compounds.txt'
        self.excludelist_compounds = '../projects/'+self.projectname+'/excludelists/compounds.txt'
        self.excludelist_reactions = '../projects/'+self.projectname+'/excludelists/reactions.txt'
        self.pathway_file = '../output/'+self.projectname+'/pathways_raw.txt'
        self.pathway_reactions_file='../output/'+self.projectname+'/pathway_reactions.txt'
        self.pathway_enzymes_file='../output/'+self.projectname+'/pathway_enzymes.txt'
        self.ranked_pathway_file = '../output/'+self.projectname+'/pathways.csv'
        self.compound_limits_file = '../projects/'+self.projectname+'/compound_parameters.csv'
        self.gephifilesfolder = '../output/'+self.projectname+'/visualization/'
        #self.pw_img_folder = '../output/'+self.projectname+'/visualization/pathways_img'
        self.pw_img_folder = '../output/'+self.projectname+'/visualization'

        # Setting the defaults in case the individual parameter files were not provided
        if os.path.exists( '../projects/'+self.projectname+'/parameters.txt'):
            self.param_file = '../projects/'+self.projectname+'/parameters.txt'
        else:
            self.param_file = '../default_parameters/parameters.txt'

        if not os.path.exists(self.compound_limits_file):
            self.compound_limits_file = '../default_parameters/compound_parameters.csv'

        if not os.path.exists(self.source_compounds_file):
            self.source_compounds_file = '../default_parameters/source_compounds.csv'

        if not os.path.exists(self.excludelist_compounds):
            self.excludelist_compounds = "../default_parameters/excludelists/compounds.txt"

        if not os.path.exists(self.excludelist_reactions):
            self.excludelist_reactions = "../default_parameters/excludelists/reactions.txt"

        # reading the data into the pandas dataframes
        self.df_network = pd.read_csv(self.network_file)
        self.df_compounds = pd.read_csv(self.compounds_file, low_memory=False)
        self.df_reactions = pd.read_csv(self.reactions_file)
        self.df_bridgit_predictions = pd.read_csv(self.bridgit_predictions_file)
        self.known_steps_all = self.df_network[self.df_network['known_reaction']==1]['UID of pair'].to_list()
        self.df_compound_properties_limits = pd.read_csv(self.compound_limits_file)

        self.colorpalette = {0: "0, 76, 163", 1: "138, 81, 165", 2:"203, 94, 153", 3: "244, 123, 137", 4: "255, 164, 126",
                             5: "255, 210, 134", 6: "255, 255, 166", 7:"205, 204, 0", 8: "255, 255, 0", 9: "255, 204, 0", 10: "255, 102, 0", 11: "204, 0, 1"}

        # creating output directories
        if not os.path.exists('../output'):
            os.mkdir('../output')
        if not os.path.exists('../output/'+self.projectname):
            os.mkdir('../output/'+self.projectname)


    def readParametersFile(self):
        dict_param = dict()
        with open(self.param_file) as f:
            for line in f:
                if not line.startswith('#') and line.strip() != '':
                    dict_param[line.split('|')[0]] = line.split('|')[1].strip()
        self.pathways_or_networkexp = dict_param['pathways_or_networkexp']
        self.source_compound_inchikey = dict_param['source']
        self.target_compound_inchikey = dict_param['target']
        self.num_shortest_pathways = int(dict_param['num_shortest_pathways'])
        self.prefer_known = int(dict_param['prefer_known'])
        self.structure_based_pairs = int(dict_param['structure_based_pairs'])
        self.use_exponential_transformation = int(dict_param['use_exponential_transformation'])
        self.use_source_compounds_list = dict_param['use_source_compounds_list']
        self.stages = [int(i) for i in dict_param['stages'].split(',')]
        self.top_pathways_to_draw = int(dict_param['top_pathways_to_draw'])

        # Network expansion parameters
        self.generations_expansion = int(dict_param['generations_expansion'])
        self.num_top_patented_compounds_select = int(dict_param['num_top_patented_compounds_select'])

        self.calculate_balance = int(dict_param['calculate_balance'])
        self.exclude_unbalanced = int(dict_param['exclude_unbalanced'])
        self.rxn_per_rp_pair = int(dict_param['rxn_per_rp_pair'])
        self.bridgit_pred_per_reaction = int(dict_param['bridgit_pred_per_reaction'])
        self.bridgit_threshold = float(dict_param['bridgit_threshold'])
        self.CAR_threshold = float(dict_param['CAR_threshold'])
        self.weight_BR = float(dict_param['weight_BR'])
        self.weight_CAR = float(dict_param['weight_CAR'])
        self.weight_PL = float(dict_param['weight_PL'])
        self.weight_KR = float(dict_param['weight_KR'])


        with open(self.excludelist_compounds) as f:
           self.excludecompounds = [int(i.strip()) for i in f.readlines() if i!= '\n']

        with open(self.excludelist_reactions) as f:
           self.excludereactions = [int(i.strip()) for i in f.readlines() if i!= '\n']

        # possible distance types: dist, dist_exp, dist_known, dist_exp_known
        if self.prefer_known == 1 and self.use_exponential_transformation == 1:
            self.distance_type = 'dist_exp_known' # known exponential distance
        elif self.prefer_known == 1 and self.use_exponential_transformation == 0:
            self.distance_type = 'dist_known' # known distance
        elif self.prefer_known == 0 and self.use_exponential_transformation == 1:
            self.distance_type = 'dist_exp' # exponential distance
        elif self.prefer_known == 0 and self.use_exponential_transformation == 0:
            self.distance_type = 'dist' # standard distance

        # filter out the reactions that have bridgIT below the threshold
        self.df_bridgit_predictions = self.df_bridgit_predictions[self.df_bridgit_predictions['score'] >= self.bridgit_threshold]
        df_bridgit_top = self.df_bridgit_predictions[self.df_bridgit_predictions['rank']==1]
        self.df_reactions = self.df_reactions.merge(df_bridgit_top, how = 'inner', on='rxnUID')

        # filter out the reactions that are in excludelists/reactions.txt
        print('filter out the reactions that are in excludelists/reactions.txt')
        self.df_reactions = self.df_reactions[~self.df_reactions.rxnUID.isin(self.excludereactions)]
        edges_uids_keep = set(self.df_reactions['UID of pair'].to_list())
        self.df_network = self.df_network[self.df_network['UID of pair'].isin(edges_uids_keep)]
        self.df_network = self.df_network[self.df_network['structure_based']<=self.structure_based_pairs]


    def findPrecursorAndTargetLCSBID(self):
        self.main_precursor = self.getLCSBIDfromINCHIKEY(self.source_compound_inchikey)
        self.main_target = self.getLCSBIDfromINCHIKEY(self.target_compound_inchikey)

    def getSourceCompounds(self):
        with open(self.source_compounds_file) as f:
            self.sourceCompounds = [self.getLCSBIDfromINCHIKEY(i.strip()) for i in f.readlines()]

    def getPathCompounds(self):
        with open(self.path_compounds_file) as f:
            self.pathCompounds = [self.getLCSBIDfromINCHIKEY(i.strip()) for i in f.readlines()]

    def getLCSBIDfromINCHIKEY(self, inchikey):
        comp_match_df = self.df_compounds[self.df_compounds['INCHIKEY']==inchikey]
        if len(comp_match_df) > 0:
            return comp_match_df['cUID'].iloc[0]
        else:
            return 0

    def getBalancedReactionsDf(self):
        """
        Exclude unbalanced reactions from the reaction DF
        :return:
        """
        if self.exclude_unbalanced == 1:
            df_balance = pd.read_csv('../data/reaction_balance.csv')
            num_reactions_before = len((set(self.df_reactions['rxnUID'].to_list())))
            self.df_reactions = self.df_reactions.merge(df_balance, how="inner", on='rxnUID')
            print(len(set(self.df_reactions['rxnUID'].to_list())), " balanced reactions selected out of ", num_reactions_before, "initial reactions")





