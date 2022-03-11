__author__ = 'anastasia'
from itertools import product
import numpy as np

class Enzymes():
    """
    class enzymes annotates the pathway predictions with enzymes predictions coming from BridgIT
    """
    def __init__(self, dat):
        self.data = dat

    def getAllReactions(self):
        """
        this function collects all the reactions LCSBIDs that have to be annotated with enzymes
        """
        reactions_all = []
        with open(self.data.pathway_reactions_file) as f:
            header = f.readline().strip().split(',')
            for line in f:
                rxn_path = line.strip().split(',')[header.index('pathway_reactions')].split('->')
                for rxn in rxn_path:
                    reactions_all.append(int(rxn))
        self.reactions_all = set(reactions_all)
        print('Total number of reactions in all pathways: ', len(self.reactions_all))

    def assignECtoReactions(self):
        """
        this function creates dictionary with all the possible ECs for each of the reactions of the pathways
        :return:
        """
        self.dict_reactions_enzymes = dict()
        for rxn in self.reactions_all:
            df_enzymes = self.data.df_bridgit_predictions[self.data.df_bridgit_predictions['rxnUID']==rxn]

            if len(df_enzymes) == 0:
                print('no bridgit predictions per reaction!', rxn)
                ec_predictions = []
            elif len(df_enzymes) > self.data.bridgit_pred_per_reaction:
                df_temp=df_enzymes.sort_values(by = ['score'], ascending=False)
                ec_predictions=df_temp['EC'].to_list()[0:self.data.bridgit_pred_per_reaction]
                scores_of_predictions=df_enzymes['score'].to_list()[0:self.data.bridgit_pred_per_reaction]
            else:
                ec_predictions=df_enzymes['EC'].to_list()
                scores_of_predictions=df_enzymes['score'].to_list()

            self.dict_reactions_enzymes[rxn] = list(zip(ec_predictions, scores_of_predictions))

    def writePathwaysAsEnzymesPw(self):
        """
        this function rewrites the pathways with the addition of the enzyme information
        :return:
        """
        with open(self.data.pathway_reactions_file) as f, open(self.data.pathway_enzymes_file, 'w') as w:
            header_line = f.readline().strip()
            header = header_line.split(',')
            w.write(header_line+',pathway_enzymes,bridgit_scores,avg_br_score\n')
            for line in f:
                #print(line)
                pw_ec = []
                rxn_path = line.strip().split(',')[header.index('pathway_reactions')].split('->')
                for rxn in rxn_path:
                    #print(self.dict_reactions_enzymes[int(rxn)])
                    pw_ec.append(self.dict_reactions_enzymes[int(rxn)])

                all_pw = product(*pw_ec)
                #print(all_pw)
                for path in all_pw:
                    scores = [i[1] for i in path]
                    avg_score = np.mean(scores)
                    w.write('{},"{}",{},{}\n'.format(line.strip(),
                                                   '=>'.join([str(i[0]) for i in path]),
                                                   '|'.join([str(i) for i in scores]),
                                                   '%.2f'%avg_score))

