__author__ = 'anastasia'
from itertools import product

class Reactions():
    def __init__(self, dat):
        self.data = dat

    def getAllPairs(self):
        pairs = []
        with open(self.data.pathway_file) as f:
            header=f.readline().strip().split(',')
            for line in f:
                path = line.strip().split(',')[header.index('pathway_intermediates')].split('->')
                for i in range(len(path)-1):
                    pairs.append((int(path[i]), int(path[i+1])))
        self.pairs_participants = set(pairs)

    def assignPairUidToPairs(self):
        self.dict_pairs_reactions = dict()
        for pair in self.pairs_participants:
            comp1, comp2 = pair[0], pair[1]
            df_pair = self.data.df_network.query('(source==@comp1&target==@comp2) or (source==@comp2&target==@comp1)')
            if len(df_pair) == 0:
                print('participants pair not found', comp1, comp2)
                pairUID = 0
            elif len(df_pair) >1:
                print('more than 1 pairUID for the same participants')
                pairUID = df_pair['UID of pair'].iloc[0]
            else:
                pairUID=df_pair['UID of pair'].iloc[0]

            df_reactionsUID = self.data.df_reactions[self.data.df_reactions['UID of pair']==pairUID]
            if len(df_reactionsUID) == 0:
                print('no reactions per pair!', pair, 'UID:', pairUID)
                reactionsUID = []
            elif len(df_reactionsUID) > self.data.rxn_per_rp_pair:
                df_temp=df_reactionsUID.sort_values(by = ['score'], ascending=False)
                reactionsUID=df_temp['rxnUID'].to_list()[0:self.data.rxn_per_rp_pair]
            else:
                reactionsUID=df_reactionsUID['rxnUID'].to_list()
            self.dict_pairs_reactions[pair] = reactionsUID

    def writePathwaysAsReactionsPw(self):
        with open(self.data.pathway_file) as f, open(self.data.pathway_reactions_file, 'w') as w:
            header_line=f.readline().strip()
            header = header_line.split(',')
            w.write(header_line+',pathway_reactions\n')
            for line in f:
                pw_rxn=[]
                path = line.strip().split(',')[header.index('pathway_intermediates')].split('->')
                for i in range(len(path)-1):
                    pair=(int(path[i]), int(path[i+1]))
                    pw_rxn.append(self.dict_pairs_reactions[pair])
                all_pw = product(*pw_rxn)
                for path in all_pw:
                    w.write('{},{}\n'.format(line.strip(), '->'.join([str(i) for i in path])))

