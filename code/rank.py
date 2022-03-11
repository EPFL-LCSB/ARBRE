__author__ = 'anastasia'
import pandas as pd

#length,cars,avr_car,known_steps,total_known_steps,percent_known,pathway_intermediates,pathway_reactions,pathway_enzymes,bridgit_scores,avg_br_score

class Rank():
    def __init__(self, dat):
        self.data = dat
        self.pw_enzymes_df = pd.read_csv(self.data.pathway_enzymes_file)
        if len(self.pw_enzymes_df) > 0:
            self.shortest_PL = min(self.pw_enzymes_df['length'].to_list())
        else:
            self.shortest_PL = 0

    def calculatePathwayScore(self, row):
        BR = row['avg_br_score'] # average BridgIT score
        CAR = row['avr_car'] # average conserved atom ratio
        PL=row['length'] # length of the pathway
        KR=row['percent_known'] # percent of known reactions in the pathway

        pathway_metric=self.data.weight_BR*BR+self.data.weight_CAR*CAR+self.data.weight_PL*(self.shortest_PL)/PL+self.data.weight_KR*KR

        return '%.3f'%pathway_metric

    def getPathwayRankForAll(self):
        if len(self.pw_enzymes_df) > 0:
            self.pw_enzymes_df['pathway_score']=self.pw_enzymes_df.apply(self.calculatePathwayScore, axis=1)
            self.pw_enzymes_df.sort_values(by='pathway_score', ascending=False, inplace=True)
            self.pw_enzymes_df.to_csv(self.data.ranked_pathway_file, sep=',', index=False)
        else:
            self.pw_enzymes_df.to_csv(self.data.ranked_pathway_file, sep=',', index=False)