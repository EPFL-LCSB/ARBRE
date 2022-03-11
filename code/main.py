__author__ = 'anastasia'
import os

from data import *
from pathways import *
from reactions import *
from enzymes import *
from rank import *
from networkexpansion import *

if __name__ == '__main__':

    dat = Data()
    dat.readParametersFile()

    # exclude unbalanced reactions if set to do so in the parameters
    if dat.exclude_unbalanced == 1:
        # calculate balance for each reaction from the compound formulas
        if dat.calculate_balance == 1:
            try:
                from balance import *
            except ImportError:
                pass
            bal = Balance()
            if not os.path.exists('../data/reaction_balance.csv'):
                bal.createBalanceFile()
        dat.getBalancedReactionsDf()

    # do network expansion if pathways_or_networkexp parameter = 'n'
    if dat.pathways_or_networkexp == 'n':
        net = Network(dat)  # load network
        nexp =  NetworkExpansion(dat, net.G) # create instance of the network expansion class and pass the data and parameters
        nexp.runExpansion() # run network expansion
        exit()

    # do pathway search if the network expansion was not requested
    if 1 in dat.stages:
        dat.findPrecursorAndTargetLCSBID()

        net = Network(dat)

        print("Size of the network loaded: ", len(net.G))

        path = Pathways(dat, net.G)
        path.findPathwaysSet()

    if 2 in dat.stages:
        rxn = Reactions(dat)
        rxn.getAllPairs()
        rxn.assignPairUidToPairs()
        rxn.writePathwaysAsReactionsPw()

    if 3 in dat.stages:
        enz=Enzymes(dat)
        enz.getAllReactions()
        enz.assignECtoReactions()
        enz.writePathwaysAsEnzymesPw()

    if 4 in dat.stages:
        rank = Rank(dat)
        rank.getPathwayRankForAll()

    if 5 in dat.stages:
        try:
            from visualize import *
            vis = Visualization(dat)
            vis.printCompoundImages()
            vis.drawPathways()
        except ImportError as e:
            print("Visualization could not be done due to the following error:")
            print(e)
            if e.msg == "No module named 'rdkit'":
                print("Please, install rdkit and use rdkit environment for running ARBRE with visualization")
                print("First install anaconda: https://docs.anaconda.com/anaconda/install/index.html")
                print("Then install rdkit as described here: https://www.rdkit.org/docs/Install.html")
            pass

    if os.path.exists(dat.pathway_file) \
        and  os.path.exists(dat.pathway_reactions_file)  \
        and os.path.exists(dat.pathway_enzymes_file) \
        and os.path.exists(dat.ranked_pathway_file):
        os.remove(dat.pathway_file)
        os.remove(dat.pathway_reactions_file)
        os.remove(dat.pathway_enzymes_file)


