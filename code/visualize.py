__author__ = 'anastasia'

import pandas as pd
import os

from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from PIL import Image

import networkx as nx
import matplotlib.pyplot as plt

class Visualization():
    def __init__(self, dat):
        self.data = dat
        self.df_top_pathways=pd.read_csv(self.data.ranked_pathway_file, nrows=self.data.top_pathways_to_draw)
        self.pw_img_folder = dat.pw_img_folder
        self.folderStructure()
        self.getAllCompounds()
        self.get_smiles_dict()

    def folderStructure(self):
        if not os.path.exists('../temp'):
            os.mkdir('../temp')
        if not os.path.exists('../temp/compound_images'):
            os.mkdir('../temp/compound_images')
        if not os.path.exists('../output/'+self.data.projectname+'/visualization'):
            os.mkdir('../output/'+self.data.projectname+'/visualization')
        if not os.path.exists(self.pw_img_folder):
            os.mkdir(self.pw_img_folder)

    def getAllCompounds(self):
        self.intermediates_all=[]
        for index, row in self.df_top_pathways.iterrows():
            self.intermediates_all.extend([int(i) for i in row['pathway_intermediates'].split('->')])

    def get_smiles_dict(self):
        self.smiles_dict = dict()
        self.patents_dict = dict()
        for intermediate in self.intermediates_all:
            df_comp = self.data.df_compounds[self.data.df_compounds['cUID']==intermediate]
            if len(df_comp) > 0:
                SMILES = df_comp['SMILES'].iloc[0]
                patents = df_comp['NUM_PATENTIDS'].iloc[0]
            else:
                SMILES = ''
                patents=0

            self.smiles_dict[intermediate]=SMILES
            self.patents_dict[intermediate]=patents

    def printCompoundImages(self):
        for cUID, smiles in self.smiles_dict.items():
            self.print_png_from_smiles(smiles, cUID)

    def print_png_from_smiles(self, smiles, refV):
        mol = AllChem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        Draw.MolToFile(mol,"../temp/compound_images/%s.png"%refV,size=(200,250))


    def print_png_from_molfile(self, comp_id, molfile_folder):
        mol = AllChem.MolFromMolFile(molfile_folder + comp_id+'.mol')
        AllChem.Compute2DCoords(mol)
        Draw.MolToFile(mol,"../temp/compound_images/%s.png"%comp_id,size=(200,250))


    # drawing the combined overview of the pathway

    def drawPathways(self):
        for index, row in self.df_top_pathways.iterrows():
            intermediates = [int(i) for i in row['pathway_intermediates'].split('->')]
            self.combine_images_for_one_pathway(index, intermediates)
            self.drawPathwayGraph(row)
            self.addPathwayGraphToImage(index)

    def combine_images_for_one_pathway(self, pathway_id, intermediates):
        compound_images_pathway = ['../temp/compound_images/' + str(refV) +'.png' for refV in intermediates]

        images = [Image.open(x) for x in compound_images_pathway]
        widths, heights = zip(*(i.size for i in images))

        total_width = sum(widths)
        max_height = max(heights)

        new_im = Image.new('RGB', (total_width, max_height))

        x_offset = 0
        for im in images:
            new_im.paste(im, (x_offset,0))
            x_offset += im.size[0]

        new_im.save(self.pw_img_folder+'/pathway_{}.jpg'.format(pathway_id))

    def drawPathwayGraph(self, row):
        G=nx.Graph()

        intermediates = [int(i) for i in row['pathway_intermediates'].split('->')]

        fig = plt.figure(figsize=(row['length']*5, 2))#figsize=(row['length']*1000, 1000), dpi=1
        ax =plt.subplot(1, 1, 1)

        node_patent_sizes = [self.patents_dict[interm]/100 for interm in intermediates]

        cars = [float(i) for i in row['cars'].split('|')]
        knowns = [int(i) for i in row['known_steps'].split('|')]
        enzymes = row['pathway_enzymes'].split('=>')
        reactions = row['pathway_reactions'].split('->')
        enz_rxn = [enzymes[i]+' ('+reactions[i]+')' for i in range(len(enzymes))]

        positions=dict()
        pos_edge_labels = dict()
        for i in range(row['length']):
            if knowns[i]==1: col = 'g'
            else: col = 'black'
            G.add_edge(intermediates[i], intermediates[i+1],color=col,weight=cars[i]*20, label=enz_rxn[i])
            positions[intermediates[i]] = (i, 0)
            pos_edge_labels[intermediates[i]] = (i, 0.02)
        positions[intermediates[i+1]] = ((i+1),0)
        pos_edge_labels[intermediates[i+1]] = ((i+1),0.02)

        nx.draw_networkx_nodes(G, pos=positions, node_size=node_patent_sizes,node_color='#FDA50F')

        edges = G.edges()
        colors = [G[u][v]['color'] for u,v in edges]
        weights = [G[u][v]['weight'] for u,v in edges]
        nx.draw_networkx_edges(G, pos=positions, width=weights, edge_color=colors)

        edge_labels = nx.get_edge_attributes(G, 'label')

        nx.draw_networkx_edge_labels(G, pos_edge_labels, edge_labels, font_size=14)

        #fig.patch.set_visible(False)
        ax.axis('off')

        fig.savefig('graph_path_temp.png', bbox_inches='tight', pad_inches=0)

    def addPathwayGraphToImage(self, pathway_id):

        imgage_pw_file = self.pw_img_folder+'/pathway_{}.jpg'.format(pathway_id)
        image_graph_file = 'graph_path_temp.png'

        image_pw = Image.open(imgage_pw_file)

        basewidth = image_pw.size[0]
        image_graph = Image.open(image_graph_file)

        wpercent = (basewidth/float(image_graph.size[0]))
        hsize = int((float(image_graph.size[1])*float(wpercent)))
        img = image_graph.resize((basewidth,hsize), Image.ANTIALIAS)
        #img.save('gra.jpg')

        images = [image_pw, img]

        widths, heights = zip(*(i.size for i in images))

        total_height = sum(heights)
        max_width = max(widths)

        new_im = Image.new('RGB', (max_width, total_height))

        y_offset = 0
        for im in images:
            new_im.paste(im, (0,y_offset))
            y_offset += im.size[1]

        new_im.save(self.pw_img_folder+'/pathway_with_graph_{}.jpg'.format(pathway_id))

        os.remove('graph_path_temp.png')
        os.remove(self.pw_img_folder+'/pathway_{}.jpg'.format(pathway_id))
