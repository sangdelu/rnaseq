import numpy as np
import pandas as pd

## Class:
class DEtable:
    def __init__(self,path):
        self.path = path
        self.FDR = 0.1
        self.logFC = 1
        self.df = None
        self.up = []
        self.down = []
    
    def set_FDR(self,FDR):
        self.FDR = FDR

    def set_logFC(self,logFC):
        self.logFC = logFC

    def load_table(self):
        try:
            self.df = pd.read_csv(self.path,sep="\t")
        except FileNotFoundError as fnf_error:
            print(fnf_error)
   
    def extract_DEG(self):      
        if self.df is not None:
            self.up = self.df[(self.df["log2FC"] >= self.logFC) & (abs(self.df["FDR"]) <= self.FDR)]["Ensembl"].tolist()
            self.down = self.df[(self.df["log2FC"] <= -self.logFC) & (abs(self.df["FDR"]) <= self.FDR)]["Ensembl"].tolist()
        else:
            print("Table has not been loaded yet.")
            
    def export_genes(self,gene_list, out_path):
        self.df[self.df["Ensemble"] in gene_list].to_csv(out_path,sep="\t",index=False)

