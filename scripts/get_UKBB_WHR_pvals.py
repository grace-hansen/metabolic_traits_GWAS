#!/usr/bin/python
import pandas as pd, argparse, gzip, subprocess
parser = argparse.ArgumentParser()
parser.add_argument("trait", help="trait for sumstats")
args = parser.parse_args()

##################
trait=args.trait

def get_pvals(trait,sex):
	sumstats=pd.read_csv("/scratch/midway2/gthansen/metabolic_traits_GWAS/%s.%s.tsv.bgz"%(trait,sex),sep='\t',compression='gzip')
	sumstat_cols=sumstats.columns
	variants=pd.read_csv("/project2/nobrega/grace/metabolic_traits_GWAS/WHR_F_UKBBvariants.txt",sep='\t')
	out=variants.merge(sumstats,on="variant")
	out.to_csv("/project2/nobrega/grace/metabolic_traits_GWAS/pvals/%s.%s_WHR_F.txt"%(trait,sex),sep='\t',index=False)

get_pvals(trait,"both_sexes")
get_pvals(trait,"female")
get_pvals(trait,"male")