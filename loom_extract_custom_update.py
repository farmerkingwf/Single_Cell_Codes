import numpy as np
import loompy
import pandas
import glob, os    
import sys

i=sys.argv[1]

print(i)

gene_check = pandas.read_csv("gene_master_list.txt", sep='\t')

sampleid_match="*"+i+"*loom"
sampleid_match=glob.glob(sampleid_match)[0]
ds = loompy.connect(sampleid_match)
#mutid=mut_table[mut_table.TapestriID==i].Variant.unique()
allid=ds.ra.id
df_merged = pandas.DataFrame()
selected_id = []
for j in allid: 
    gene=ds.ra.amplicon[ds.ra.id==j][0]
    batch=gene.split("_")[0]
#    gene=gene.split("_")[1]
    if batch == 'nan':
       gene = 'FLT3'
    elif batch =='customMYE':
       m_id=j.split(":",1)[0]
       if m_id == 'chr2':
          gene = 'DNMT3A'
       if m_id == 'chr4':
          gene = 'TET2'
       if m_id == 'chr17':
          m_id2 = j.split(":",2)[1]
          if int(m_id2) < 58740172:
            gene = 'TP53'
          else:
            gene = 'PPM1D'
    else:
       gene=gene.split("_")[1]

    if (gene in gene_check.loc[gene_check['SampleID'] == i,'gene'].values) or (gene == 'FLT3'):

#        gene=j.split(":")[0]
#        gene=gene.split("_")[0]
#        m_id=j.split(":",1)[1]
#        m_id=m_id.replace("chr","")
        
#        df=pandas.DataFrame(ds.layer[""][np.flatnonzero(np.core.defchararray.find(ds.ra.id,m_id)!=-1),:])[:1]
#        df_ad=pandas.DataFrame(ds.layer["AD"][np.flatnonzero(np.core.defchararray.find(ds.ra.id,m_id)!=-1),:])[:1]
#        df_dp=pandas.DataFrame(ds.layer["DP"][np.flatnonzero(np.core.defchararray.find(ds.ra.id,m_id)!=-1),:])[:1]
#        df_gq=pandas.DataFrame(ds.layer["GQ"][np.flatnonzero(np.core.defchararray.find(ds.ra.id,m_id)!=-1),:])[:1]
      selected_id.append(j)
      df=pandas.DataFrame(ds.layer[""][ds.ra.id == j,:])[:1]
      df_ad=pandas.DataFrame(ds.layer["AD"][ds.ra.id == j,:])[:1]
      df_dp=pandas.DataFrame(ds.layer["DP"][ds.ra.id == j,:])[:1]
      df_gq=pandas.DataFrame(ds.layer["GQ"][ds.ra.id == j,:])[:1]
      if gene == "SRSF2" or gene == "RUNX1" or gene == "DNMT3A" or gene == "MYC" or gene == "ZRSR2":
          df[df_dp <5] =3
          df[(df_ad <2) & (df>0)] =3
          df[(df_dp <100) &  (df_ad/df_dp <0.15) & (df>0)] =3
          df[(df_dp >=100) &  (df_ad/df_dp <0.10) & (df>0)] =3

      else:
          df[df_dp <10] =3
          df[(df_ad <3) & (df>0)] =3
          df[(df_dp <100) &  (df_ad/df_dp <0.15) & (df>0) ] =3
          df[(df_dp >=100) &  (df_ad/df_dp <0.10) & (df>0)] =3
      df_merged = pandas.concat([df_merged,df], axis=0)
outname = i+".selectedvariants.genotype_modified.txt"
#    outname = i+".LTPmutations.genotype_modified.txt"
#    df_merged.index = mutid
df_merged.index = selected_id
df_merged.to_csv(outname, sep="\t", header=False)
ds.close()
