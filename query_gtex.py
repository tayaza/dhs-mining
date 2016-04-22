#!/usr/bin/python
import argparse
import requests
import json
from sets import Set

tissues = Set(["Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes","Cells_Transformed_fibroblasts","Colon_Sigmoid","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood"])

def getGTExResponse():
	reqList = []
	reqList.append({"snpId":"rs12740374","gencodeId":"ENSG00000134243.7","tissueName":"Liver"})
	reqList.append({"snpId":"rs599839","gencodeId":"ENSG00000134243.7","tissueName":"Liver"})
	gtexResp = requests.post("http://gtexportal.org/api/v6/dyneqtl?v=clversion", json=reqList)
	results = gtexResp.json()["result"]
	print results
	for d in results:
		print d["snpId"]
		print d["geneSymbol"]
		print d["tissueId"]
	#resultsJson = json.loads(results)

if __name__ == "__main__":
	#parser = argparse.ArgumentParser(description="")
	#parser.add_argument("-i","--input")
	#args = parser.parse_args()
	getGTExResponse()