#!/usr/bin/python
import pandas as pd
import glob

adult_age = [
	#month
	'2-month-old stage (mouse)',
	'3-month-old stage (mouse)',
	'4-month-old stage (mouse)',
	'5-month-old stage (mouse)',
	'15-month-old stage (mouse)',
	'20-month-old stage and over (mouse)',
	#week
	'8-week-old stage (mouse)',
	'11-week-old stage (mouse)',
	'14-week-old stage (mouse)',
	'16-week-old stage (mouse)',
	'24-week-old stage (mouse)',
	#other
	'life cycle',
	'post-juvenile',
	'prime adult stage'
	]

def data_processing(dir)
	files = glob.glob(dir + '/*.gz')
	dict = {}
	for file in files:
		df = pd.read_table(file)
		dict[file] = df
		df.head(5)
	matrix = pd.concat(dict)
	matrix.to_csv("Figure.3B-4A_Mus_musculus_RNA-Seq_experiments_libraries_ver.230417.csv")
	return matrix

def data_filtering(matrix):
	partial_matrix = matrix.loc[:,["Experiment ID","Gene ID","Strain","Anatomical entity name","Stage name","TPM"]].fillna("unknown")
	#C57BL/6
	partial_matrix = partial_matrix[partial_matrix["Strain"].str.contains("C57BL/6")]
	#Adult
	partial_matrix = partial_matrix[partial_matrix["Stage name"].isin(adult_age)]
	#Write
	partial_matrix.to_csv("Figure.4A_adult_C57BL6_expr_table.csv")
	"""table_generation"""
	tmp_matrix = partial_matrix.groupby(["Gene ID","Anatomical entity name"]).mean()
	matrix_table = tmp_matrix.reset_index().pivot(values="TPM", index="Gene ID", columns="Anatomical entity name")
	matrix_table.to_csv("Figure.3B-4A_adult_C57BL6_expr_matrix.csv")
	"""list_of_dataset"""
	list = np.unique(partial_matrix["Experiment ID"].tolist())
	with open("Figure.3B-4A_List_of_Dataset.txt", "w") as f:
		_ = f.write(f'{list}')
	f.close()

if __name__ == "__main__":
	print('start processing..')
    dir = 'Experiments_Mus_musculus_RNA-Seq_read_counts_TPM_FPKM'
    """Data Processing"""
	matrix = data_processing(dir = dir)
	"""filtering"""
	data_filtering(matrix)
	print('done.')

