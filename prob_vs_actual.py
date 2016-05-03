import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import pickle

sns.set_context('talk')

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--predict_file',required=True)
parser.add_argument('-d', '--data_file', required=True)
parser.add_argument('-y', '--y_column', required=True)

args = parser.parse_args()
with open(args.predict_file, 'r') as f:
    prediction_df = pd.read_csv(f)
with open(args.data_file, 'rb') as f:
    data_df = pickle.load(f)
prediction_df = prediction_df.sort_values('name')
data_df = data_df.sort_values('name')
y_column = [data_df.loc[data_df['name']==n, args.y_column].values[0] for n in prediction_df['name']]
prediction_df[args.y_column] = y_column
positive = prediction_df[prediction_df.columns[1]] == 1
negative = prediction_df[prediction_df.columns[1]] == -1
plt.plot(prediction_df[positive][args.y_column], prediction_df[positive]['pi'], 'go')
plt.plot(prediction_df[negative][args.y_column], prediction_df[negative]['pi'], 'ro')
plt.margins(0.02)
plt.xlabel(args.y_column)
plt.ylabel('predicted probability')
plt.savefig('plots/' + args.y_column + '_prob_vs_actual.pdf')
print prediction_df

