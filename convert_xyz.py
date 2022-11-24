import os
import sys
from argparse import ArgumentParser
import warnings
warnings.filterwarnings("error")

import pandas as pd

from tqdm import tqdm

from openbabel import openbabel as ob
from openbabel import pybel as pb
ob.obErrorLog.SetOutputLevel(ob.obError)

parser = ArgumentParser(description='Convert smiles to xyz')
parser.add_argument('--input', help='input file containing smiles and molecule names')
parser.add_argument('--output', help='output directory')

args = parser.parse_args()

df = pd.read_csv(args.input)

#for i, row in tqdm(df.iterrows(), total=df.shape[0]):
message_counter = 0
for i, row in df.iterrows():
    mol = pb.readstring('smi', row['smiles'])
    mol.addh()
    mol.make3D()
    mol.write('xyz', os.path.join(args.output, str(row['num']) + '.xyz'), overwrite=True)
    if ob.obErrorLog.GetWarningMessageCount() > message_counter:
        print('\nProcessing molecule {}. Molecule name and SMILES: {} {}'.format(row['num'], row['name'], row['smiles']))
        print(ob.obErrorLog.GetMessagesOfLevel(ob.obWarning)[-1])
        ob.obErrorLog.ClearLog()
        message_counter = ob.obErrorLog.GetWarningMessageCount()