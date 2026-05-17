from modeller import *
from modeller.automodel import *

log.verbose()

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(
    env,
    alnfile='chain_alignment.ali',
    knowns='3BC3_protein_monomer',
    sequence='3BC3_A',
    assess_methods=(assess.DOPE, assess.GA341),
)

a.starting_model = 1
a.ending_model = 1

a.make()

score_file = 'model_scores.tsv'
with open(score_file, 'w', encoding='utf-8') as handle:
    handle.write("model_name\tdope_score\tga341_score\n")
    for model in a.outputs:
        if model.get('failure') is None:
            model_name = model.get('name')
            dope_score = model.get('DOPE score')
            ga341_score = model.get('GA341 score')
            handle.write(f"{model_name}\t{dope_score}\t{ga341_score}\n")
