from modeller import *
from modeller.automodel import *

log.none()

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(
    env,
    alnfile='3FDN_A.ali',
    knowns='3FDN_protein',
    sequence='3FDN_A',
    assess_methods=(assess.DOPE, assess.GA341),
)

a.starting_model = 1
a.ending_model = 2

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
