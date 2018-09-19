import ihm
import os

# If we're running from an SGE job, override the from_pubmed_id() function
# to return a cached value, since we don't have network access (needed to
# query PubMed directly)

def mock_from_pubmed(cls, pubmed_id):
    return ihm.Citation(
            pmid=30190596,
            title='Architecture of Pol II(G) and molecular mechanism of '
                  'transcription regulation by Gdown1.',
            journal='Nat Struct Mol Biol', volume=25, page_range=(859,867),
            year=2018, doi='10.1038/s41594-018-0118-5', authors=[
                'Jishage M', 'Yu X', 'Shi Y', 'Ganesan SJ', 'Chen WY',
                'Sali A', 'Chait BT', 'Asturias FJ', 'Roeder RG'])

if 'JOB_ID' in os.environ:
    ihm.Citation.from_pubmed_id = classmethod(mock_from_pubmed)
