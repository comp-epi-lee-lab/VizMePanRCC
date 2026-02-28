import pandas as pd
from pathlib import Path

pickle_path = Path(r'data/pickle_file.pk1')
table = pd.read_pickle(pickle_path)

i = 0
for column in table.columns[17:]:
    i += 1
print(i)