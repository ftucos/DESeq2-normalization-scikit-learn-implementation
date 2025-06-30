# DESeq2-sklearn (playground)

*A small, exploratory repo that ports the DESeq2 size-factor normalisation from R to a scikit-learn–style transformer.*

------

## What this repo contains

| Folder / file                                    | Purpose                                                      |
| ------------------------------------------------ | ------------------------------------------------------------ |
| **`01-simulate_data_and_compute_sizefactors.R`** | Simulates a toy RNA-seq dataset (20 samples × 2 000 genes) and computes size factors with a **patched** `DESeq2::estimateSizeFactorsForMatrix()` (geometric-mean rescaling **disabled**). |
| **`02-test_python_implementation.ipynb`**        | Pure-Python `DESeq2` transformer (`fit`, `transform`, `get_size_factors`) that plugs straight into scikit-learn pipelines. Built on top of the RNAnorm code base. |

Use in Python:

```
from src.deseq2_norm import DESeq2

deseq = DESeq2().fit(X_train)
train_norm = deseq.transform(X_train)
test_norm  = deseq.transform(X_test)   # uses training *reference pseudosample*
```

------

## Important caveats ⚠️

- **No geometric-mean rescaling** – size factors remain on the raw training library-size scale.
- **Sparse matrices only partially tested** – dense `DataFrame` are the primary target.
- **Still very much a work-in-progress** – things will move around and break; kick the tyres on your own data before relying on it.