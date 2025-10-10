<h1>
  <img align="right" alt="Logo" width="120" height="120" src="logo.png">
  <br><br>
  Masato
</h1>

Masato is a modular amplicon sequencing analysis tookit. It highlights transparent and 
composible command line interface that helps users assemble and understand their pipeline.

NOTE: This package is renamed from [easy_amplicon](https://github.com/whatever60/easy_amplicon). The old repository has been archived.

## Try in Colab
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1SPe_950g6cfB1PjDa8zKjgH4xgiiWQL_?usp=sharing)

Thie Colab notebook demonstrates basic usage of the package, example data included.


## Preprocessing

### Primers

- 16SV4 
  ```
  Fwd (2+19bp)
      GTGTGYCAGCMGCCGCGGTAA
      └┘└─────────────────┘
    Linker  16SV4-515F-Y

  Rev (2+20bp)
      CCGGACTACNVGGGTWTCTAAT
      └┘└──────────────────┘
    Linker  16SV4-515R-B
  ```
- 16SV3V4

- 16SV4V5
#### A note about the fesibility of pair merging