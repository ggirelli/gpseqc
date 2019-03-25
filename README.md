gpseqc v2.3.5
===

A Python3 package that provides tools to estimate the 3D spatial nuclear centrality of genomic regions and compare different centrality rankings.

* Read the (public) [documentation](https://ggirelli.github.io/gpseqc/) for more details.
* Read the (private) [documentation](https://github.com/ggirelli/gpseqc/wiki) for more details.  
*Once the repo goes public, private docs will be merged with public ones.*

Installation
---

To **install**, run the following:

```
git clone http://github.com/ggirelli/gpseqc
cd gpseqc
sudo -H pip3 install .
```

To **test** your installation, run: `pytest-3 --pyargs gpseqc`.

To **uninstall** run the following from within the repository folder:

```
sudo -H pip3 uninstall gpseqc
```

To **update**, first uninstall, and then run the following from within the repository folder.

```
git pull
sudo -H pip3 install .
```

### Additional dependencies

Most of the required dependencies are automatically install by `pip3`. Still, some require some manual steps. Specifically, `gpseqc` requires the following packages:

* `tkinter`

That on Ubuntu can be easily installed with:

```bash
sudo apt install python3-tk
```

Usage
---

#### Estimate centrality

The `gpseqc_estimate` script allows to estimate regional nuclear centrality based on a multi-condition GPSeq experiment. Run `gpseqc_estimate -h` for more details.

#### Compare centrality ranks

The `gpseqc_compare` script allows to compare different regional centrality ranks. Run `gpseqc_compare -h` for more details.

Contributing
---

We welcome any contributions to `GPSeqC`. Please, refer to the [contribution guidelines](https://ggirelli.github.io/gpseqc/contributing) if this is your first time contributing! Also, check out our [code of conduct](https://ggirelli.github.io/gpseqc/code_of_conduct).

License
---

```
MIT License
Copyright (c) 2017-18 Gabriele Girelli
```
