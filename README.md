# Ramachandran

Ramachandran plot visualizes energetically allowed regions for backbone dihedral angles ψ against φ of amino acid residues in protein structure.

## Dependencies

* Python 3.5+
* Numpy 1.20.1+
* Matplotlib 3.3.4+

## Usages

It is recommended to use Docker container to run the program.

### Build Docker Image

```
$ docker build -f docker/ramachandran.Dockerfile --no-cache --tag=ramachandran:0.0.1 .
```

### Run Docker Container

```
$ docker run -it --rm --gpus device=0 -v $(pwd):/mnt ramachandran:0.0.1
```

### Install

```
$ pip install ramachandran
```

## TODOs

- [ ] Add PDB crawling tool to compute the favorable and acceptable regions.
- [ ] Add favorable and acceptable regions to the plot.
- [ ] Add mixture of Gaussians for visualizing the clusters of the dihedral angles.
- [ ] Allow the user to show residue name and number in the plot.

## FAQs

### Feature Requests

Please post feature requests on GitHub Issues.

## References

* [Ramachandran Plot](https://en.wikipedia.org/wiki/Ramachandran_plot)
